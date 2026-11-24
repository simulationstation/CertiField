#!/usr/bin/env python3
"""Grids for the certified level-set programme in degree B (Kernel/SkyrmeLevelSetDegree.wl).

    python make_grids.py B MU PRESET [OUTDIR]        e.g.   python make_grids.py 2 0 quick

writes  OUTDIR/b<B>_mu<MU>_<preset>/  (default OUTDIR: ../../Data/SkyrmeDegree) with
    taus.f64   the levels tau_0 < ... < tau_N
    vref.f64   the volume of the extremal of the relaxed problem at each level
    ug.f64     the logarithmic offsets u of the volume nodes,  eug.f64 = exp(ug)
    offs.i64   index of the lowest node kept at each level
    ctg.f64    cost of the extremal beyond each level (used only by the path test)
    meta.json  parameters
as raw little-endian IEEE doubles / 64-bit integers.

The data are heuristic.  The nodes of level i are 0 and vref[i] * eug[offs[i]:].  They place the cells on which the
certified programme builds its piecewise affine lower bounds; any positive, increasing data give a valid bound.
The layout follows the description of the grids in the level-set note (levels with steps min(dtmax, cg / lambda),
lambda = d ln vref / d tau; nodes at fixed offsets u = ln(v / vref), spacing du for -U1 <= u <= U2, growing
geometrically outside; lowest node at radius eps1 sin tau).
"""
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

import relaxed as R

PRESETS = {
    'quick': dict(du=0.04, cg=0.04),
    'standard': dict(du=0.02, cg=0.02),
    'fine': dict(du=0.01, cg=0.01),
}
COMMON = dict(U1=4.0, U2=2.0, growth=1.05, dumax=0.2, ulo=-32.0, uhi=22.0, t1=0.03, lN=0.02, eps1=0.1, dtmax=0.004,
              m1=8.0, m2=3.0, nmin=6, table=[2000, 900, 1e-10])
MU = {'0': 0.0, '1': 1.0}


def u_grid(du, U1, U2, growth, dumax, ulo, uhi):
    up, e = [0.0], du
    while up[-1] < uhi:
        up.append(up[-1] + e)
        if up[-1] >= U2:
            e = min(e * growth, dumax)
    dn, e = [0.0], du
    while dn[-1] > ulo:
        dn.append(dn[-1] - e)
        if dn[-1] <= -U1:
            e = min(e * growth, dumax)
    return np.array(dn[:0:-1] + up)


def reference_path(res):
    """ln r(tau) of the extremal and lambda = d ln v / d tau, interpolated in s = ln(tau / (pi - tau))"""
    t, r = res['t'], res['r']
    keep = np.concatenate([[True], np.diff(t) > 0])
    t, r = t[keep], r[keep]
    sv = np.log(t / (math.pi - t))
    lr = PchipInterpolator(sv, np.log(r))
    s0, s1 = sv[0], sv[-1]
    dlr = lr.derivative()

    def lam(tau):
        s = math.log(tau / (math.pi - tau))
        return 3.0 * float(dlr(min(max(s, s0), s1))) * (1.0 / tau + 1.0 / (math.pi - tau))

    return lam


def extremal_at(res, taus):
    """radius and accumulated cost of the extremal at the given levels (integrated between consecutive levels)"""
    B, mu, par = res['B'], res['mu'], res['par']
    first = R.shoot(B, par, mu, eps=math.pi - taus[0])              # from tau0 to the first level
    y = [float(first['r'][-1]), float(first['kap'][-1]), float(first['cost'][-1])]
    assert abs(first['t'][-1] - taus[0]) < 1e-12 and not first['early']

    def rhs(tau, y):
        r, kap, _ = y
        r = max(r, 1e-12)
        V = R.Vpot(mu, tau)
        kap = max(kap, V, 1e-300)
        c = B * R.FOURPI * math.sin(tau) ** 2
        D = c / (R.FOURPI * r * r)
        hs, hk, I, lam = R.envelope(kap, D)
        if I < 0:
            ds = R.D_star(kap)
            vp = c * R.h_k(kap, ds) / ds
            return [vp / (R.FOURPI * r * r), -2 * mu * mu * math.sin(tau), c * R.h(kap, ds) / ds - (kap - V) * vp]
        return [hk, -2 * mu * mu * math.sin(tau) - 2.0 / r * I, R.FOURPI * r * r * (hs - (kap - V) * hk)]

    rr, cc = [y[0]], [y[2]]
    for a, b in zip(taus[:-1], taus[1:]):
        sol = solve_ivp(rhs, (float(a), float(b)), y, method='LSODA', rtol=1e-10, atol=1e-13)
        assert sol.status == 0
        y = [float(x) for x in sol.y[:, -1]]
        rr.append(y[0])
        cc.append(y[2])
    return np.array(rr), np.array(cc)


def make(B, mu_name, preset, outdir):
    mu = MU[mu_name]
    p = dict(COMMON)
    p.update(PRESETS[preset])
    res = R.solve(B, mu)
    assert res is not None and not res['early'], "no extremal found"
    lam = reference_path(res)
    tN = math.pi - p['lN']
    taus = [p['t1']]
    while taus[-1] < tN - 1e-12:
        dt = min(p['dtmax'], p['cg'] / max(lam(taus[-1]), 1e-6))
        taus.append(min(taus[-1] + dt, tN))
    taus = np.array(taus)
    N = len(taus) - 1
    rref, cost = extremal_at(res, taus)
    vref = R.FOURPI / 3 * rref ** 3
    assert np.all(np.diff(vref) > 0) and np.all(np.diff(taus) > 0)
    ctg = res['total'] - cost
    ug = u_grid(p['du'], p['U1'], p['U2'], p['growth'], p['dumax'], p['ulo'], p['uhi'])
    i0 = int(np.argmin(np.abs(ug)))
    assert ug[i0] == 0.0
    ulow = 3 * np.log(p['eps1'] * np.sin(taus) / rref)
    offs = np.clip(np.searchsorted(ug, ulow, side='left'), 0, i0 - 1)
    offs = np.maximum.accumulate(offs[::-1])[::-1].astype('<i8')
    d = os.path.join(outdir, f"b{B}_mu{mu_name}_{preset}")
    os.makedirs(d, exist_ok=True)
    for name, arr in (('taus', taus), ('vref', vref), ('ug', ug), ('eug', np.exp(ug)), ('ctg', ctg)):
        np.asarray(arr, dtype='<f8').tofile(os.path.join(d, name + '.f64'))
    offs.tofile(os.path.join(d, 'offs.i64'))
    meta = {
        "model": "SkyrmeDegreeTwo" if B == 2 else f"Skyrme, degree {B}",
        "certificate": "level-set lower bound in degree B (docs/SKYRME_DEGREE_TWO.md, Theorem 2)",
        "degree": B, "mu": mu, "mu_exact": mu_name, "preset": preset, "params": p,
        "N": N, "n_u": int(len(ug)), "i0": i0,
        "files": {"taus": "taus.f64", "vref": "vref.f64", "eug": "eug.f64", "ug": "ug.f64", "offs": "offs.i64", "ctg": "ctg.f64"},
        "note": "Grid data are heuristic; they affect the tightness of the bound, not its validity.",
        "extremal": {"parameter": res['par'], "value_float": res['value_units'],
                     "note": "floating-point value of the relaxed problem along its extremal, in units of 12 pi^2; not certified"},
    }
    with open(os.path.join(d, 'meta.json'), 'w', newline='\n') as f:
        json.dump(meta, f, indent=1)
        f.write('\n')
    print(f"{d}: N = {N} levels, {len(ug)} offsets, extremal value {res['value_units']:.7f}", flush=True)


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    here = os.path.dirname(os.path.abspath(__file__))
    out = sys.argv[4] if len(sys.argv) > 4 else os.path.join(here, '..', '..', 'Data', 'SkyrmeDegree')
    make(int(sys.argv[1]), sys.argv[2], sys.argv[3], os.path.normpath(out))
