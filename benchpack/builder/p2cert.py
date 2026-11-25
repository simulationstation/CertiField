"""P2 certificates (L-bracket): find the fields (float), evaluate them exactly, store, re-verify.

usage:
    python p2cert.py build  STEM [preset]      build and store certificates/STEM.{json,npz}
    python p2cert.py verify STEM               re-evaluate a stored certificate with this builder's evaluator
All certified numbers of P2 are exact rationals (no transcendental constants occur).
"""
import limits  # noqa: F401  (thread and process limits; must come before numpy)
import json
import math
import os
import shutil
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import p2float as pf
from p2float import Setup2, graded_breaks, solve_spd
import p2rig
import certio as cio
from splines import Spline1D

CST = pf.C1111 / pf.PATCH_AREA / 2                     # (1/2) C_1111 / |S_p|


def params(h=Fr(1, 8), L=10, m=2, ps=6, pu=5, Lc=None, Lp=4, umul0=1):
    Lc = L if Lc is None else Lc
    t = [(0, (-1, 1), L, m), (-1, (1,), Lc, m), (1, (-1,), Lc, m), (Fr(1, 10), (-1, 1), Lp, m)]
    br = graded_breaks(-1, 1, h, t, extra=[Fr(1, 10)])
    return dict(bx=br, by=br, ps=ps, pu=pu,
                smult={'x': {0: ps - 1, Fr(1, 10): ps - 1}, 'y': {0: ps - 1, Fr(1, 10): ps - 1}},
                umult={'x': {0: umul0, Fr(1, 10): pu}, 'y': {0: umul0, Fr(1, 10): pu}})


def setup_from_params(par):
    return Setup2(par['bx'], par['by'], par['ps'], par['bx'], par['by'], par['pu'], smult=par['smult'], umult=par['umult'])


def setup_from_meta(meta):
    ds, ss = meta['displacement_space'], meta['stress_space']
    S = Setup2.__new__(Setup2)
    S.sx = Spline1D([Fr(k) for k in ss['knots_x']], ss['degree'])
    S.sy = Spline1D([Fr(k) for k in ss['knots_y']], ss['degree'])
    S.ux = Spline1D([Fr(k) for k in ds['knots_x']], ds['degree'])
    S.uy = Spline1D([Fr(k) for k in ds['knots_y']], ds['degree'])
    S.ps, S.pu = ss['degree'], ds['degree']
    return S


def check_ties(S, fields):
    free = S.stress_free()
    for name, f in fields.items():
        if 'A' in f:
            for M in (f['A'], f['B']):
                assert np.all(M[0, :] == 0) and np.all(M[-1, :] == 0), name
                assert np.all(np.isfinite(M)), name
        else:
            assert np.all(f['C'][~free] == 0), name
            assert np.all(np.isfinite(f['C'])) and np.all(np.isfinite(f['par'])), name
    for Sx in (S.sx, S.sy):
        assert sum(1 for k in Sx.t if k == 0) == S.ps - 1


def bounds(out, kappa):
    """exact rational bounds from the exact functionals"""
    res = dict(K_up=out['v0']['a'], K_low=out['tau0']['F'] ** 2 / out['tau0']['Q'])
    if kappa is not None:
        k = Fr(kappa)
        vp, vm, tp, tm = out['vp'], out['vm'], out['tp'], out['tm']
        res['m_up_plus'] = vp['a'] / 2 - k * vp['l']
        res['m_up_minus'] = vm['a'] / 2 + k * vm['l']
        res['m_low_plus'] = tp['F'] - tp['Q'] / 2 - k * tp['l'] - k * k * CST
        res['m_low_minus'] = tm['F'] - tm['Q'] / 2 + k * tm['l'] - k * k * CST
        res['Q2_low'] = (res['m_low_minus'] - res['m_up_plus']) / (2 * k)
        res['Q2_up'] = (res['m_up_minus'] - res['m_low_plus']) / (2 * k)
    return res


def solve(S, verbose=True):
    t0 = time.time()
    A, f, lt = S.assemble_dual()
    Pr = S.stress_ties()
    Z = solve_spd(A, np.stack([f, -lt], 1), Pr, protect=3)
    K, b, k0, l0, lv = S.assemble_primal()
    Pv = S.disp_ties()
    V = solve_spd(K, np.stack([-b, lv], 1), Pv)
    cst = float(CST)

    def mup(k):
        v = V[:, 0] + k * V[:, 1]
        return 0.5 * (v @ (K @ v) + 2 * b @ v + k0) - k * (l0 + lv @ v)

    def mlow(k):
        z = Z[:, 0] + k * Z[:, 1]
        return f @ z - 0.5 * z @ (A @ z) - k * lt @ z - k * k * cst

    U0, L0 = mup(0.0), mlow(0.0)
    U2 = (mup(1e-2) + mup(-1e-2) - 2 * U0) / 2e-4
    L2 = (mlow(1e-2) + mlow(-1e-2) - 2 * L0) / 2e-4
    kf = math.sqrt(max(U0 - L0, 1e-16) / (U2 - L2))
    if verbose:
        print(f"  float: K in [{2 * L0:.10f}, {2 * U0:.10f}], adjoint gap {U2 - L2:.3e}, kappa {kf:.3e} "
              f"({time.time() - t0:.0f}s)", flush=True)
    return V, Z, kf


def _fingerprint(par):
    return json.dumps({k: (str(v) if not isinstance(v, list) else [str(x) for x in v]) for k, v in sorted(par.items())})


def build(stem, par, verbose=True):
    """resumable: rerun the same command; the float solution and the finished chunks of the exact evaluation are
    kept in benchpack/work/STEM and reused"""
    t0 = time.time()
    S = setup_from_params(par)
    wd = os.path.join(limits.WORK_DIR, stem) if stem else None
    if verbose:
        print(f"  {S.sx.ne} x {S.sy.ne} cells, degrees {S.ps}/{S.pu}"
              + (f"; to resume after a stop rerun the same command (checkpoints in benchpack/work/{stem})" if wd else ''),
              flush=True)
    sol = os.path.join(wd, 'solve.npz') if wd else None
    fp = _fingerprint(par)
    if sol and os.path.exists(sol) and str(np.load(sol)['fp']) == fp:
        d = np.load(sol)
        V, Z, kf = d['V'], d['Z'], float(d['kf'])
        if verbose:
            print("  float solution loaded from the checkpoint", flush=True)
    else:
        V, Z, kf = solve(S, verbose)
        if sol:
            if os.path.isdir(wd):                       # stale checkpoints of another solution
                shutil.rmtree(wd)
            os.makedirs(wd, exist_ok=True)
            np.savez(sol + '.tmp.npz', V=V, Z=Z, kf=kf, fp=fp)
            os.replace(sol + '.tmp.npz', sol)
    kappa = cio.kappa_dyadic(kf)
    kk = float(kappa)
    nA = S.ux.n * S.uy.n
    fields = {}
    for name, v in (('v0', V[:, 0]), ('vp', V[:, 0] + kk * V[:, 1]), ('vm', V[:, 0] - kk * V[:, 1])):
        fields[name] = dict(A=v[:nA].reshape(S.ux.n, S.uy.n).copy(), B=v[nA:].reshape(S.ux.n, S.uy.n).copy())
    for name, z in (('tau0', Z[:, 0]), ('tp', Z[:, 0] + kk * Z[:, 1]), ('tm', Z[:, 0] - kk * Z[:, 1])):
        fields[name] = dict(par=z[:3].copy(), C=z[3:].reshape(S.sx.n, S.sy.n).copy())
    check_ties(S, fields)
    out = p2rig.evaluate(S, fields, checkpoint=os.path.join(wd, 'chunks') if wd else None)
    res = bounds(out, kappa)
    K_lo, K_hi, Q_lo, Q_hi = res['K_low'], res['K_up'], res['Q2_low'], res['Q2_up']
    assert K_lo <= K_hi and Q_lo <= Q_hi
    cK, cQ = cio.claim_strings(K_lo, K_hi), cio.claim_strings(Q_lo, Q_hi)
    meta = dict(
        format='benchpack-certificate-1', problem='P2', variant='standard', model='full',
        displacement_space=dict(degree=S.pu, knots_x=[cio.fstr(t) for t in S.ux.t], knots_y=[cio.fstr(t) for t in S.uy.t]),
        stress_space=dict(degree=S.ps, knots_x=[cio.fstr(t) for t in S.sx.t], knots_y=[cio.fstr(t) for t in S.sy.t]),
        fields=dict(Q1=dict(v='v0', tau='tau0'),
                    Q2=dict(kappa=cio.fstr(kappa), v_plus='vp', v_minus='vm', tau_plus='tp', tau_minus='tm')),
        claims=dict(Q1=dict(lower=cK[0], upper=cK[1]), Q2=dict(lower=cQ[0], upper=cQ[1])),
        computed=dict(
            fields={n: {k: float(v) for k, v in d.items()} for n, d in out.items()},
            bounds={k: float(v) for k, v in res.items()},
            K=[float(K_lo), float(K_hi)], Q2=[float(Q_lo), float(Q_hi)],
            note='all values above are exact rationals rounded to the nearest float'),
        label=cio.LABEL,
    )
    arrays = {}
    for name, fdat in fields.items():
        if 'A' in fdat:
            arrays[name + '_A'], arrays[name + '_B'] = fdat['A'], fdat['B']
        else:
            arrays[name + '_C'], arrays[name + '_par'] = fdat['C'], fdat['par']
    if stem is not None:
        cio.save(stem, meta, arrays)
    meta['_fields'], meta['_res'] = fields, res
    if verbose:
        rw = lambda lo, hi: float((hi - lo) / abs(hi))
        print(f"  CERTIFIED  {float(K_lo):.12f} <= K  <= {float(K_hi):.12f}   rel width {rw(K_lo, K_hi):.2e}")
        print(f"  CERTIFIED  {float(Q_lo):.10f} <= Q2 <= {float(Q_hi):.10f}   rel width {rw(Q_lo, Q_hi):.2e}")
        print(f"  claims: K in [{cK[0]}, {cK[1]}],  Q2 in [{cQ[0]}, {cQ[1]}]   ({time.time() - t0:.0f}s)", flush=True)
    return meta


def verify(stem, verbose=True):
    meta, arr = cio.load(stem)
    S = setup_from_meta(meta)
    fields = {}
    for name in ('v0', 'vp', 'vm'):
        if name + '_A' in arr:
            fields[name] = dict(A=arr[name + '_A'], B=arr[name + '_B'])
    for name in ('tau0', 'tp', 'tm'):
        if name + '_C' in arr:
            fields[name] = dict(C=arr[name + '_C'], par=arr[name + '_par'])
    check_ties(S, fields)
    kappa = Fr(meta['fields']['Q2']['kappa']) if 'Q2' in meta['fields'] else None
    res = bounds(p2rig.evaluate(S, fields), kappa)
    c = meta['claims']
    ok = Fr(c['Q1']['lower']) <= res['K_low'] and res['K_up'] <= Fr(c['Q1']['upper'])
    if kappa is not None:
        ok = ok and Fr(c['Q2']['lower']) <= res['Q2_low'] and res['Q2_up'] <= Fr(c['Q2']['upper'])
    if verbose:
        print(f"{stem}: K in [{float(res['K_low']):.12f}, {float(res['K_up']):.12f}]"
              + (f", Q2 in [{float(res['Q2_low']):.10f}, {float(res['Q2_up']):.10f}]" if kappa else '')
              + f"; claims {c}: {'OK' if ok else 'FAILED'}")
    return ok


PRESETS = {
    'small': dict(h=Fr(1, 4), L=3, m=1, ps=4, pu=3, Lc=3, Lp=2),
}

if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    if sys.argv[1] == 'build':
        stem = sys.argv[2]
        kw = dict(PRESETS[sys.argv[3]]) if len(sys.argv) > 3 and sys.argv[3] in PRESETS else {}
        for a in sys.argv[3:]:
            if '=' in a:
                k, v = a.split('=')
                kw[k] = Fr(v) if k == 'h' else int(v)
        build(stem, params(**kw))
    elif sys.argv[1] == 'verify':
        sys.exit(0 if verify(sys.argv[2]) else 1)
