"""P1 / P1c certificates: find the fields (float), evaluate them rigorously, store, re-verify.

usage:
    python p1cert.py build  STEM N ps pu [frictionless|clamped]     build and store certificates/STEM.{json,npz}
    python p1cert.py verify STEM                                    re-evaluate a stored certificate (this builder's
                                                                    evaluator; the independent checker is separate)
"""
import limits  # noqa: F401  (thread and process limits; must come before numpy)
import math
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from flint import arb

import p1float as pf
from p1float import Setup, solve_spd, uniform_breaks
from p1rig import Rig
from exact import q
import certio as cio

KHALF = pf.D.P1_QUARTER['K_from_quarter_energy']          # K = 4 a_q(u,u)


def setup_from_params(par):
    return Setup(par['sbx'], par['sby'], par['ps'], par['ubx'], par['uby'], par['pu'], gexp=par.get('gexp'),
                 gbase=par.get('gbase', Fr(1, 2)), grip=par.get('grip', 'frictionless'), hole=par.get('hole', True),
                 smult=par.get('smult'), umult=par.get('umult'))


def setup_from_meta(meta):
    ds, ss = meta['displacement_space'], meta['stress_space']
    S = Setup.__new__(Setup)
    from splines import Spline1D
    S.sx = Spline1D([Fr(k) for k in ss['knots_x']], ss['degree'])
    S.sy = Spline1D([Fr(k) for k in ss['knots_y']], ss['degree'])
    S.ux = Spline1D([Fr(k) for k in ds['knots_x']], ds['degree'])
    S.uy = Spline1D([Fr(k) for k in ds['knots_y']], ds['degree'])
    S.ps, S.pu = ss['degree'], ds['degree']
    S.gexp = tuple(ss['g_exponents'])
    S.gbase = Fr(ss['g_base'])
    S.hole = meta.get('hole', True)
    S.grip = meta['variant']
    S.nqs, S.nqu = S.ps + 9, S.pu + 1
    S.theta0 = math.pi / 2 - 0.1
    return S


def split_stress(S, z):
    ng = len(S.gexp)
    return z[:ng].copy(), z[ng:].reshape(S.sx.n, S.sy.n).copy()


def split_disp(S, v):
    nA = S.ux.n * S.uy.n
    return v[:nA].reshape(S.ux.n, S.uy.n).copy(), v[nA:].reshape(S.ux.n, S.uy.n).copy()


def check_ties(S, fields):
    """exact checks of the coefficient conditions of FORMAT.md 3.1 / 3.2"""
    for name, f in fields.items():
        if 'A' in f:
            A, B = f['A'], f['B']
            assert np.all(A[0, :] == 0) and np.all(A[-1, :] == 0) and np.all(B[:, 0] == 0), name
            if S.grip == 'clamped':
                assert np.all(B[-1, :] == 0), name
        else:
            C = f['C']
            assert np.array_equal(C[0, :], C[1, :]), name
            if S.grip == 'frictionless':
                assert np.array_equal(C[-1, :], C[-2, :]), name
            assert np.array_equal(C[:, 0], C[:, 1]), name
            assert np.all(C[:, -1] == 0) and np.all(C[:, -2] == 0), name
            assert np.all(np.isfinite(C)) and np.all(np.isfinite(f['par'])), name


def evaluate(S, fields, kappa, verbose=True):
    """rigorous balls of every functional and of the bounds"""
    R = Rig(S)
    out = {}
    t0 = time.time()
    for name, f in fields.items():
        if 'A' in f:
            out[name] = dict(a=R.a(f['A'], f['B']), l=R.l_v(f['A'], f['B']))
        else:
            F = R.F(f['par'])
            out[name] = dict(Q=R.Q(f['C'], f['par']), F=arb(q(F)), Fexact=F, l=R.l_tau(f['C'], f['par']))
        if verbose:
            print(f"    {name}: evaluated ({time.time() - t0:.0f}s)", flush=True)
    res = {}
    v0, t0_ = out['v0'], out['tau0']
    res['K_up'] = arb(q(KHALF)) * v0['a']
    res['K_low'] = t0_['F'] * t0_['F'] / t0_['Q']
    if kappa is not None:
        k = arb(q(kappa))
        cst = arb(q(pf.C1111 / pf.PH_AREA)) / 2
        vp, vm, tp, tm = out['vp'], out['vm'], out['tp'], out['tm']
        res['m_up_plus'] = vp['a'] / 2 - k * vp['l']
        res['m_up_minus'] = vm['a'] / 2 + k * vm['l']
        res['m_low_plus'] = tp['F'] / 2 - tp['Q'] / 2 - k * tp['l'] - k * k * cst
        res['m_low_minus'] = tm['F'] / 2 - tm['Q'] / 2 + k * tm['l'] - k * k * cst
        res['Q2_low'] = (res['m_low_minus'] - res['m_up_plus']) / (2 * k)
        res['Q2_up'] = (res['m_up_minus'] - res['m_low_plus']) / (2 * k)
    return out, res


def build(stem, par, verbose=True):
    t0 = time.time()
    S = setup_from_params(par)
    A, f = S.assemble_dual()
    Pr = S.stress_ties()
    lt = S.patch_dual()
    ng = len(S.gexp)
    Z = solve_spd(A, np.stack([0.5 * f, -lt], 1), Pr, protect=ng)
    K, b, k0 = S.assemble_primal()
    Pv = S.disp_ties()
    l0, lv = S.patch_primal()
    V = solve_spd(K, np.stack([-b, lv], 1), Pv)
    cst = float(pf.C1111 / pf.PH_AREA)

    def mup(k):
        v = V[:, 0] + k * V[:, 1]
        return 0.5 * (v @ (K @ v) + 2 * b @ v + k0) - k * (l0 + lv @ v)

    def mlow(k):
        z = Z[:, 0] + k * Z[:, 1]
        return 0.5 * f @ z - 0.5 * z @ (A @ z) - k * lt @ z - 0.5 * k * k * cst

    U0, L0 = mup(0.0), mlow(0.0)
    U2 = (mup(1e-2) + mup(-1e-2) - 2 * U0) / 2e-4
    L2 = (mlow(1e-2) + mlow(-1e-2) - 2 * L0) / 2e-4
    kf = math.sqrt(max(U0 - L0, 1e-16) / (U2 - L2))
    kappa = cio.kappa_dyadic(kf)
    kk = float(kappa)
    if verbose:
        print(f"  float: K in [{8 * L0:.12f}, {8 * U0:.12f}], kappa {kf:.3e} -> {kappa} ({time.time() - t0:.0f}s)", flush=True)
    fields = {}
    for name, v in (('v0', V[:, 0]), ('vp', V[:, 0] + kk * V[:, 1]), ('vm', V[:, 0] - kk * V[:, 1])):
        Am, Bm = split_disp(S, v)
        fields[name] = dict(A=Am, B=Bm)
    for name, z in (('tau0', Z[:, 0]), ('tp', Z[:, 0] + kk * Z[:, 1]), ('tm', Z[:, 0] - kk * Z[:, 1])):
        g, C = split_stress(S, z)
        fields[name] = dict(par=g, C=C)
    check_ties(S, fields)
    out, res = evaluate(S, fields, kappa, verbose)
    K_lo, K_hi = cio.ball_lo(res['K_low']), cio.ball_hi(res['K_up'])
    Q2_lo, Q2_hi = cio.ball_lo(res['Q2_low']), cio.ball_hi(res['Q2_up'])
    assert K_lo <= K_hi and Q2_lo <= Q2_hi
    cK, cQ = cio.claim_strings(K_lo, K_hi), cio.claim_strings(Q2_lo, Q2_hi)
    meta = dict(
        format='benchpack-certificate-1', problem='P1', variant=S.grip, model='quarter',
        displacement_space=dict(degree=S.pu, knots_x=[cio.fstr(t) for t in S.ux.t], knots_y=[cio.fstr(t) for t in S.uy.t]),
        stress_space=dict(degree=S.ps, knots_x=[cio.fstr(t) for t in S.sx.t], knots_y=[cio.fstr(t) for t in S.sy.t],
                          g_exponents=list(S.gexp), g_base=cio.fstr(S.gbase)),
        fields=dict(Q1=dict(v='v0', tau='tau0'),
                    Q2=dict(kappa=cio.fstr(kappa), v_plus='vp', v_minus='vm', tau_plus='tp', tau_minus='tm')),
        claims=dict(Q1=dict(lower=cK[0], upper=cK[1]), Q2=dict(lower=cQ[0], upper=cQ[1])),
        computed=dict(
            fields={n: {k: (cio.encl(v) if k != 'Fexact' else cio.fstr(v)) for k, v in d.items()} for n, d in out.items()},
            bounds={k: cio.encl(v) for k, v in res.items()},
            K=[float(K_lo), float(K_hi)], Q2=[float(Q2_lo), float(Q2_hi)]),
        label=cio.LABEL,
    )
    arrays = {}
    for name, fdat in fields.items():
        if 'A' in fdat:
            arrays[name + '_A'], arrays[name + '_B'] = fdat['A'], fdat['B']
        else:
            arrays[name + '_C'], arrays[name + '_par'] = fdat['C'], fdat['par']
    if not S.hole:
        meta['hole'] = False                       # test certificates only (plate without hole)
    if stem is not None:
        cio.save(stem, meta, arrays)
    meta['_fields'], meta['_res'] = fields, res
    if verbose:
        rw = lambda lo, hi: float((hi - lo) / abs(hi))
        print(f"  CERTIFIED  {float(K_lo):.12f} <= K  <= {float(K_hi):.12f}   rel width {rw(K_lo, K_hi):.2e}")
        print(f"  CERTIFIED  {float(Q2_lo):.10f} <= Q2 <= {float(Q2_hi):.10f}   rel width {rw(Q2_lo, Q2_hi):.2e}")
        print(f"  claims: K in [{cK[0]}, {cK[1]}],  Q2 in [{cQ[0]}, {cQ[1]}]   ({time.time() - t0:.0f}s)")
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
    out, res = evaluate(S, fields, kappa, verbose=False)
    ok = Fr(meta['claims']['Q1']['lower']) <= cio.ball_lo(res['K_low']) and cio.ball_hi(res['K_up']) <= Fr(meta['claims']['Q1']['upper'])
    if kappa is not None:
        ok = ok and Fr(meta['claims']['Q2']['lower']) <= cio.ball_lo(res['Q2_low']) and \
            cio.ball_hi(res['Q2_up']) <= Fr(meta['claims']['Q2']['upper'])
    if verbose:
        print(f"{stem}: K in [{float(cio.ball_lo(res['K_low'])):.12f}, {float(cio.ball_hi(res['K_up'])):.12f}]"
              + (f", Q2 in [{float(cio.ball_lo(res['Q2_low'])):.10f}, {float(cio.ball_hi(res['Q2_up'])):.10f}]" if kappa else '')
              + f"; claims {meta['claims']}: {'OK' if ok else 'FAILED'}")
    return ok


PRESETS = {
    # name: (N, ps, pu)
}


def params(N, ps, pu, grip='frictionless', Lg=0, m=1):
    """uniform N x N cells; Lg > 0 adds geometric layers (m per halving) towards x = 1 and y = 1 (the corner where a
    clamped grip meets the free edge, P1c)"""
    br = uniform_breaks(N)
    if Lg:
        from p2float import graded_breaks
        br = graded_breaks(0, 1, Fr(1, N), [(1, (-1,), Lg, m)])
    return dict(sbx=br, sby=br, ps=ps, ubx=br, uby=br, pu=pu, grip=grip)


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    if sys.argv[1] == 'build':
        stem, N, ps, pu = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        grip = sys.argv[6] if len(sys.argv) > 6 else 'frictionless'
        Lg = int(sys.argv[7]) if len(sys.argv) > 7 else 0
        m = int(sys.argv[8]) if len(sys.argv) > 8 else 1
        build(stem, params(N, ps, pu, grip, Lg, m))
    elif sys.argv[1] == 'verify':
        sys.exit(0 if verify(sys.argv[2]) else 1)
