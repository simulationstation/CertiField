"""Tests of the P2 builder (python tests_p2.py [STEM]).

1. Known answers of the exact evaluator: the lift alone, v = (0, -(x+1)/2), has a(v,v) = (5/13)(1/4)|Omega| = 15/52
   and l(v) = 0; the particular stress alone (Psi = 0) has Q, F and l given in closed form by sympy (an independent
   code path), on the L-shaped domain.
2. Exact identities of random admissible stress fields: equilibrium, zero traction on the four free edges, traction
   continuity across cell edges, F = -beta (the Airy part of Psi carries no vertical grip force).
3. The exact evaluator agrees with the float assembly on random admissible fields.
4. Falsification against the stored certificate: perturbed certified fields and other kappa never cross.
5. The stored certificate re-evaluates to its claims.
"""
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp_

import p2cert as pc
import p2float as pf
import p2rig
import certio as cio
from exact import q, X, Y, rect_integral
from p2rig import _Cache, _spline_cell, CHI

rng = np.random.default_rng(11)
T0 = time.time()
STEM = sys.argv[1] if len(sys.argv) > 1 else None


def head(s):
    print(f"[{time.time() - T0:5.0f}s] {s}", flush=True)




def main():
    global rng
    S = pc.setup_from_params(pc.params(h=Fr(1, 4), L=3, m=1, ps=4, pu=3, Lc=2, Lp=2))

    # ---------------------------------------------------------------------------------------------------------------------
    head("1. known answers of the exact evaluator")
    zA = np.zeros((S.ux.n, S.uy.n))
    al, be, ga = 0.375, -1.25, 0.5
    fields = dict(v0=dict(A=zA, B=zA), tau0=dict(C=np.zeros((S.sx.n, S.sy.n)), par=np.array([al, be, ga])))
    out = p2rig.evaluate(S, fields, nproc=2)
    assert out['v0']['a'] == Fr(15, 52) and out['v0']['l'] == 0, out['v0']
    x, y = sp_.symbols('x y')
    chi = 1 - 3 * y ** 2 + 2 * y ** 3
    phi = chi * (sp_.Rational(al) + sp_.Rational(be) * x + sp_.Rational(ga) * y)
    txx, txy = sp_.diff(phi, y, 2), -sp_.diff(phi, x, y)
    Qs = sp_.integrate(txx ** 2 + sp_.Rational(13, 5) * txy ** 2, (x, -1, 1), (y, 0, 1))      # zero stress for y < 0
    ls = 100 * sp_.integrate(txx, (x, 0, sp_.Rational(1, 10)), (y, 0, sp_.Rational(1, 10)))
    assert out['tau0']['Q'] == Fr(str(Qs)) and out['tau0']['l'] == Fr(str(ls)) and out['tau0']['F'] == -Fr(be)
    print(f"   a(lift) = {out['v0']['a']}, Q(particular) = {out['tau0']['Q']}, l = {out['tau0']['l']}, F = {out['tau0']['F']}")
    head("   ok")

    # ---------------------------------------------------------------------------------------------------------------------
    head("2. exact identities of random admissible stress fields")
    Pr = S.stress_ties()
    z = Pr @ rng.standard_normal(Pr.shape[1])
    par, C = z[:3], z[3:].reshape(S.sx.n, S.sy.n)
    csx, csy = _Cache(S.sx, X), _Cache(S.sy, Y)
    php = CHI * (q(par[0]) + q(par[1]) * X + q(par[2]) * Y)
    pcs = {}
    for (ex, ey) in S.cells(S.sx, S.sy):
        phi = _spline_cell(csx, csy, S.sx, S.sy, C, ex, ey) + (php if S.sy.a[ey] >= 0 else 0)
        pcs[(ex, ey)] = (phi, phi.derivative('y').derivative('y'), phi.derivative('x').derivative('x'),
                         -phi.derivative('x').derivative('y'))
    nfree = 0
    for (ex, ey), (phi, txx, tyy, txy) in pcs.items():
        x0, x1, y0, y1 = S.sx.a[ex], S.sx.b[ex], S.sy.a[ey], S.sy.b[ey]
        assert (txx.derivative('x') + txy.derivative('y')).is_zero() and (txy.derivative('x') + tyy.derivative('y')).is_zero()
        free_edges = []
        if y1 == 1:
            free_edges.append(('y', 1))
        if y0 == -1 and x1 <= 0:
            free_edges.append(('y', -1))
        if x1 == 0 and y1 <= 0:
            free_edges.append(('x', 0))
        if y0 == 0 and x0 >= 0:
            free_edges.append(('y', 0))
        for var, val in free_edges:
            if var == 'y':      # traction (tau_xy, tau_yy)
                assert txy.subs({'y': q(val)}).is_zero() and tyy.subs({'y': q(val)}).is_zero(), (ex, ey, var, val)
            else:               # traction (tau_xx, tau_xy)
                assert txx.subs({'x': q(val)}).is_zero() and txy.subs({'x': q(val)}).is_zero(), (ex, ey, var, val)
            nfree += 1
        if (ex + 1, ey) in pcs:
            n = pcs[(ex + 1, ey)]
            assert (txx - n[1]).subs({'x': q(x1)}).is_zero() and (txy - n[3]).subs({'x': q(x1)}).is_zero()
        if (ex, ey + 1) in pcs:
            n = pcs[(ex, ey + 1)]
            assert (txy - n[3]).subs({'y': q(y1)}).is_zero() and (tyy - n[2]).subs({'y': q(y1)}).is_zero()
    # vertical grip force: F = -int_0^1 tau_xy(1, y) dy = [phi_x(1, y)]_0^1 must equal -beta
    ex = S.sx.ne - 1
    Fv = Fr(0)
    for ey in range(S.sy.ne):
        if S.sy.a[ey] >= 0:
            px = pcs[(ex, ey)][0].derivative('x').subs({'x': q(1)})
            for yy, sgn in ((S.sy.b[ey], 1), (S.sy.a[ey], -1)):
                v = px.subs({'y': q(yy)})
                Fv += sgn * (Fr(int(v.coefficient(0).p), int(v.coefficient(0).q)) if not v.is_zero() else 0)
    assert Fv == -Fr(float(par[1]))
    print(f"   all identities hold exactly ({len(pcs)} cells, {nfree} free cell edges)")
    head("   ok")

    # ---------------------------------------------------------------------------------------------------------------------
    head("3. exact evaluator against the float assembly (random admissible fields)")
    A, f, lt = S.assemble_dual()
    K_, b, k0, l0, lv = S.assemble_primal()
    Pv = S.disp_ties()
    v = Pv @ rng.standard_normal(Pv.shape[1])
    nA = S.ux.n * S.uy.n
    fields = dict(v0=dict(A=v[:nA].reshape(S.ux.n, S.uy.n), B=v[nA:].reshape(S.ux.n, S.uy.n)), tau0=dict(C=C, par=par))
    out = p2rig.evaluate(S, fields, nproc=2)
    pairs = [(out['tau0']['Q'], z @ A @ z), (out['tau0']['l'], lt @ z), (out['tau0']['F'], f @ z),
             (out['v0']['a'], v @ K_ @ v + 2 * b @ v + k0), (out['v0']['l'], l0 + lv @ v)]
    for r, fl in pairs:
        assert abs(float(r) - fl) <= 1e-9 * max(1.0, abs(float(r))), (float(r), fl)
    print("   " + "  ".join(f"{float(r):.12g}" for r, _ in pairs))
    head("   ok")

    if STEM:
        head(f"4. falsification against the stored certificate {STEM}")
        meta, arr = cio.load(STEM)
        S = pc.setup_from_meta(meta)
        Klo, Khi = Fr(meta['claims']['Q1']['lower']), Fr(meta['claims']['Q1']['upper'])
        Qlo, Qhi = Fr(meta['claims']['Q2']['lower']), Fr(meta['claims']['Q2']['upper'])
        kap = Fr(meta['fields']['Q2']['kappa'])
        Pr, Pv = S.stress_ties(), S.disp_ties()
        A, f, lt = S.assemble_dual()
        K_, b, k0, l0, lv = S.assemble_primal()
        vec = lambda n: np.concatenate([arr[n + '_A'].ravel(), arr[n + '_B'].ravel()])
        tec = lambda n: np.concatenate([arr[n + '_par'], arr[n + '_C'].ravel()])
        v0, vp, vm, z0, zp, zm = vec('v0'), vec('vp'), vec('vm'), tec('tau0'), tec('tp'), tec('tm')
        Qn, an = z0 @ A @ z0, v0 @ K_ @ v0
        nA = S.ux.n * S.uy.n

        def pert(base, P, M, nrm, eps):
            d = P @ rng.standard_normal(P.shape[1])
            return base + eps * d * np.sqrt(nrm / max(d @ M @ d, 1e-300))

        def dsp(v):
            return dict(A=v[:nA].reshape(S.ux.n, S.uy.n), B=v[nA:].reshape(S.ux.n, S.uy.n))

        def stf(z):
            return dict(par=z[:3], C=z[3:].reshape(S.sx.n, S.sy.n))

        nfail = 0
        for eps in (0.0, 1e-6, 1e-4, 1e-2):
            t = Fr(float(rng.uniform(0.25, 4.0))).limit_denominator(64)
            tf = float(t)
            fl = dict(v0=dsp(pert(v0, Pv, K_, an, eps)), tau0=stf(pert(z0, Pr, A, Qn, eps)),
                      vp=dsp(pert(v0 + tf * (vp - v0), Pv, K_, an, eps)), vm=dsp(pert(v0 - tf * (v0 - vm), Pv, K_, an, eps)),
                      tp=stf(pert(z0 + tf * (zp - z0), Pr, A, Qn, eps)), tm=stf(pert(z0 - tf * (z0 - zm), Pr, A, Qn, eps)))
            res = pc.bounds(p2rig.evaluate(S, fl), kap * t)
            ok = res['K_low'] <= Khi and res['K_up'] >= Klo and res['Q2_low'] <= Qhi and res['Q2_up'] >= Qlo
            print(f"   eps={eps:g}, kappa'={tf:.3f} kappa: K in [{float(res['K_low']):.10f}, {float(res['K_up']):.10f}]  "
                  f"Q2 in [{float(res['Q2_low']):.8f}, {float(res['Q2_up']):.8f}]  {'no crossing' if ok else 'CROSSING'}")
            nfail += not ok
        assert nfail == 0
        head("   ok")
        head("5. re-evaluation of the stored certificate")
        assert pc.verify(STEM)
    head("ALL P2 TESTS PASSED")


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    main()
