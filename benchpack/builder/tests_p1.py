"""Tests of the P1 builder (python tests_p1.py; about two minutes).

1. Known answers: the plate WITHOUT hole, frictionless grips: K = 1 and Q2 = 1/2 exactly; both certified intervals
   must contain them (and be narrow).
2. Exact identities of random admissible stress fields (python-flint, exact rationals): zero traction on the hole
   (divisibility by rho), on y = 1, zero shear on the mirror lines and the frictionless grip, traction continuity
   across cell edges, zero force of the Airy part.
3. The rigorous evaluator agrees with the float assembly (independent code paths) on random fields.
4. Falsification: random admissible fields, and random perturbations of the optimal ones, never give bounds that
   cross the certified interval of the stored certificate.
"""
import limits  # noqa: F401  (thread and process limits; must come before numpy)
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from flint import arb, fmpq

import p1cert as pc
import p1float as pf
from p1rig import Rig, RHO
from exact import q, X, Y
import certio as cio

sys.stdout.reconfigure(encoding='utf-8')
rng = np.random.default_rng(7)
T0 = time.time()
STEM = 'P1_frictionless_N16_p9'


def head(s):
    print(f"[{time.time() - T0:5.0f}s] {s}", flush=True)


# ---------------------------------------------------------------------------------------------------------------------
head("1. known answers: plate without hole, K = 1, Q2 = 1/2")
par = pc.params(4, 4, 2)
par['hole'] = False
meta = pc.build(None, par, verbose=False)
K, Q2 = meta['computed']['K'], meta['computed']['Q2']
print(f"   {K[0]!r} <= K <= {K[1]!r};   {Q2[0]!r} <= Q2 <= {Q2[1]!r}")
assert cio.ball_lo(meta['_res']['K_low']) <= 1 <= cio.ball_hi(meta['_res']['K_up'])
assert cio.ball_lo(meta['_res']['Q2_low']) <= Fr(1, 2) <= cio.ball_hi(meta['_res']['Q2_up'])
assert K[1] - K[0] < 1e-11 and Q2[1] - Q2[0] < 1e-4
head("   ok")

# ---------------------------------------------------------------------------------------------------------------------
head("2. exact identities of random admissible stress fields (N = 4, ps = 5)")
for grip in ('frictionless', 'clamped'):
    S = pc.setup_from_params(pc.params(4, 5, 3, grip))
    R = Rig(S)
    Pr = S.stress_ties()
    z = Pr @ rng.standard_normal(Pr.shape[1])
    gamma, C = pc.split_stress(S, z)
    gp = R.g_poly(gamma)
    pcs = {}
    for ex, ey, kind in R.scells:
        pcs[(ex, ey)] = R.stress_cell(C, gamma, ex, ey, gp)
    rx, ry = RHO.derivative('x'), RHO.derivative('y')
    nchk = 0
    for (ex, ey), (txx, tyy, txy) in pcs.items():
        x0, x1, y0, y1 = S.sx.a[ex], S.sx.b[ex], S.sy.a[ey], S.sy.b[ey]
        # equilibrium
        assert (txx.derivative('x') + txy.derivative('y')).is_zero() and (txy.derivative('x') + tyy.derivative('y')).is_zero()
        # hole
        if x0 * x0 + y0 * y0 < Fr(1, 4):
            for comp in (txx * rx + txy * ry, txy * rx + tyy * ry):
                qq, rr = divmod(comp, RHO)
                assert rr.is_zero(), "hole traction not divisible by rho"
            nchk += 1
        # free edge y = 1, mirror lines, grip
        if y1 == 1:
            assert txy.subs({'y': q(1)}).is_zero() and tyy.subs({'y': q(1)}).is_zero()
        if y0 == 0:
            assert txy.subs({'y': q(0)}).is_zero()
        if x0 == 0:
            assert txy.subs({'x': q(0)}).is_zero()
        if x1 == 1 and grip == 'frictionless':
            assert txy.subs({'x': q(1)}).is_zero()
        # continuity across the right and the upper edge
        if (ex + 1, ey) in pcs:
            n = pcs[(ex + 1, ey)]
            assert (txx - n[0]).subs({'x': q(x1)}).is_zero() and (txy - n[2]).subs({'x': q(x1)}).is_zero()
        if (ex, ey + 1) in pcs:
            n = pcs[(ex, ey + 1)]
            assert (txy - n[2]).subs({'y': q(y1)}).is_zero() and (tyy - n[1]).subs({'y': q(y1)}).is_zero()
    # Airy part carries no force through x = 1: int_0^1 phi_yy(1, y) dy = [phi_y(1, y)]_0^1 = 0
    ex = S.sx.ne - 1
    tot = fmpq(0)
    for ey in range(S.sy.ne):
        phi = R.airy_phi(C, ex, ey).derivative('y').subs({'x': q(1)})
        tot += phi.subs({'y': q(S.sy.b[ey])}).coefficient(0) if not phi.subs({'y': q(S.sy.b[ey])}).is_zero() else 0
        tot -= phi.subs({'y': q(S.sy.a[ey])}).coefficient(0) if not phi.subs({'y': q(S.sy.a[ey])}).is_zero() else 0
    assert tot == 0
    print(f"   [{grip}] all identities hold exactly ({len(pcs)} cells, {nchk} cut by the circle)")
head("   ok")

# ---------------------------------------------------------------------------------------------------------------------
head("3. rigorous evaluator against the float assembly (random admissible fields, N = 8)")
S = pc.setup_from_params(pc.params(8, 6, 5))
R = Rig(S)
A, f = S.assemble_dual()
K_, b, k0 = S.assemble_primal()
lt = S.patch_dual()
l0, lv = S.patch_primal()
for trial in range(2):
    z = S.stress_ties() @ rng.standard_normal(S.stress_ties().shape[1])
    v = S.disp_ties() @ rng.standard_normal(S.disp_ties().shape[1])
    g, C = pc.split_stress(S, z)
    Am, Bm = pc.split_disp(S, v)
    Qr, ar = R.Q(C, g), R.a(Am, Bm)
    ltr, lvr = R.l_tau(C, g), R.l_v(Am, Bm)
    pairs = [(float(Qr.mid()), z @ A @ z), (float(ar.mid()), v @ K_ @ v + 2 * b @ v + k0),
             (float(ltr.mid()), lt @ z), (float(lvr.mid()), l0 + lv @ v), (float(R.F(g)), f @ z)]
    for r, fl in pairs:
        assert abs(r - fl) <= 1e-9 * max(1.0, abs(r)), (r, fl)
    print("   " + "  ".join(f"{r:.12g}" for r, _ in pairs))
head("   ok")

# ---------------------------------------------------------------------------------------------------------------------
head(f"4. falsification against the stored certificate {STEM}")
meta, arr = cio.load(STEM)
S = pc.setup_from_meta(meta)
R = Rig(S)
Klo, Khi = Fr(meta['claims']['Q1']['lower']), Fr(meta['claims']['Q1']['upper'])
Qlo, Qhi = Fr(meta['claims']['Q2']['lower']), Fr(meta['claims']['Q2']['upper'])
kap = Fr(meta['fields']['Q2']['kappa'])
cst = arb(q(pf.C1111 / pf.PH_AREA)) / 2
Pr, Pv = S.stress_ties(), S.disp_ties()
A, f = S.assemble_dual()
K_, b, k0 = S.assemble_primal()
vec = lambda n: np.concatenate([arr[n + '_A'].ravel(), arr[n + '_B'].ravel()])
tec = lambda n: np.concatenate([arr[n + '_par'], arr[n + '_C'].ravel()])
v0, vp, vm, z0, zp, zm = vec('v0'), vec('vp'), vec('vm'), tec('tau0'), tec('tp'), tec('tm')
Qn, an = z0 @ A @ z0, v0 @ K_ @ v0


def pert(base, P, M, nrm, eps):
    """base + a random admissible perturbation of relative energy eps"""
    d = P @ rng.standard_normal(P.shape[1])
    return base + eps * d * np.sqrt(nrm / max(d @ M @ d, 1e-300))


nfail = 0
for eps in (0.0, 1e-7, 1e-6, 1e-5, 1e-4, 1e-2):
    zr, vr = pert(z0, Pr, A, Qn, eps), pert(v0, Pv, K_, an, eps)
    g, C = pc.split_stress(S, zr)
    Am, Bm = pc.split_disp(S, vr)
    F = arb(q(R.F(g)))
    Kl = cio.ball_hi(F * F / R.Q(C, g))                        # any admissible tau: K >= F^2/Q
    Ku = cio.ball_lo(4 * R.a(Am, Bm))                          # any admissible v:   K <= 4 a_q
    # Q2 at another kappa' (fields recombined from the stored ones, then perturbed)
    t = float(rng.uniform(0.25, 4.0))
    kp = arb(q(kap * Fr(t)))
    fp, fm = [pert(v0 + s * t * (vp - v0), Pv, K_, an, eps) for s in (1, -1)]
    gp_, gm_ = [pert(z0 + s * t * (zp - z0), Pr, A, Qn, eps) for s in (1, -1)]
    (Ap, Bp), (Am2, Bm2) = pc.split_disp(S, fp), pc.split_disp(S, fm)
    (gp, Cp), (gm, Cm) = pc.split_stress(S, gp_), pc.split_stress(S, gm_)
    mup_p = R.a(Ap, Bp) / 2 - kp * R.l_v(Ap, Bp)
    mup_m = R.a(Am2, Bm2) / 2 + kp * R.l_v(Am2, Bm2)
    mlow_p = arb(q(R.F(gp))) / 2 - R.Q(Cp, gp) / 2 - kp * R.l_tau(Cp, gp) - kp * kp * cst
    mlow_m = arb(q(R.F(gm))) / 2 - R.Q(Cm, gm) / 2 + kp * R.l_tau(Cm, gm) - kp * kp * cst
    q2l, q2u = cio.ball_hi((mlow_m - mup_p) / (2 * kp)), cio.ball_lo((mup_m - mlow_p) / (2 * kp))
    ok = Kl <= Khi and Ku >= Klo and q2l <= Qhi and q2u >= Qlo
    print(f"   eps={eps:g}, kappa'={t:.2f} kappa: K in [{float(Kl):.12f}, {float(Ku):.12f}]  "
          f"Q2 in [{float(q2l):.8f}, {float(q2u):.8f}]  {'no crossing' if ok else 'CROSSING'}")
    nfail += not ok
assert nfail == 0
head("   ok")

head("5. re-evaluation of the stored certificate")
assert pc.verify(STEM)
head("ALL P1 TESTS PASSED")
