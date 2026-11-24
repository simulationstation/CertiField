#!/usr/bin/env python3
"""tests.py -- checks of the construction behind the 3D bounds (run: python3 tests.py; about 5 minutes)"""
import time
from fractions import Fraction as Fr
import numpy as np
import curved3d as c3
from curved3d import *

t0 = time.time()
rng = np.random.default_rng(0)
H = Fr(1, 2)


def problem(N, grip, L):
    return Problem(N, grip=grip, particular='step', column=H, kflat=(3 * N // 4 if grip == 'frictionless' else None), L_max=L)


def full(t6):
    S = np.zeros(t6.shape[:-1] + (3, 3))
    for s, (i, j) in enumerate(COMPS):
        S[..., i, j] = t6[..., s]
        S[..., j, i] = t6[..., s]
    return S


def tau_at(P, cfull, pts):
    """float stress of inc(rho^2 Psi) at arbitrary points (assembly code path)"""
    pts = np.asarray(pts, float)
    e = np.minimum((pts * P.N).astype(int), P.N - 1)
    t = pts * P.N - e
    tabs = [float_tables_at(P.N, P.funcs, e[:, k], t[:, k]) for k in range(3)]
    T, taup = dual_core_float(P, tabs[0][0], tabs[1][0], tabs[2][0], tabs[0][1], tabs[1][1], tabs[2][1])
    bx = np.zeros((len(pts), 8), dtype=np.int64)
    bx[:, :3] = e
    cl = c3._local_coeffs(cfull, bx).reshape(len(pts), 162)
    return np.einsum('bsn,bn->bs', T[:, 0], cl), taup[:, 0]


print("1. Gauss rules: rational enclosures against numpy")
for n in (3, 7):
    nodes, weights = gauss_enclosures(n)
    x, w = np.polynomial.legendre.leggauss(n)
    assert max(abs(float(a.lo) - (1 + b) / 2) for a, b in zip(nodes, x)) < 1e-15
    assert max(abs(float(a.lo) - b / 2) for a, b in zip(weights, w)) < 1e-15
    assert all(a.lo <= a.hi and a.hi - a.lo < Fr(1, 10 ** 30) for a in nodes + weights)
print("   ok")

for grip in ('clamped', 'frictionless'):
    N = 8
    P = problem(N, grip, 4)
    Pp, own = psi_ties(N, grip, owners=True)
    Pv, ownv = v_ties(N, grip, owners=True)
    cp = (Pp @ rng.standard_normal(Pp.shape[1])).reshape(6, P.M, P.M, P.M)
    cv = (Pv @ rng.standard_normal(Pv.shape[1])).reshape(3, P.M, P.M, P.M)
    print(f"2. [{grip}] random potential: equilibrium and boundary conditions of inc(rho^2 Psi)")
    assert ties_hold(Pp, own, cp) and ties_hold(Pv, ownv, cv)
    pts = (rng.integers(0, N, size=(300, 3)) + rng.uniform(0.1, 0.9, size=(300, 3))) / N      # away from cell faces
    S0 = np.abs(tau_at(P, cp, pts)[0]).max()
    errs = []
    for h in (1e-4, 1e-5):
        div = np.zeros((len(pts), 3))
        for k in range(3):
            dp = np.zeros(3)
            dp[k] = h
            div += (full(tau_at(P, cp, pts + dp)[0])[:, :, k] - full(tau_at(P, cp, pts - dp)[0])[:, :, k]) / (2 * h)
        errs.append(np.abs(div).max())
        print(f"   max|tau| {S0:.1f}; max|div tau| by central differences, step {h:g}: {errs[-1]:.2e}")
    assert errs[1] < 1e-6 * S0 and errs[1] < 0.05 * errs[0], "the difference quotient must vanish like the step squared"
    u = np.abs(rng.standard_normal((400, 3)))
    u /= np.linalg.norm(u, axis=1)[:, None]
    S = full(tau_at(P, cp, 0.5 * u)[0])
    tr = np.abs(np.einsum('bij,bj->bi', S, u)).max()
    print(f"   sphere: max|tau| {np.abs(S).max():.2f}, max|tau n| {tr:.1e}")
    assert tr < 1e-11
    for k in (1, 2):
        pf = rng.uniform(0.02, 0.98, size=(300, 3))
        pf[:, k] = 1 - 1e-12
        S = full(tau_at(P, cp, pf)[0])
        assert np.abs(S[:, :, k]).max() < 1e-6
    for k in range(3):
        pf = rng.uniform(0.02, 0.98, size=(300, 3))
        pf[:, k] = 1e-12
        S = full(tau_at(P, cp, pf)[0])
        assert np.abs(np.delete(S[:, :, k], k, axis=1)).max() < 1e-6
    if grip == 'frictionless':
        pf = rng.uniform(0.02, 0.98, size=(300, 3))
        pf[:, 0] = 1 - 1e-12
        S = full(tau_at(P, cp, pf)[0])
        sh = np.abs(S[:, 0, 1:]).max()
        print(f"   grip face: max|tau| {np.abs(S).max():.1f}, max shear traction {sh:.1e}")
        assert sh < 1e-6
    print("   free faces, mirror planes: tractions vanish  ok")

    print(f"3. [{grip}] two code paths give the same stress; the Gauss rules are exact for the integrands")
    bx = P.boxes(2)[::7][:200]
    tb = P.tables(7)
    Bx, By, Bz, fnx, fny, fnz, W = c3._gather_float(P, tb, bx)
    T, taup = dual_core_float(P, Bx, By, Bz, fnx, fny, fnz)
    cl = c3._local_coeffs(cp, bx).reshape(len(bx), 162)
    tA = np.einsum('bQsn,bn->bQs', T, cl) + taup
    tB = c3.dual_stress(P, tb, bx, cp, rigorous=False)
    dmax = max(np.abs(np.broadcast_to(tB[s], (len(bx), 7, 7, 7)).reshape(len(bx), -1) - tA[:, :, s]).max() for s in range(6))
    print(f"   max difference between the assembly path and the certificate path: {dmax:.1e} (stress scale {np.abs(tA).max():.1f})")
    assert dmax < 1e-9 * np.abs(tA).max()
    q7 = dual_compliance(P, bx, cp, rigorous=False, nq=7)
    q9 = dual_compliance(P, bx, cp, rigorous=False, nq=9)
    e3 = primal_energy(P, bx, cv, rigorous=False, nq=3)
    e5 = primal_energy(P, bx, cv, rigorous=False, nq=5)
    print(f"   stress integral with 7 and 9 points: relative difference {abs(q7 - q9) / q9:.1e};  strain integral with 3 and 5 points: {abs(e3 - e5) / e5:.1e}")
    assert abs(q7 - q9) < 1e-11 * q9 and abs(e3 - e5) < 1e-11 * e5
    qi = dual_compliance(P, bx, cp, rigorous=True)
    ei = primal_energy(P, bx, cv, rigorous=True)
    print(f"   interval upper bounds exceed the float values by {(qi - q7) / q7:.1e} and {(ei - e3) / e3:.1e} (relative)")
    assert 0 <= (qi - q7) / q7 < 1e-6 and 0 <= (ei - e3) / e3 < 1e-9

    print(f"4. [{grip}] Green identity with random fields:  int_Omega tau:eps(v) = resultant x grip displacement")
    for L in (2, 3, 4):
        bxa = P.boxes(L)
        res = {}
        for name, cf, sc in (("particular", None, 1.0), ("potential", cp, 0.0)):
            tin = tcut = mag = 0.0
            for s0 in range(0, len(bxa), 2048):
                b_ = bxa[s0:s0 + 2048]
                tau = c3.dual_stress(P, tb, b_, cf, rigorous=False, scale_p=sc)
                e = c3.primal_strain(P, tb, b_, cv, rigorous=False)
                w = c3._weights_hi(tb, b_, False)
                dens = 0.0
                for k in range(6):
                    if tau[k] is not None:
                        dens = dens + tau[k] * e[k]
                dens = np.broadcast_to(dens, w.shape)
                val = np.sum((dens * w).reshape(len(b_), -1), axis=1)
                mag += np.sum(np.abs(dens) * w)
                tin += val[b_[:, 7] == 0].sum()
                tcut += val[b_[:, 7] == 1].sum()
            res[name] = (tin, tin + tcut, mag)
        a, b = res["particular"], res["potential"]
        print(f"   depth {L}: particular part {a[0]:.9f} (expected 0.5);  potential part between {min(b[0], b[1]):+.2e} and {max(b[0], b[1]):+.2e}"
              f" (expected 0; integral of |.| is {b[2]:.0f})")
        assert abs(a[0] - 0.5) < 1e-9 and abs(a[1] - 0.5) < 1e-9
    assert abs(b[0]) < 1e-5 * b[2] and abs(b[1]) < 1e-5 * b[2]
    print(f"   ok   {time.time() - t0:.0f}s", flush=True)

print("5. known answers")
P = Problem(8, a=0, grip='frictionless', particular='uniform', L_max=2, kflat=6)
r = certify(P, L_solve=0, L_cert=0, verbose=False)
print(f"   solid cube between frictionless grips (exact K = 2):   {r['K_low']:.13f} <= K <= {r['K_up']:.13f}")
assert r['K_low'] <= 2 <= r['K_up'] and r['K_up'] - r['K_low'] < 1e-9
R = Fr(3, 4)
pi_lo, pi_hi = Fr(314159265358979, 10 ** 14), Fr(314159265358980, 10 ** 14)
exact = float(np.pi * R * R / 2)
for L in (2, 3):
    P = Problem(8, a=R, kappa=(0, 1, 1), solid_inside=True, grip='frictionless', particular='uniform',
                area=(pi_lo * R * R, pi_hi * R * R), L_max=L)
    r = certify(P, L_solve=2, L_cert=L, verbose=False)
    print(f"   round bar, radius 3/4 (exact K = pi R^2/2 = {exact:.6f}), depth {L}:   {r['K_low']:.6f} <= K <= {r['K_up']:.6f}")
    assert r['K_low'] <= exact <= r['K_up']
print(f"ALL TESTS PASSED   {time.time() - t0:.0f}s")
