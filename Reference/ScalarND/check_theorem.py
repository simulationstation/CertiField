#!/usr/bin/env python3
"""Falsification tests of the lower bound (theorem + implementation).  For problems whose exact minimum is known
independently, the certified value  LB(S)  must stay below it for EVERY field S: optimised fields, randomly perturbed
fields, and fields produced by an optimiser that is allowed to cheat between the samples.  Any LB(S) above the
reference would expose an error in the proof or in the rigorous evaluation."""
import sys, time, numpy as np
from fractions import Fraction as Fr
import scalarnd as S
from ref1d import wall_energy
from refseries import screened_min

rng = np.random.default_rng(20261123)
W4 = [Fr(1, 4), 0, Fr(-1, 2), 0, Fr(1, 4)]
fails = 0
worst = -np.inf


def attack(tag, P, cal, ref, nper=60, amps=(1e-3, 1e-2, 1e-1, 1.0), Qs=16):
    global fails, worst
    x0 = cal.pack()
    lb0 = S.lower_bound(P, cal, Qs=Qs)
    top = lb0
    for a in amps:
        for _ in range(nper):
            sc = np.abs(x0).max() + 1e-3
            x = x0 + a * sc * rng.standard_normal(x0.size) * (rng.random(x0.size) < rng.choice([0.02, 0.2, 1.0]))
            cal.unpack(x)
            lb = S.lower_bound(P, cal, Qs=Qs)
            top = max(top, lb)
    cal.unpack(x0)
    ok = top <= ref
    fails += (not ok)
    worst = max(worst, top - ref)
    print(f"  {tag:34s} optimised LB {lb0:+.8f}  highest LB over perturbations {top:+.8f}  reference {ref:+.8f}  "
          f"{'ok' if ok else '**** VIOLATION ****'}", flush=True)


t0 = time.time()
print("A. 1D double well on (-3,3), exact minimum from the first integral")
ref = wall_energy(3)[0]
P = S.Problem(box=[(-3, 3)], N=(48,), k=Fr(1, 2), W=W4, bc=[('D', 'D')], g=lambda X: np.where(X < 0, -1.0, 1.0))
for K, Qo, curv in ((16, 4, True), (16, 1, False), (8, 1, False), (32, 2, False)):
    cal = S.Calibration(P, -1, 1, K, affine=False)
    S.optimise_calibration(P, cal, Qo=Qo, taus=(1e-2, 1e-3, 1e-4, 1e-5), iters=400, curvature=curv)
    attack(f"1D K={K} Qo={Qo} curvature={curv}", P, cal, ref, Qs=16)
    attack(f"   same field, verifier Qs=2", P, cal, ref, nper=10, Qs=2)
print("   mixed field (affine + spline optimised together)")
cal = S.Calibration(P, -1, 1, 8, affine=True)
S.optimise_calibration(P, cal, Qo=2, taus=(1e-2, 1e-3, 1e-4), iters=600, curvature=False)
attack("1D mixed K=8", P, cal, ref)

print("B. 2D and 3D double-well walls (exact minimum = area x 1D value)")
P = S.Problem(box=[(-3, 3), (0, 1)], N=(24, 3), k=Fr(1, 2), W=W4, bc=[('D', 'D'), ('N', 'N')], g=lambda X, Y: np.where(X < 0, -1.0, 1.0))
cal = S.Calibration(P, -1, 1, 8, affine=False)
S.optimise_calibration(P, cal, Qo=1, taus=(1e-2, 1e-3, 1e-4), iters=500, curvature=False)
attack("2D wall K=8 (cheating optimiser)", P, cal, ref, nper=30)
P = S.Problem(box=[(-3, 3), (0, 1), (0, 2)], N=(16, 2, 2), k=Fr(1, 2), W=W4, bc=[('D', 'D'), ('N', 'N'), ('N', 'N')], g=lambda X, Y, Z: np.where(X < 0, -1.0, 1.0))
cal = S.Calibration(P, -1, 1, 8, affine=False)
S.optimise_calibration(P, cal, Qo=1, taus=(1e-2, 1e-3, 1e-4), iters=500, curvature=False)
attack("3D wall K=8, area 2", P, cal, 2 * ref, nper=20)

print("C. convex problems with series solutions (affine fields, random perturbations)")
for n, N in ((2, 12), (3, 6)):
    ref = screened_min(n, 0.5, 1.0, 1.0)
    P = S.Problem(box=[(0, 1)] * n, N=(N,) * n, k=Fr(1, 2), W=[0, -1, 1], bc=[('D', 'D')] * n, g=lambda *X: 0 * X[0])
    out, U, cal = S.certify_convex(P)
    attack(f"{n}D screened Poisson, convex mode", P, cal, ref, nper=40, amps=(1e-3, 1e-2, 1e-1))
    cal2 = S.Calibration(P, Fr(-1, 10), Fr(3, 10), 4, affine=True)
    cal2.S0 = [s.copy() for s in cal.S0]
    S.optimise_calibration(P, cal2, Qo=1, taus=(1e-3, 1e-4), iters=300, curvature=False)
    attack(f"{n}D screened Poisson, mixed K=4", P, cal2, ref, nper=20, amps=(1e-3, 1e-2, 1e-1))

print("D. non-constant boundary data: every Q1 competitor on every grid must lie above the bound")
lam = 16
W16 = [Fr(lam, 4), 0, Fr(-lam, 2), 0, Fr(lam, 4)]
mk = lambda N: S.Problem(box=[(-1, 1), (-1, 1)], N=(N, N), k=Fr(1, 2), W=W16, bc=[('D', 'D'), ('D', 'D')], g=lambda X, Y: X)
best = np.inf
for N in (16, 32, 64, 128):
    Pf = mk(N)
    X, Y = np.meshgrid(*Pf.axes, indexing='ij')
    for trial in range(4 if N <= 32 else 1):
        U = S.solve_candidate(Pf, X + (0.5 * rng.standard_normal(X.shape) if trial else 0.0))
        best = min(best, S.upper_bound(Pf, U)[1])
P = mk(16)
cal = S.Calibration(P, -1, 1, 8, affine=False)
S.optimise_calibration(P, cal, Qo=1, taus=(1e-2, 1e-3, 1e-4), iters=600, scale=4, curvature=False)
attack("2D ramp (g = x), K=8", P, cal, best, nper=20)

print(f"\n{'ALL TESTS PASSED' if fails == 0 else 'FAILURES: %d' % fails}   largest (LB - reference) = {worst:.3e}   ({time.time() - t0:.0f} s)")
sys.exit(1 if fails else 0)
