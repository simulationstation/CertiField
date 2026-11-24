#!/usr/bin/env python3
"""Known answers and falsification tests for the elasticity reference (docs/ELASTICITY.md).

A. Problems whose minimum is known exactly as a rational number.  The grids, the coefficients and the compliance
   are chosen dyadic, so that the stored nodal values are exactly the data of the problem; the certified interval
   must then contain the exact minimum, and here both bounds must also equal it to rounding, because the exact
   stress and the exact displacement lie in the discrete spaces.
     E-T1  affine data (homogeneous strain), all faces Dirichlet; and uniaxial tension with free sides
     E-T2  manufactured quadratic displacements with a constant body force; pure bending with free sides
B. Falsification (E-T3).  Every statically admissible stress field and every kinematically admissible displacement
   must satisfy  LB(sigma) <= E[u].  Random fields, perturbed optimal fields and optimal fields of other grids are
   paired; a single pair with LB > E would expose an error in the theorem or in the evaluation.
"""
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import elasticity as E

rng = np.random.default_rng(20261123)
fails = 0
lam, mu = Fr(1), Fr(1)                       # compliance S = [[3/8, -1/8, 0], [-1/8, 3/8, 0], [0, 0, 1]]: dyadic


def known(tag, P, exact, tol=1e-10):
    global fails
    out, U, st = E.certify(P)
    inside = Fr(out['lower']) <= exact <= Fr(out['upper'])
    tight = float(exact) - out['lower'] <= tol * (1 + abs(float(exact))) and out['upper'] - float(exact) <= tol * (1 + abs(float(exact)))
    fails += not (inside and tight)
    print(f"  {tag:46s} {out['lower']!r} <= {float(exact)!r} <= {out['upper']!r}   "
          f"{'ok' if inside and tight else ('**** NOT ENCLOSED ****' if not inside else 'NOT TIGHT')}", flush=True)


print("A. known minima (exact rational values)")
D = E.Problem(box=[(0, 1), (0, 1)], N=(1, 1), lam=lam, mu=mu).D
S = E.Problem(box=[(0, 1), (0, 1)], N=(1, 1), lam=lam, mu=mu).S
# E-T1a: u = A x + c
e0 = [Fr(1, 8), Fr(1, 4), Fr(-3, 16) + Fr(1, 2)]
ex = 2 * sum(D[m][n] * e0[m] * e0[n] for m in range(3) for n in range(3)) / 2
for p in (1, 2):
    known(f"E-T1 affine data, p={p}", E.Problem(box=[(0, 2), (0, 1)], N=(4, 2), lam=lam, mu=mu, degree=p,
                                                g=lambda X, Y: (0.125 * X - 0.1875 * Y + 0.5, 0.5 * X + 0.25 * Y - 0.25)), ex)
# E-T1b: stress (s, 0, 0), free top and bottom: u = (S00 s x, S10 s y)
s = Fr(1, 2)
ex = 2 * s * s * S[0][0] / 2
for p in (1, 2):
    known(f"E-T1 uniaxial tension, free sides, p={p}", E.Problem(box=[(0, 2), (0, 1)], N=(4, 2), lam=lam, mu=mu, degree=p, bc=(('D', 'D'), ('N', 'N')),
                                                                 g=lambda X, Y: (float(S[0][0] * s) * X, float(S[1][0] * s) * Y)), ex)
# E-T2a: u1 = x^2/4 + x y/8, u2 = y^2/8 - x y/4 + x/2 on the unit square; b = -div(D e(u))
#   e = (x/2 + y/8, y/4 - x/4, x/8 - y/4 + 1/2)
ecoef = [(Fr(0), Fr(1, 2), Fr(1, 8)), (Fr(0), Fr(-1, 4), Fr(1, 4)), (Fr(1, 2), Fr(1, 8), Fr(-1, 4))]      # (const, x, y)
scoef = [tuple(sum(D[m][n] * ecoef[n][k] for n in range(3)) for k in range(3)) for m in range(3)]
b = (-(scoef[0][1] + scoef[2][2]), -(scoef[2][1] + scoef[1][2]))


def integ_unit(c1, c2):
    """int over the unit square of (a0 + a1 x + a2 y)(b0 + b1 x + b2 y)"""
    mom = {(0, 0): Fr(1), (1, 0): Fr(1, 2), (0, 1): Fr(1, 2), (2, 0): Fr(1, 3), (0, 2): Fr(1, 3), (1, 1): Fr(1, 4)}
    pw = [(0, 0), (1, 0), (0, 1)]
    return sum(c1[i] * c2[j] * mom[(pw[i][0] + pw[j][0], pw[i][1] + pw[j][1])] for i in range(3) for j in range(3))


strain_energy = sum(D[m][n] * integ_unit(ecoef[m], ecoef[n]) for m in range(3) for n in range(3)) / 2
# int u1 = 1/12 + 1/32, int u2 = 1/24 - 1/16 + 1/4
work = b[0] * (Fr(1, 12) + Fr(1, 32)) + b[1] * (Fr(1, 24) - Fr(1, 16) + Fr(1, 4))
known("E-T2 quadratic displacement, body force, p=2", E.Problem(box=[(0, 1), (0, 1)], N=(2, 2), lam=lam, mu=mu, degree=2, b=b,
      g=lambda X, Y: (X * X / 4 + X * Y / 8, Y * Y / 8 - X * Y / 4 + X / 2)), strain_energy - work)
# E-T2b: pure bending, stress (y, 0, 0) on (0, 2) x (-1/2, 1/2), free top and bottom:
#   u1 = S00 x y, u2 = S10 y^2/2 - S00 x^2/2;  minimum = 1/2 S00 int y^2 = S00/12
known("E-T2 pure bending, free sides, p=2", E.Problem(box=[(0, 2), (Fr(-1, 2), Fr(1, 2))], N=(4, 2), lam=lam, mu=mu, degree=2, bc=(('D', 'D'), ('N', 'N')),
      g=lambda X, Y: (float(S[0][0]) * X * Y, float(S[1][0]) * Y * Y / 2 - float(S[0][0]) * X * X / 2)), S[0][0] / 12)

print("B. falsification: LB(sigma) <= E[u] for every admissible pair")
Daniso = [[Fr(4), Fr(1), Fr(1, 2)], [Fr(1), Fr(3), Fr(-1, 4)], [Fr(1, 2), Fr(-1, 4), Fr(2)]]
gq = lambda X, Y: (0.1 * X - 0.15 * Y + 0.3 * X * Y, 0.25 * X + 0.2 * Y * Y - 0.1)
problems = [
    ("block: clamped, pulled, free sides", lambda N, p: E.Problem(box=[(0, 1), (0, 1)], N=N, lam=Fr(3, 2), mu=Fr(1), degree=p, bc=(('D', 'D'), ('N', 'N')), g=lambda X, Y: (0.1 * X, 0 * Y))),
    ("all Dirichlet, body force", lambda N, p: E.Problem(box=[(0, 2), (0, 1)], N=N, lam=Fr(3, 2), mu=Fr(1), degree=p, b=(Fr(1, 2), Fr(-1, 3)), g=gq)),
    ("three free faces, anisotropic, body force", lambda N, p: E.Problem(box=[(Fr(1, 3), 1), (0, 1)], N=N, D=Daniso, degree=p, bc=(('D', 'N'), ('N', 'N')), b=(Fr(1, 5), Fr(-2, 3)), g=gq)),
    ("two adjacent free faces", lambda N, p: E.Problem(box=[(0, 1), (0, 1)], N=N, lam=Fr(1), mu=Fr(1), degree=p, bc=(('N', 'D'), ('N', 'D')), b=(Fr(1), Fr(1)), g=gq)),
]
closest = np.inf
t0 = time.time()
for tag, mk in problems:
    lbs, ubs = [], []
    for N in ((2, 2), (4, 4), (8, 8), (16, 16)):
        for p in (1, 2):
            P = mk(N, p)
            out, U, st = E.certify(P)
            if p == 1:                                   # the stress field does not depend on p when the data are in Q_1 ...
                lbs.append(out['lower'])
            ubs.append(out['upper'])
            if N[0] <= 8:
                for amp in (1e-3, 1e-2, 1e-1, 1.0):
                    for _ in range(4):
                        st2 = E.Stress(P, st.X + amp * rng.standard_normal(st.X.shape) * (rng.random(st.X.shape) < rng.choice([0.05, 0.5, 1.0])),
                                       st.par + amp * rng.standard_normal(st.par.shape))
                        if p == 1 or tag.startswith("block"):
                            lbs.append(E.lower_bound(P, st2))
                        U2 = U + amp * rng.standard_normal(U.shape) * (rng.random(U.shape) < rng.choice([0.05, 0.5, 1.0]))
                        U2[:, P.fixed] = P.gval[:, P.fixed]
                        ubs.append(E.upper_bound(P, U2)[1])
    # the data of these problems are Q_1 only for the block; for the others the datum depends on (N, p), so the pairs
    # are formed within one problem instance only
    if tag.startswith("block"):
        top, bot = max(lbs), min(ubs)
        ok = top <= bot
        closest = min(closest, bot - top)
        print(f"  {tag:44s} highest LB {top:+.8f}   lowest E {bot:+.8f}   {'ok' if ok else '**** VIOLATION ****'}", flush=True)
        fails += not ok
    else:
        worst = np.inf
        for N in ((2, 2), (4, 4), (8, 8)):
            for p in (1, 2):
                P = mk(N, p)
                out, U, st = E.certify(P)
                lbs2, ubs2 = [out['lower']], [out['upper']]
                for amp in (1e-3, 1e-2, 1e-1, 1.0):
                    for _ in range(6):
                        st2 = E.Stress(P, st.X + amp * rng.standard_normal(st.X.shape), st.par + amp * rng.standard_normal(st.par.shape))
                        lbs2.append(E.lower_bound(P, st2))
                        U2 = U + amp * rng.standard_normal(U.shape)
                        U2[:, P.fixed] = P.gval[:, P.fixed]
                        ubs2.append(E.upper_bound(P, U2)[1])
                worst = min(worst, min(ubs2) - max(lbs2))
        ok = worst >= 0
        closest = min(closest, worst)
        print(f"  {tag:44s} smallest (E - LB) over all pairs {worst:+.3e}   {'ok' if ok else '**** VIOLATION ****'}", flush=True)
        fails += not ok

print(f"\n{'ALL TESTS PASSED' if fails == 0 else 'FAILURES: %d' % fails}   smallest (E - LB) = {closest:.3e}   ({time.time() - t0:.0f} s)")
sys.exit(1 if fails else 0)
