#!/usr/bin/env python3
"""Known answers and falsification tests for the three-dimensional elasticity reference (docs/ELASTICITY_3D.md).

A. Problems whose minimum is known exactly as a rational number.  The grids, the coefficients and the compliance
   are chosen dyadic, so that the stored nodal values are exactly the data of the problem; the certified interval
   must then contain the exact minimum, and here both bounds must also equal it to rounding, because the exact
   stress and the exact displacement lie in the discrete spaces.
     E-T1  affine data (homogeneous strain), all faces Dirichlet; uniaxial tension of a bar with four free faces
     E-T2  a manufactured quadratic displacement with a constant body force; pure bending of a bar with four free
           faces
B. Falsification (E-T3).  Every statically admissible stress field and every kinematically admissible displacement
   must satisfy  LB(sigma) <= E[u].  Random fields, perturbed optimal fields and optimal fields of other grids are
   paired; a single pair with LB > E would expose an error in the theorem or in the evaluation.
C. The number of polynomial states for every arrangement of traction-free faces (r = 6 for a ring, 0 otherwise),
   with and without body force.

usage: python check_theorem3d.py [--quick]
"""
import itertools
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import elasticity3d as E

sys.stdout.reconfigure(encoding="utf-8")
quick = '--quick' in sys.argv
rng = np.random.default_rng(20261123)
fails = 0
lam, mu = Fr(2), Fr(1)                       # compliance: 3/8 on the diagonal, -1/8 beside it, 1 for the shears: dyadic
DD, DN, ND, NN = ('D', 'D'), ('D', 'N'), ('N', 'D'), ('N', 'N')
RING = (DD, NN, NN)
T0 = time.time()


def known(tag, P, exact, tol=1e-9):
    global fails
    out, U, st = E.certify(P)
    inside = Fr(out['lower']) <= exact <= Fr(out['upper'])
    tight = float(exact) - out['lower'] <= tol * (1 + abs(float(exact))) and out['upper'] - float(exact) <= tol * (1 + abs(float(exact)))
    fails += not (inside and tight)
    print(f"  {tag:52s} {out['lower']!r} <= {float(exact)!r} <= {out['upper']!r}   "
          f"{'ok' if inside and tight else ('**** NOT ENCLOSED ****' if not inside else 'NOT TIGHT')}", flush=True)


# ---- exact polynomials in three variables, for the manufactured solution -------------------------------------------------
def pder(p, var):
    r = {}
    for e, v in p.items():
        if e[var]:
            ne = list(e)
            ne[var] -= 1
            r[tuple(ne)] = r.get(tuple(ne), Fr(0)) + e[var] * v
    return r


def padd(p, q, c=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, Fr(0)) + c * v
    return r


def pint_prod(p, q, box):
    tot = Fr(0)
    for a, v in p.items():
        for c, w in q.items():
            t = v * w
            for d in range(3):
                n = a[d] + c[d] + 1
                t *= (box[d][1] ** n - box[d][0] ** n) / n
            tot += t
    return tot


def pfloat(p):
    return lambda X, Y, Z: sum(float(v) * X ** e[0] * Y ** e[1] * Z ** e[2] for e, v in p.items())


print("A. known minima (exact rational values)")
P0 = E.Problem3D(box=[(0, 1)] * 3, N=(1, 1, 1), lam=lam, mu=mu)
D, S = P0.D, P0.S
# E-T1a: u = A x + c
A = [[Fr(1, 8), Fr(-3, 16), Fr(1, 4)], [Fr(1, 2), Fr(1, 4), Fr(-1, 8)], [Fr(-1, 4), Fr(3, 8), Fr(1, 16)]]
e0 = [A[0][0], A[1][1], A[2][2], A[2][1] + A[1][2], A[2][0] + A[0][2], A[1][0] + A[0][1]]
ex = 2 * sum(D[m][n] * e0[m] * e0[n] for m in range(6) for n in range(6)) / 2
g = lambda X, Y, Z: tuple(float(A[a][0]) * X + float(A[a][1]) * Y + float(A[a][2]) * Z + 0.25 * a for a in range(3))
for p, N, fn in ((1, (2, 1, 1), 'beltrami'), (2, (2, 2, 1), 'beltrami'), (1, (2, 1, 2), 'maxwell'), (1, (1, 2, 2), 'morera')):
    known(f"E-T1 affine data, p={p}, N={N}, {fn}", E.Problem3D(box=[(0, 2), (0, 1), (0, 1)], N=N, lam=lam, mu=mu, degree=p, g=g, functions=fn), ex)
# E-T1b: stress (s, 0, 0, 0, 0, 0), four free faces: u = s (S00 x, S10 y, S20 z)
s = Fr(1, 2)
ex = 2 * s * s * S[0][0] / 2
g = lambda X, Y, Z: (float(S[0][0] * s) * X, float(S[1][0] * s) * Y, float(S[2][0] * s) * Z)
for p, N in ((1, (2, 1, 1)), (2, (2, 2, 2))):
    known(f"E-T1 uniaxial tension, four free faces, p={p}, N={N}", E.Problem3D(box=[(0, 2), (0, 1), (0, 1)], N=N, lam=lam, mu=mu, degree=p, bc=RING, g=g), ex)
# the same bar along the two other axes
known("E-T1 uniaxial tension along x2, p=1", E.Problem3D(box=[(0, 1), (0, 2), (0, 1)], N=(1, 2, 1), lam=lam, mu=mu, degree=1, bc=(NN, DD, NN),
      g=lambda X, Y, Z: (float(S[1][0] * s) * X, float(S[0][0] * s) * Y, float(S[1][0] * s) * Z)), ex)
known("E-T1 uniaxial tension along x3, p=1", E.Problem3D(box=[(0, 1), (0, 1), (0, 2)], N=(1, 1, 2), lam=lam, mu=mu, degree=1, bc=(NN, NN, DD),
      g=lambda X, Y, Z: (float(S[1][0] * s) * X, float(S[1][0] * s) * Y, float(S[0][0] * s) * Z)), ex)
# E-T2a: a quadratic displacement on the unit cube, b = -div(D e(u)) constant
u = [{(2, 0, 0): Fr(1, 4), (1, 1, 0): Fr(1, 8), (0, 0, 2): Fr(1, 8), (0, 1, 1): Fr(-1, 4), (0, 0, 1): Fr(1, 2)},
     {(0, 2, 0): Fr(1, 8), (1, 1, 0): Fr(-1, 4), (1, 0, 1): Fr(1, 2), (2, 0, 0): Fr(-1, 8), (1, 0, 0): Fr(1, 2)},
     {(0, 0, 2): Fr(-1, 4), (1, 0, 1): Fr(1, 8), (0, 1, 1): Fr(1, 4), (0, 2, 0): Fr(1, 8), (0, 1, 0): Fr(-1, 4)}]
e = [pder(u[0], 0), pder(u[1], 1), pder(u[2], 2), padd(pder(u[1], 2), pder(u[2], 1)), padd(pder(u[0], 2), pder(u[2], 0)),
     padd(pder(u[0], 1), pder(u[1], 0))]
sig = [dict() for _ in range(6)]
for m in range(6):
    for n in range(6):
        if D[m][n] != 0:
            sig[m] = padd(sig[m], e[n], D[m][n])
VG = [[0, 5, 4], [5, 1, 3], [4, 3, 2]]
b = []
for i in range(3):
    dv = {}
    for j in range(3):
        dv = padd(dv, pder(sig[VG[i][j]], j))
    assert all(k == (0, 0, 0) for k, v in dv.items() if v != 0)
    b.append(-dv.get((0, 0, 0), Fr(0)))
unit = [(0, 1)] * 3
energy = sum(D[m][n] * pint_prod(e[m], e[n], unit) for m in range(6) for n in range(6)) / 2 \
    - sum(b[a] * pint_prod(u[a], {(0, 0, 0): Fr(1)}, unit) for a in range(3))
uf = [pfloat(x) for x in u]
print(f"  (body force of the manufactured solution: {[str(x) for x in b]}; exact energy {energy})")
known("E-T2 quadratic displacement, body force, p=2", E.Problem3D(box=unit, N=(2, 2, 2), lam=lam, mu=mu, degree=2, b=b,
      g=lambda X, Y, Z: (uf[0](X, Y, Z), uf[1](X, Y, Z), uf[2](X, Y, Z))), energy)
# E-T2b: pure bending, stress (y, 0, ...) on (0, 2) x (-1/2, 1/2)^2, four free faces:
#   u = (S00 x y, -S00 x^2/2 + S10 (y^2 - z^2)/2, S10 y z);  minimum = 1/2 S00 int y^2 = S00/12
gb = lambda X, Y, Z: (float(S[0][0]) * X * Y, -float(S[0][0]) * X * X / 2 + float(S[1][0]) * (Y * Y - Z * Z) / 2, float(S[1][0]) * Y * Z)
known("E-T2 pure bending, four free faces, p=2", E.Problem3D(box=[(0, 2), (Fr(-1, 2), Fr(1, 2)), (Fr(-1, 2), Fr(1, 2))], N=(2, 2, 2), lam=lam, mu=mu,
      degree=2, bc=RING, g=gb), S[0][0] / 12)
print(f"  ({time.time() - T0:.0f} s)")

print("B. falsification: LB(sigma) <= E[u] for every admissible pair")
lamf, muf = Fr(3, 2), Fr(1)
Daniso = [[(lamf + 2 * muf if i == j else lamf) if (i < 3 and j < 3) else (muf if i == j else Fr(0)) for j in range(6)] for i in range(6)]
for (i, j), v in {(0, 3): Fr(1, 4), (0, 5): Fr(-1, 8), (1, 4): Fr(1, 5), (2, 5): Fr(1, 3), (3, 4): Fr(1, 8), (4, 5): Fr(-1, 6), (0, 1): Fr(1, 4)}.items():
    Daniso[i][j] += v
    Daniso[j][i] += v
gq = lambda X, Y, Z: (0.1 * X - 0.15 * Y + 0.3 * X * Z, 0.25 * X + 0.2 * Y * Z - 0.1, 0.05 * Z - 0.2 * X * Y + 0.1 * Y)
closest = np.inf


def perturbed(P, U, st, amps, reps):
    """lower bounds of perturbed stress fields and upper bounds of perturbed candidates"""
    lbs, ubs = [], []
    for amp in amps:
        for _ in range(reps):
            dens = rng.choice([0.05, 0.5, 1.0])
            st2 = E.Stress(P, st.X + amp * rng.standard_normal(st.X.shape) * (rng.random(st.X.shape) < dens),
                           st.par + amp * rng.standard_normal(st.par.shape))
            lbs.append(E.lower_bound(P, st2))
            U2 = U + amp * rng.standard_normal(U.shape) * (rng.random(U.shape) < dens)
            U2[:, P.fixed] = P.gval[:, P.fixed]
            ubs.append(E.upper_bound(P, U2)[1])
    return lbs, ubs


# 1. the block of E-T4: the datum is trilinear, so it is the same problem on every grid and for both degrees,
#    and every lower bound can be compared with every upper bound
lbs, ubs = [], []
grids = [((1, 1, 1), 'beltrami', 'resultants'), ((2, 2, 2), 'beltrami', 'resultants'), ((2, 1, 2), 'beltrami', 'all'),
         ((2, 2, 2), 'maxwell', 'all'), ((1, 2, 2), 'morera', 'all')]
if not quick:
    grids += [((3, 3, 3), 'beltrami', 'resultants'), ((3, 2, 2), 'maxwell', 'resultants')]
for N, fn, rule in grids:
    for p in (1, 2):
        P = E.Problem3D(box=[(0, 1)] * 3, N=N, lam=lamf, mu=muf, degree=p, bc=RING, g=lambda X, Y, Z: (0.1 * X, 0 * Y, 0 * Z),
                        functions=fn, states=rule)
        out, U, st = E.certify(P)
        ubs.append(out['upper'])
        if p == 1:
            lbs.append(out['lower'])
            l2, u2 = perturbed(P, U, st, (1e-3, 1e-2, 1e-1, 1.0), 1 if quick else 2)
            lbs += l2
            ubs += u2
        else:
            ubs += perturbed(P, U, st, (1e-3, 1e-1), 1)[1]
top, bot = max(lbs), min(ubs)
closest = min(closest, bot - top)
fails += not (top <= bot)
print(f"  block of E-T4, all grids and fields paired         highest LB {top:+.8f}   lowest E {bot:+.8f}   "
      f"({len(lbs)} stress fields, {len(ubs)} candidates)   {'ok' if top <= bot else '**** VIOLATION ****'}   ({time.time() - T0:.0f} s)", flush=True)

# 2. problems whose datum depends on the grid: pairs within one instance
problems = [
    ("all Dirichlet, body force", dict(box=[(0, 2), (0, 1), (0, 1)], lam=lamf, mu=muf, b=(Fr(1, 2), Fr(-1, 3), Fr(1, 5)), g=gq), [(2, 1, 1), (2, 2, 2)]),
    ("five free faces, anisotropic, body force", dict(box=[(Fr(1, 3), 1), (0, 1), (0, 1)], D=Daniso, bc=(DN, NN, NN), b=(Fr(1, 5), Fr(-2, 3), Fr(1, 4)), g=gq),
     [(2, 1, 1), (2, 2, 2)]),
    ("two opposite free faces", dict(box=[(0, 1)] * 3, lam=lamf, mu=muf, bc=(DD, NN, DD), b=(Fr(0), Fr(1), Fr(0)), g=gq), [(1, 2, 1), (2, 2, 2)]),
    ("ring, twisted and sheared end, body force", dict(box=[(0, 1)] * 3, lam=lamf, mu=muf, bc=RING, b=(Fr(1, 4), Fr(0), Fr(-1, 2)),
                                                       g=lambda X, Y, Z: (0.02 * X * Y, 0.1 * X - 0.1 * X * (Z - 0.5), 0.1 * X * (Y - 0.5))), [(2, 1, 1), (2, 2, 2)]),
]
for tag, kw, Ns in problems:
    worst = np.inf
    for N in (Ns[:1] if quick else Ns):
        for p in (1, 2):
            P = E.Problem3D(N=N, degree=p, **kw)
            out, U, st = E.certify(P)
            l2, u2 = perturbed(P, U, st, (1e-3, 1e-2, 1e-1, 1.0), 1 if (quick or p == 2) else 2)
            worst = min(worst, min(u2 + [out['upper']]) - max(l2 + [out['lower']]))
    closest = min(closest, worst)
    fails += not (worst >= 0)
    print(f"  {tag:50s} smallest (E - LB) over all pairs {worst:+.3e}   {'ok' if worst >= 0 else '**** VIOLATION ****'}   ({time.time() - T0:.0f} s)", flush=True)

print("C. number of polynomial states by arrangement of the traction-free faces (expected: 6 for a ring, 0 otherwise)")
counts = {}
bad = 0
for pattern in itertools.product(('D', 'N'), repeat=6):
    bc = (pattern[0:2], pattern[2:4], pattern[4:6])
    if 'D' not in pattern:
        continue
    for bf in ((0, 0, 0), (Fr(1, 2), Fr(-1, 3), Fr(1, 5))):
        P = E.Problem3D(box=[(0, 2), (Fr(-1, 2), Fr(1, 2)), (Fr(1, 3), 1)], N=(1, 1, 1), lam=lamf, mu=muf, bc=bc, b=bf)
        psi_b, psis = P.states()
        ring = P.ring_axis() is not None
        bad += (len(psis) != (6 if ring else 0))
        key = (pattern.count('N'), ring, len(psis), len(P._all_states[1]))
        counts[key] = counts.get(key, 0) + 1
for (nfree, ring, r, rall), c in sorted(counts.items()):
    print(f"  {nfree} free faces{' (ring)' if ring else '       '}: r = {r}   (independent of the polynomial null sextuples: {rall})   {c} cases")
fails += bad
print(f"  a polynomial particular field exists in all {sum(counts.values())} cases; counts {'as expected' if bad == 0 else 'NOT AS EXPECTED'}")

print(f"\n{'ALL TESTS PASSED' if fails == 0 else 'FAILURES: %d' % fails}   smallest (E - LB) = {closest:.3e}   ({time.time() - T0:.0f} s)")
sys.exit(1 if fails else 0)
