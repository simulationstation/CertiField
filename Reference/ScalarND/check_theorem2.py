#!/usr/bin/env python3
"""Known answers for the extensions of docs/SCALAR_ND_2.md (Section 7 there), with fields and candidates written
down by hand, so that only the rigorous evaluators are exercised (no solver, no optimiser).
In every case the exact minimum must lie between the certified lower bound and the certified upper bound.
usage: python check_theorem2.py        (seconds)

Not yet here: falsification with random fields, the calibration-mode cases, tables.  See PORTING.md, Section 10."""
import math
import sys
from fractions import Fraction as Fr

import numpy as np

import scalarnd as S

fails = 0


def report(tag, lb, exact, ub, tol_lb=None):
    global fails
    inside = Fr(lb) <= exact <= Fr(ub)
    tight = tol_lb is None or float(exact) - lb <= tol_lb
    fails += not (inside and tight)
    print(f"  {tag:58s} {lb!r} <= {float(exact)!r} <= {ub!r}   "
          f"{'ok' if inside and tight else ('**** NOT ENCLOSED ****' if not inside else 'LOWER BOUND NOT TIGHT')}", flush=True)


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    k = Fr(1, 2)

    print("1. source f = 4k, u = 0 on x = 0 and x = 1, natural elsewhere:  min E = -k/3 per unit cross-section")
    for N, bc2, lin_kind in (((8, 4), ('N', 'N'), 'cells'), ((16, 2), ('P', 'P'), 'nodes'), ((8, 2, 2), ('N', 'N'), 'cells')):
        n = len(N)
        bc = [('D', 'D')] + [bc2] * (n - 1)
        lin = (-4 * float(k)) * (np.ones(N) if lin_kind == 'cells' else np.ones([x + 1 for x in N]))
        P = S.Problem(box=[(0, 1)] * n, N=N, k=k, W=[0], bc=bc, lin=lin)
        cal = S.Calibration(P, 0, 0, 0)
        x = P.axes[0]
        shape = [-1] + [1] * (n - 1)
        cal.S0[0] = (2 * float(k) * (1 - 2 * x)).reshape(shape) * np.ones(P.face_shape(0))      # the exact flux 2k u'
        X = np.meshgrid(*P.axes, indexing='ij')
        U = X[0] * (1 - X[0])
        lb, parts = S.lower_bound(P, cal, return_parts=True)
        ub = S.upper_bound(P, U)[1]
        report(f"N={N} sides {bc2[0]} a on {lin_kind}  (rho {parts['rho_norm']:.1e})", lb, -k / 3, ub, tol_lb=1e-13)

    print("2. source f = c x:  min E = -c^2/(180 k) per unit cross-section")
    c = Fr(3, 2)
    for N in ((8, 2), (32, 2), (128, 2)):
        P = S.Problem(box=[(0, 1), (0, 1)], N=N, k=k, W=[0], bc=[('D', 'D'), ('N', 'N')], lin=lambda x, y: -float(c) * x)
        cal = S.Calibration(P, 0, 0, 0)
        x = P.axes[0]
        cal.S0[0] = (float(c) * (1 - 3 * x * x) / 6).reshape(-1, 1) * np.ones(P.face_shape(0))    # 2k u*' with 2k = 1
        X = np.meshgrid(*P.axes, indexing='ij')
        U = float(c) * X[0] * (1 - X[0] ** 2) / (12 * float(k))
        lb, parts = S.lower_bound(P, cal, return_parts=True)
        ub = S.upper_bound(P, U)[1]
        exact = -c * c / (180 * k)
        report(f"N={N}  (R1 {parts['R1']:.1e}, rho {parts['rho_norm']:.1e})", lb, exact, ub)
        print(f"      relative distance of the bounds from the minimum: lower {float((exact - Fr(lb)) / abs(exact)):.2e}, upper {float((Fr(ub) - exact) / abs(exact)):.2e}")

    print("3. convex potential with a prescribed mean, no Dirichlet face:  min E = |Omega| W(m)")
    for W, m, bcs in (([0, 0, 1], Fr(1, 4), [('N', 'N'), ('N', 'N')]), ([0, -1, Fr(1, 2), 0, Fr(1, 4)], Fr(1, 2), [('P', 'P'), ('N', 'N')])):
        P = S.Problem(box=[(0, 2), (0, 1)], N=(4, 4), k=k, W=W, bc=bcs, mean=m)
        cal = S.Calibration(P, 0, 0, 0)
        cal.lam = float(sum(j * Fr(cf) * m ** (j - 1) for j, cf in enumerate(W) if j >= 1))       # W'(m)
        U = float(m) * np.ones(P.fixed.shape)
        lb = S.lower_bound(P, cal)
        ub = S.upper_bound(P, U)[1]
        exact = 2 * sum(Fr(cf) * m ** j for j, cf in enumerate(W))
        report(f"W={[str(x) for x in W]} m={m} bc={bcs[0][0]}", lb, exact, ub, tol_lb=1e-6)

    print("4. source with a prescribed mean and no Dirichlet face (Proposition 16 b):  E = int k|u'|^2 - f u,")
    print("   f = cos(pi x) replaced by its nodal interpolant; the constant state is not the minimiser; here only LB <= UB")
    P = S.Problem(box=[(0, 1), (0, 1)], N=(16, 2), k=k, W=[0], bc=[('N', 'N'), ('N', 'N')], mean=Fr(1, 3),
                  lin=lambda x, y: -np.cos(math.pi * x))
    cal = S.Calibration(P, 0, 0, 0)
    x = P.axes[0]
    cal.S0[0] = (-np.sin(math.pi * x) / math.pi).reshape(-1, 1) * np.ones(P.face_shape(0))       # sigma' = -cos(pi x) = a
    cal.S0[0][0] = 0.0
    cal.S0[0][-1] = 0.0
    X = np.meshgrid(*P.axes, indexing='ij')
    U = 1 / 3 + np.cos(math.pi * X[0]) / (2 * float(k) * math.pi ** 2)                             # the minimiser for f = cos
    lb, parts = S.lower_bound(P, cal, return_parts=True)
    ub = S.upper_bound(P, U)[1]
    ok = lb <= ub
    fails += not ok
    print(f"  {lb!r} <= inf E <= {ub!r}   gap {100 * (ub - lb) / abs(ub):.2f} %   (rho {parts['rho_norm']:.1e}, R1 {parts['R1']:.1e})   {'ok' if ok else '**** CROSSED ****'}")

    print("5. competitor energy with entropy terms: the constant state u = m has energy |Omega| W(m)")
    chi = Fr(5, 2)
    pot = S.PolyLog([0, chi, -chi], 1, 1, (0, 1))
    for m in (Fr(1, 2), Fr(1, 8)):
        P = S.Problem(box=[(0, 1), (0, 1)], N=(4, 4), k=k, W=pot, bc=[('N', 'N'), ('N', 'N')], mean=m)
        lo, hi = S.upper_bound(P, float(m) * np.ones(P.fixed.shape))
        ref = float(m) * math.log(float(m)) + (1 - float(m)) * math.log(1 - float(m)) + float(chi * m * (1 - m))
        ok = lo <= ref + 1e-12 and ref - 1e-12 <= hi and hi - lo < 1e-10
        fails += not ok
        print(f"  m={m}: enclosure [{lo!r}, {hi!r}]   floating-point value {ref!r}   {'ok' if ok else '**** MISMATCH ****'}")

    print('ALL TESTS PASSED' if fails == 0 else f'FAILURES: {fails}')
    sys.exit(1 if fails else 0)
