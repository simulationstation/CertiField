#!/usr/bin/env python3
"""Exact cross-check of the three-dimensional elasticity reference on tiny grids, in rational arithmetic
(fractions.Fraction), with no floating point and no code shared with the evaluators of elasticity3d.py.

For stress fields and candidates in the format of elasticity3d.py (computed ones and random ones), the functions
below rebuild both fields cell by cell as polynomials in three variables and check, exactly:

  1. static admissibility (the hypotheses of Theorem 1 of docs/ELASTICITY_3D.md):
       div sigma + b = 0 in every cell, the traction is the same from both sides of every interior cell face, and
       it vanishes identically on the traction-free faces;
  2. kinematic admissibility: the candidate is continuous;
  3. the exact values LB and E of the two bounds lie inside the computed enclosures:
       interval lower bound <= LB,   lower end <= E <= upper end,   and they agree to rounding;
  4. the identity behind the theorem,   E - LB = 1/2 int (e - S s)^T D (e - S s)   (so E - LB >= 0), which holds
     only if the integration by parts of the proof is right for these fields.

The stress formulas below are typed in from the note, not taken from the table that elasticity3d.py generates.
The only data taken from elasticity3d.py are the problem, the nodal values, and the exact values that it prescribes
at the nodes of the traction-free faces; that these make the traction vanish is what check 1 verifies.
"""
import itertools
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import elasticity3d as E

F = lambda x: Fr(float(x))


# ---- polynomials in three variables: dict {(i, j, k): coefficient} (global coordinates) ------------------------------
def padd(p, q, c=1):
    r = dict(p)
    for k, v in q.items():
        nv = r.get(k, Fr(0)) + c * v
        if nv == 0:
            r.pop(k, None)
        else:
            r[k] = nv
    return r


def pder(p, var):
    r = {}
    for e, v in p.items():
        if e[var]:
            ne = list(e)
            ne[var] -= 1
            r[tuple(ne)] = r.get(tuple(ne), Fr(0)) + e[var] * v
    return {k: v for k, v in r.items() if v != 0}


def prestrict(p, var, value):
    """the polynomial on the plane {x_var = value}: the exponent of x_var is set to zero"""
    r = {}
    for e, v in p.items():
        ne = list(e)
        ne[var] = 0
        r[tuple(ne)] = r.get(tuple(ne), Fr(0)) + v * value ** e[var]
    return {k: v for k, v in r.items() if v != 0}


def moments(lo, hi, nmax=16):
    return [(hi ** (e + 1) - lo ** (e + 1)) / (e + 1) for e in range(nmax)]


def pint_prod(p, q, mom, skip=None):
    """exact integral of p q over a box (or over a face: skip = the normal direction, whose exponents must be zero)"""
    tot = Fr(0)
    dirs = [d for d in range(3) if d != skip]
    for a, v in p.items():
        for c, w in q.items():
            t = v * w
            for d in dirs:
                t *= mom[d][a[d] + c[d]]
            tot += t
    return tot


def pint(p, mom):
    return sum(v * mom[0][e[0]] * mom[1][e[1]] * mom[2][e[2]] for e, v in p.items())


def p1mul(p, q):
    r = {}
    for a, v in p.items():
        for c, w in q.items():
            r[a + c] = r.get(a + c, Fr(0)) + v * w
    return r


def poly1(coeffs, x0, h):
    """the one-variable polynomial sum c_k ((x - x0)/h)^k as {power: coefficient} of x"""
    r = {}
    base = {0: Fr(1)}
    lin = {0: -x0 / h, 1: 1 / h}
    for c in coeffs:
        for k, v in base.items():
            r[k] = r.get(k, Fr(0)) + c * v
        base = p1mul(base, lin)
    return r


def tensor(px, py, pz, c=1):
    out = {}
    for a, v in px.items():
        for b_, w in py.items():
            vw = v * w * c
            if vw == 0:
                continue
            for d, y in pz.items():
                if y != 0:
                    out[(a, b_, d)] = vw * y
    return out


HERM = [[1, 0, -3, 2], [0, 1, -2, 1], [0, 0, 3, -2], [0, 0, -1, 1]]
LAGR = {1: [[1, -1], [0, 1]], 2: [[1, -3, 2], [0, 4, -4], [0, -1, 2]]}
NAMES = ['11', '22', '33', '23', '13', '12']


def d2(p, i, j):
    return pder(pder(p, i), j)


def beltrami(Phi):
    """the six stress components (order 11, 22, 33, 23, 13, 12) from the formulas of docs/ELASTICITY_3D.md, Section 3.2"""
    P = lambda nm: Phi.get(nm, {})
    s11 = padd(padd(d2(P('33'), 1, 1), d2(P('22'), 2, 2)), d2(P('23'), 1, 2), -2)
    s22 = padd(padd(d2(P('11'), 2, 2), d2(P('33'), 0, 0)), d2(P('13'), 0, 2), -2)
    s33 = padd(padd(d2(P('22'), 0, 0), d2(P('11'), 1, 1)), d2(P('12'), 0, 1), -2)
    s23 = padd(padd(padd(d2(P('12'), 0, 2), d2(P('13'), 0, 1)), d2(P('11'), 1, 2), -1), d2(P('23'), 0, 0), -1)
    s13 = padd(padd(padd(d2(P('23'), 0, 1), d2(P('12'), 1, 2)), d2(P('22'), 0, 2), -1), d2(P('13'), 1, 1), -1)
    s12 = padd(padd(padd(d2(P('13'), 1, 2), d2(P('23'), 0, 2)), d2(P('33'), 0, 1), -1), d2(P('12'), 2, 2), -1)
    return [s11, s22, s33, s23, s13, s12]


VOIGT = [[0, 5, 4], [5, 1, 3], [4, 3, 2]]


def exact_nodal_values(P, st):
    """X[f][t][i][j][k] as Fractions: stored floats, and at the prescribed entries the exact values"""
    mask, Tb, T = P.prescribed()
    X = np.empty(st.X.shape, dtype=object)
    for idx in np.ndindex(st.X.shape):
        X[idx] = F(st.X[idx])
    for q, idx in enumerate(np.argwhere(mask)):
        X[tuple(idx)] = Tb[q] + sum(F(st.par[k]) * T[k][q] for k in range(len(T)))
    return X


def exact_fields(P, U, st):
    """per cell: the stress (six polynomials) and the displacement (three polynomials)"""
    X = exact_nodal_values(P, st)
    h, p = P.h, P.p
    cells = {}
    for c in itertools.product(*[range(n) for n in P.N]):
        x0 = [P.vnodes[d][c[d]] for d in range(3)]
        H1 = [[poly1([Fr(v) for v in HERM[k]], x0[d], h[d]) for k in range(4)] for d in range(3)]
        Phi = {}
        for fi, f in enumerate(P.funcs):
            acc = {}
            for k in itertools.product(range(4), repeat=3):
                t = (k[0] % 2) + 2 * (k[1] % 2) + 4 * (k[2] % 2)
                val = X[fi, t, c[0] + k[0] // 2, c[1] + k[1] // 2, c[2] + k[2] // 2]
                if val != 0:
                    val = val * h[0] ** (k[0] % 2) * h[1] ** (k[1] % 2) * h[2] ** (k[2] % 2)
                    acc = padd(acc, tensor(H1[0][k[0]], H1[1][k[1]], H1[2][k[2]], val))
            Phi[NAMES[f]] = acc
        s = beltrami(Phi)
        s[0] = padd(s[0], {(1, 0, 0): -P.b[0]} if P.b[0] != 0 else {})
        s[1] = padd(s[1], {(0, 1, 0): -P.b[1]} if P.b[1] != 0 else {})
        s[2] = padd(s[2], {(0, 0, 1): -P.b[2]} if P.b[2] != 0 else {})
        L1 = [[poly1([Fr(v) for v in LAGR[p][q]], x0[d], h[d]) for q in range(p + 1)] for d in range(3)]
        u = []
        for a in range(3):
            ua = {}
            for q in itertools.product(range(p + 1), repeat=3):
                val = F(U[a, p * c[0] + q[0], p * c[1] + q[1], p * c[2] + q[2]])
                if val != 0:
                    ua = padd(ua, tensor(L1[0][q[0]], L1[1][q[1]], L1[2][q[2]], val))
            u.append(ua)
        cells[c] = dict(s=s, u=u, lo=x0, hi=[x0[d] + h[d] for d in range(3)])
    return cells


def check(P, U, st, tag):
    t0 = time.time()
    cells = exact_fields(P, U, st)
    S, D, b = P.S, P.D, P.b
    eq = cont = free = True
    en = comp = defect = bnd = Fr(0)
    for c, d in cells.items():
        s, u = d['s'], d['u']
        # 1a. equilibrium
        for i in range(3):
            r = {(0, 0, 0): b[i]} if b[i] != 0 else {}
            for j in range(3):
                r = padd(r, pder(s[VOIGT[i][j]], j))
            eq &= (r == {})
        mom = [moments(d['lo'][k], d['hi'][k]) for k in range(3)]
        for i in range(3):
            # 1b, 2. the neighbour across the upper face in direction i
            nb = tuple(c[k] + (1 if k == i else 0) for k in range(3))
            if nb in cells:
                for j in range(3):
                    cont &= (prestrict(s[VOIGT[i][j]], i, d['hi'][i]) == prestrict(cells[nb]['s'][VOIGT[i][j]], i, d['hi'][i]))
                    cont &= (prestrict(u[j], i, d['hi'][i]) == prestrict(cells[nb]['u'][j], i, d['hi'][i]))
            # 1c and the boundary term
            for side, val, onb in ((0, d['lo'][i], c[i] == 0), (1, d['hi'][i], c[i] == P.N[i] - 1)):
                if not onb:
                    continue
                if P.bc[i][side] == 'N':
                    for j in range(3):
                        free &= (prestrict(s[VOIGT[i][j]], i, val) == {})
                else:
                    sign = -1 if side == 0 else 1
                    for j in range(3):                 # on Gamma_D the trace of the candidate is the datum
                        bnd += sign * pint_prod(prestrict(s[VOIGT[i][j]], i, val), prestrict(u[j], i, val), mom, skip=i)
        # 3. the two functionals and the defect
        e = [pder(u[0], 0), pder(u[1], 1), pder(u[2], 2), padd(pder(u[1], 2), pder(u[2], 1)), padd(pder(u[0], 2), pder(u[2], 0)),
             padd(pder(u[0], 1), pder(u[1], 0))]
        r = [dict(x) for x in e]
        for m in range(6):
            for n in range(6):
                if S[m][n] != 0:
                    r[m] = padd(r[m], s[n], -S[m][n])
        for m in range(6):
            for n in range(6):
                if D[m][n] != 0:
                    en += D[m][n] * pint_prod(e[m], e[n], mom) / 2
                    defect += D[m][n] * pint_prod(r[m], r[n], mom) / 2
                if S[m][n] != 0:
                    comp += S[m][n] * pint_prod(s[m], s[n], mom) / 2
        en -= sum(b[a] * pint(u[a], mom) for a in range(3) if b[a] != 0)
    LB = bnd - comp
    lb_ia = E.lower_bound(P, st)
    ub_ia = E.upper_bound(P, U)
    ident = (en - LB == defect)
    good = (eq and cont and free and Fr(lb_ia) <= LB and Fr(ub_ia[0]) <= en <= Fr(ub_ia[1]) and ident and defect >= 0
            and float(LB) - lb_ia < 1e-8 * (1 + abs(lb_ia)) and ub_ia[1] - ub_ia[0] < 1e-8 * (1 + abs(ub_ia[1])))
    print(f"  {tag:40s} equilibrium {'exact' if eq else 'VIOLATED'}  tractions/continuity {'exact' if cont else 'VIOLATED'}  "
          f"free faces {'exact' if free else 'VIOLATED'}  LB - interval = {float(LB - Fr(lb_ia)):.1e}  "
          f"E in enclosure {Fr(ub_ia[0]) <= en <= Fr(ub_ia[1])}  identity {'exact' if ident else 'VIOLATED'}  "
          f"E - LB = {float(en - LB):.3e}  {'ok' if good else 'MISMATCH'}  ({time.time() - t0:.0f} s)", flush=True)
    return good


def anisotropic():
    """a symmetric positive definite 6x6 stiffness with all couplings"""
    lam, mu = Fr(3, 2), Fr(1)
    D = [[(lam + 2 * mu if i == j else lam) if (i < 3 and j < 3) else (mu if i == j else Fr(0)) for j in range(6)] for i in range(6)]
    pert = {(0, 3): Fr(1, 4), (0, 5): Fr(-1, 8), (1, 4): Fr(1, 5), (2, 5): Fr(1, 3), (3, 4): Fr(1, 8), (4, 5): Fr(-1, 6), (1, 3): Fr(1, 7),
            (0, 1): Fr(1, 4), (3, 5): Fr(1, 9)}
    for (i, j), v in pert.items():
        D[i][j] += v
        D[j][i] += v
    return D


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding="utf-8")
    quick = '--quick' in sys.argv
    rng = np.random.default_rng(7)
    ok = True
    lam, mu = Fr(3, 2), Fr(1)
    gq = lambda X, Y, Z: (0.1 * X - 0.15 * Y + 0.3 * X * Z, 0.25 * X + 0.2 * Y * Z - 0.1, 0.05 * Z - 0.2 * X * Y + 0.1 * Y)
    DD, DN, ND, NN = ('D', 'D'), ('D', 'N'), ('N', 'D'), ('N', 'N')
    cases = [
        ('all Dirichlet, p=1', dict(box=[(0, 2), (0, 1), (0, 1)], N=(2, 1, 1), lam=lam, mu=mu, degree=1, g=gq)),
        ('all Dirichlet, p=2, body force', dict(box=[(0, 2), (0, 1), (Fr(1, 3), 1)], N=(1, 2, 1), lam=lam, mu=mu, degree=2, g=gq,
                                                 b=(Fr(1, 2), Fr(-1, 3), Fr(1, 5)))),
        ('ring (four free faces), p=1', dict(box=[(0, 1)] * 3, N=(2, 1, 1), lam=lam, mu=mu, degree=1, bc=(DD, NN, NN), g=gq)),
        ('ring, body force, anisotropic, p=2', dict(box=[(-1, 1), (Fr(1, 2), 2), (0, 1)], N=(1, 1, 2), D=anisotropic(), degree=2, bc=(NN, DD, NN), g=gq,
                                                     b=(Fr(1, 3), Fr(1, 4), Fr(-1, 2)))),
        ('one free face, anisotropic', dict(box=[(0, 1), (0, 2), (0, 1)], N=(1, 2, 1), D=anisotropic(), degree=1, bc=(DD, DN, DD), g=gq, b=(Fr(0), Fr(-1), Fr(0)))),
        ('two opposite free faces, body force', dict(box=[(Fr(1, 3), 1), (0, 1), (0, 1)], N=(2, 1, 1), lam=lam, mu=mu, degree=2, bc=(NN, DD, DD), g=gq,
                                                      b=(Fr(1, 5), Fr(-2, 3), Fr(1, 2)))),
        ('three free faces at a corner', dict(box=[(0, 1)] * 3, N=(1, 1, 2), lam=lam, mu=mu, degree=1, bc=(ND, ND, ND), g=gq, b=(Fr(1), Fr(1), Fr(1)))),
        ('five free faces, body force', dict(box=[(0, 1), (0, 1), (0, 2)], N=(1, 1, 2), lam=lam, mu=mu, degree=2, bc=(NN, NN, DN), g=gq, b=(Fr(0), Fr(0), Fr(-1)))),
        ('ring, all states kept', dict(box=[(0, 1)] * 3, N=(1, 2, 1), lam=lam, mu=mu, degree=1, bc=(DD, NN, NN), g=gq, states='all', b=(Fr(1, 2), Fr(0), Fr(1, 4)))),
        ('Maxwell functions only, ring', dict(box=[(0, 1)] * 3, N=(1, 1, 2), lam=lam, mu=mu, degree=1, bc=(DD, NN, NN), g=gq, functions='maxwell', states='all')),
        ('Morera functions only, one free face', dict(box=[(0, 1)] * 3, N=(2, 1, 1), lam=lam, mu=mu, degree=1, bc=(DD, DD, DN), g=gq, functions='morera')),
        ('2 x 2 x 2 cells, two adjacent free faces', dict(box=[(0, 1), (0, 1), (-1, 1)], N=(2, 2, 2), lam=lam, mu=mu, degree=1, bc=(DN, DN, DD), g=gq,
                                                           b=(Fr(1, 2), Fr(1, 3), Fr(0)))),
    ]
    if quick:
        cases = cases[:4]
    for tag, kw in cases:
        P = E.Problem3D(**kw)
        out, U, st = E.certify(P)
        ok &= check(P, U, st, tag)
        # the same with random fields: every nodal array and every parameter vector is admissible
        U2 = U + 0.3 * rng.standard_normal(U.shape)
        U2[:, P.fixed] = P.gval[:, P.fixed]
        st2 = E.Stress(P, st.X + 0.3 * rng.standard_normal(st.X.shape), st.par + 0.3 * rng.standard_normal(st.par.shape))
        ok &= check(P, U2, st2, '   random fields')
    print('ALL EXACT CHECKS PASSED' if ok else 'EXACT CHECK FAILED')
    sys.exit(0 if ok else 1)
