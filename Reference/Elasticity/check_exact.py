#!/usr/bin/env python3
"""Exact cross-check of the elasticity reference on small grids, in rational arithmetic (fractions.Fraction), with no
floating point and no code shared with the evaluators of elasticity.py.

For the stress field and the candidate that elasticity.py produces, the functions below rebuild both fields cell by
cell as polynomials in two variables and check, exactly:

  1. static admissibility (the hypotheses of Theorem 1 of docs/ELASTICITY.md):
       div sigma + b = 0 in every cell, the normal traction is the same from both sides of every interior cell
       face, and it vanishes identically on the traction-free faces;
  2. kinematic admissibility: the candidate is continuous and equals the Q_p interpolant of the datum on Gamma_D;
  3. the exact values LB and E of the two bounds lie inside the computed enclosures:
       interval lower bound <= LB,   lower end <= E <= upper end,   and they agree to rounding;
  4. the identity behind the theorem,   E - LB = 1/2 int (e - S s)^T D (e - S s)   (so E - LB >= 0), which holds
     only if the integration by parts of the proof is right for these fields.
"""
import itertools
import sys
from fractions import Fraction as Fr

import numpy as np

import elasticity as E

F = lambda x: Fr(float(x))


# ---- polynomials in two variables: dict {(i, j): coefficient} of x^i y^j (global coordinates) --------------------------
def padd(p, q, c=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, Fr(0)) + c * v
    return {k: v for k, v in r.items() if v != 0}


def pmul(p, q):
    r = {}
    for (a, b), v in p.items():
        for (c, d), w in q.items():
            r[(a + c, b + d)] = r.get((a + c, b + d), Fr(0)) + v * w
    return {k: v for k, v in r.items() if v != 0}


def pscale(p, c):
    return {k: v * c for k, v in p.items() if v * c != 0}


def pder(p, var):
    r = {}
    for (a, b), v in p.items():
        if var == 0 and a:
            r[(a - 1, b)] = r.get((a - 1, b), Fr(0)) + a * v
        if var == 1 and b:
            r[(a, b - 1)] = r.get((a, b - 1), Fr(0)) + b * v
    return r


def pint(p, x0, x1, y0, y1):
    return sum(v * (x1 ** (a + 1) - x0 ** (a + 1)) / (a + 1) * (y1 ** (b + 1) - y0 ** (b + 1)) / (b + 1) for (a, b), v in p.items())


def prestrict(p, var, value):
    """the polynomial of the other variable on the line {x_var = value}, as {power: coefficient}"""
    r = {}
    for (a, b), v in p.items():
        k, c = (b, v * value ** a) if var == 0 else (a, v * value ** b)
        r[k] = r.get(k, Fr(0)) + c
    return {k: v for k, v in r.items() if v != 0}


def p1int(p, t0, t1):
    return sum(v * (t1 ** (k + 1) - t0 ** (k + 1)) / (k + 1) for k, v in p.items())


def p1mul(p, q):
    r = {}
    for a, v in p.items():
        for c, w in q.items():
            r[a + c] = r.get(a + c, Fr(0)) + v * w
    return {k: v for k, v in r.items() if v != 0}


def poly1(coeffs, x0, h):
    """the one-variable polynomial sum c_k ((x - x0)/h)^k as {(power,): coefficient} of x"""
    r = {0: Fr(0)}
    base = {0: Fr(1)}
    lin = {0: -x0 / h, 1: 1 / h}
    for c in coeffs:
        for k, v in base.items():
            r[k] = r.get(k, Fr(0)) + c * v
        base = p1mul(base, lin)
    return r


def tensor(px, py):
    return {(a, b): v * w for a, v in px.items() for b, w in py.items() if v * w != 0}


HERM = [[1, 0, -3, 2], [0, 1, -2, 1], [0, 0, 3, -2], [0, 0, -1, 1]]
LAGR = {1: [[1, -1], [0, 1]], 2: [[1, -3, 2], [0, 4, -4], [0, -1, 2]]}


def exact_fields(P, U, st):
    """per cell: the Airy function, the stress (s11, s22, s12) and the displacement (u1, u2) as exact polynomials"""
    nodes, Tb, T = P.traces()
    X = [[[F(st.X[t, i, j]) for j in range(P.N[1] + 1)] for i in range(P.N[0] + 1)] for t in range(4)]
    for t, i, j in E._pinned(P):
        X[t][i][j] = Fr(0)
    for q, (i, j) in enumerate(nodes):
        for t in range(4):
            X[t][i][j] = Tb[q][t] + sum(F(st.par[k]) * T[k][q][t] for k in range(len(T)))
    h, p = P.h, P.p
    cells = {}
    for c in itertools.product(range(P.N[0]), range(P.N[1])):
        x0, y0 = P.vnodes[0][c[0]], P.vnodes[1][c[1]]
        phi = {}
        for k1 in range(4):
            for k2 in range(4):
                e1, d1, e2, d2 = k1 // 2, k1 % 2, k2 // 2, k2 % 2
                val = X[d1 + 2 * d2][c[0] + e1][c[1] + e2] * h[0] ** d1 * h[1] ** d2
                if val != 0:
                    phi = padd(phi, tensor(poly1([Fr(v) for v in HERM[k1]], x0, h[0]), poly1([Fr(v) for v in HERM[k2]], y0, h[1])), val)
        s11 = padd(pder(pder(phi, 1), 1), {(1, 0): -P.b[0]})
        s22 = padd(pder(pder(phi, 0), 0), {(0, 1): -P.b[1]})
        s12 = pscale(pder(pder(phi, 0), 1), Fr(-1))
        u = []
        for a in range(2):
            ua = {}
            for i in range(p + 1):
                for j in range(p + 1):
                    val = F(U[a, p * c[0] + i, p * c[1] + j])
                    if val != 0:
                        ua = padd(ua, tensor(poly1([Fr(v) for v in LAGR[p][i]], x0, h[0]), poly1([Fr(v) for v in LAGR[p][j]], y0, h[1])), val)
            u.append(ua)
        cells[c] = dict(s=(s11, s22, s12), u=u, box=(x0, x0 + h[0], y0, y0 + h[1]))
    return cells


def traction(s, i):
    """the traction on a face with normal +e_i: (s11, s12) or (s12, s22)"""
    return (s[0], s[2]) if i == 0 else (s[2], s[1])


def check(P, U, st, tag):
    cells = exact_fields(P, U, st)
    ok = True
    N1, N2 = P.N
    # 1. static admissibility
    for c, d in cells.items():
        s = d['s']
        r1 = padd(padd(pder(s[0], 0), pder(s[2], 1)), {(0, 0): P.b[0]})
        r2 = padd(padd(pder(s[2], 0), pder(s[1], 1)), {(0, 0): P.b[1]})
        ok &= (r1 == {} and r2 == {})
    eq = ok
    cont = True
    for c, d in cells.items():
        x0, x1, y0, y1 = d['box']
        for i, nb, val in ((0, (c[0] + 1, c[1]), x1), (1, (c[0], c[1] + 1), y1)):
            if nb in cells:
                for ta, tb in zip(traction(d['s'], i), traction(cells[nb]['s'], i)):
                    cont &= (prestrict(ta, i, val) == prestrict(tb, i, val))
                for ua, ub in zip(d['u'], cells[nb]['u']):
                    cont &= (prestrict(ua, i, val) == prestrict(ub, i, val))
    free = True
    for c, d in cells.items():
        x0, x1, y0, y1 = d['box']
        for i, side, val, onb in ((0, 0, x0, c[0] == 0), (0, 1, x1, c[0] == N1 - 1), (1, 0, y0, c[1] == 0), (1, 1, y1, c[1] == N2 - 1)):
            if onb and P.bc[i][side] == 'N':
                for ta in traction(d['s'], i):
                    free &= (prestrict(ta, i, val) == {})
    # 3. the exact bounds
    S, D = P.S, P.D
    en = Fr(0)
    comp = Fr(0)
    defect = Fr(0)
    for c, d in cells.items():
        x0, x1, y0, y1 = d['box']
        u, s = d['u'], d['s']
        e = [pder(u[0], 0), pder(u[1], 1), padd(pder(u[0], 1), pder(u[1], 0))]
        for m in range(3):
            for n in range(3):
                if D[m][n] != 0:
                    en += D[m][n] * pint(pmul(e[m], e[n]), x0, x1, y0, y1) / 2
                if S[m][n] != 0:
                    comp += S[m][n] * pint(pmul(s[m], s[n]), x0, x1, y0, y1) / 2
        en -= P.b[0] * pint(u[0], x0, x1, y0, y1) + P.b[1] * pint(u[1], x0, x1, y0, y1)
        # e - S s
        r = [padd(e[m], {}, 1) for m in range(3)]
        for m in range(3):
            for n in range(3):
                if S[m][n] != 0:
                    r[m] = padd(r[m], s[n], -S[m][n])
        for m in range(3):
            for n in range(3):
                if D[m][n] != 0:
                    defect += D[m][n] * pint(pmul(r[m], r[n]), x0, x1, y0, y1) / 2
    bnd = Fr(0)
    for c, d in cells.items():
        x0, x1, y0, y1 = d['box']
        for i, side, val, onb in ((0, 0, x0, c[0] == 0), (0, 1, x1, c[0] == N1 - 1), (1, 0, y0, c[1] == 0), (1, 1, y1, c[1] == N2 - 1)):
            if onb and P.bc[i][side] == 'D':
                sign = -1 if side == 0 else 1
                t0, t1 = (y0, y1) if i == 0 else (x0, x1)
                tr = traction(d['s'], i)
                # on Gamma_D the trace of the candidate is the datum
                bnd += sign * sum(p1int(p1mul(prestrict(tr[a], i, val), prestrict(d['u'][a], i, val)), t0, t1) for a in range(2))
    LB = bnd - comp
    lb_ia = E.lower_bound(P, st)
    ub_ia = E.upper_bound(P, U)
    good = (eq and cont and free and Fr(lb_ia) <= LB and Fr(ub_ia[0]) <= en <= Fr(ub_ia[1]) and en - LB == defect and defect >= 0
            and float(LB) - lb_ia < 1e-9 * (1 + abs(lb_ia)) and ub_ia[1] - ub_ia[0] < 1e-9 * (1 + abs(ub_ia[1])))
    print(f"  {tag:34s} equilibrium {'exact' if eq else 'VIOLATED'}  tractions/continuity {'exact' if cont else 'VIOLATED'}  "
          f"free faces {'exact' if free else 'VIOLATED'}  LB - interval = {float(LB - Fr(lb_ia)):.1e}  "
          f"E in enclosure {Fr(ub_ia[0]) <= en <= Fr(ub_ia[1])}  identity {'exact' if en - LB == defect else 'VIOLATED'}  "
          f"gap {float(en - LB):.3e}  {'ok' if good else 'MISMATCH'}")
    return good


if __name__ == '__main__':
    rng = np.random.default_rng(5)
    ok = True
    lam, mu = Fr(3, 2), Fr(1)
    Daniso = [[Fr(4), Fr(1), Fr(1, 2)], [Fr(1), Fr(3), Fr(-1, 4)], [Fr(1, 2), Fr(-1, 4), Fr(2)]]
    g1 = lambda X, Y: (0.1 * X - 0.15 * Y + 0.3 * X * Y, 0.25 * X + 0.2 * Y * Y - 0.1)
    cases = [
        ('all Dirichlet, p=1', dict(box=[(0, 2), (0, 1)], N=(3, 2), lam=lam, mu=mu, degree=1, g=g1)),
        ('all Dirichlet, p=2, body force', dict(box=[(0, 2), (0, 1)], N=(2, 2), lam=lam, mu=mu, degree=2, g=g1, b=(Fr(1, 2), Fr(-1, 3)))),
        ('free top and bottom, p=1', dict(box=[(0, 1), (0, 1)], N=(3, 3), lam=lam, mu=mu, degree=1, bc=(('D', 'D'), ('N', 'N')), g=lambda X, Y: (0.1 * X, 0.05 * X * Y))),
        ('free top/bottom, body force, p=2', dict(box=[(-1, 1), (Fr(1, 2), 2)], N=(2, 3), lam=lam, mu=mu, degree=2, bc=(('D', 'D'), ('N', 'N')), g=g1, b=(Fr(1, 3), Fr(1, 4)))),
        ('one free face, anisotropic', dict(box=[(0, 1), (0, 2)], N=(2, 3), D=Daniso, degree=1, bc=(('D', 'D'), ('D', 'N')), g=g1, b=(Fr(0), Fr(-1)))),
        ('three free faces, body force', dict(box=[(Fr(1, 3), 1), (0, 1)], N=(3, 2), D=Daniso, degree=2, bc=(('D', 'N'), ('N', 'N')), g=g1, b=(Fr(1, 5), Fr(-2, 3)))),
        ('two adjacent free faces', dict(box=[(0, 1), (0, 1)], N=(2, 2), lam=lam, mu=mu, degree=1, bc=(('N', 'D'), ('N', 'D')), g=g1, b=(Fr(1), Fr(1)))),
    ]
    for tag, kw in cases:
        P = E.Problem(**kw)
        out, U, st = E.certify(P)
        ok &= check(P, U, st, tag)
        # the same with random fields: every nodal array and every parameter vector is admissible
        U2 = U + 0.3 * rng.standard_normal(U.shape)
        U2[:, P.fixed] = P.gval[:, P.fixed]
        st2 = E.Stress(P, st.X + 0.3 * rng.standard_normal(st.X.shape), st.par + 0.3 * rng.standard_normal(st.par.shape))
        ok &= check(P, U2, st2, '   random fields')
    print('ALL EXACT CHECKS PASSED' if ok else 'EXACT CHECK FAILED')
    sys.exit(0 if ok else 1)
