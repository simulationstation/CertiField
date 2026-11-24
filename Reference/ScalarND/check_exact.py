#!/usr/bin/env python3
"""Cross-check of the rigorous evaluators against an independent implementation in exact rational arithmetic.

The functions below recompute the certified lower-bound formula of NOTE.md (Proposition 2) directly from its
definition, cell by cell and vertex by vertex, with fractions.Fraction throughout: no floating point, no interval
library, no code shared with scalarnd.lower_bound.  For the same field S the interval-arithmetic value must not exceed
the exact one (outward rounding can only lose), and the two must agree to rounding accuracy.
Also checked the same way: the energy of the Q1 competitor (exact rational value inside the computed enclosure)."""
import itertools, sys
from fractions import Fraction as Fr
import numpy as np
import scalarnd as S

F = lambda x: Fr(float(x))


def poly(c, t):
    return sum(cf * t ** q for q, cf in enumerate(c))


def ipoly(c, lo, hi):
    """exact rational interval Horner"""
    a, b = c[-1], c[-1]
    for cf in reversed(c[:-1]):
        prods = [a * lo, a * hi, b * lo, b * hi]
        a, b = min(prods) + cf, max(prods) + cf
    return a, b


def exact_lower_calibration(P, cal, Qs, tails):
    """flux - sum_c |c| 2^-n sum_v M_c(v) in exact arithmetic.  Spline part only, constant datum on each face.
    tails = (inf of W on t <= tlo, inf of W on t >= thi), supplied exactly."""
    n, K = P.n, cal.K
    tlo, thi = cal.tlo, cal.thi
    Dl = (thi - tlo) / K
    k = P.k
    h = P.h

    def piece(t):                                             # knot interval and local coordinate of a clamped point
        z = min(max((t - tlo) / Dl, Fr(0)), Fr(K))
        kap = min(z.numerator // z.denominator, K - 1)
        return kap, z - kap

    def s_val(cf, t):
        kap, th = piece(t)
        return F(cf[kap]) * (1 - th) ** 2 / 2 + F(cf[kap + 1]) * (-2 * th * th + 2 * th + 1) / 2 + F(cf[kap + 2]) * th * th / 2

    flux = Fr(0)
    for i in range(n):
        for side in (0, 1):
            if P.bc[i][side] != 'D':
                continue
            sgn = -1 if side == 0 else 1
            area = P.cellvol / h[i]
            for f in itertools.product(*[range(P.N[j]) for j in range(n) if j != i]):
                node = list(f)
                node.insert(i, 0 if side == 0 else P.N[i])
                gvals = set()
                for ev in itertools.product((0, 1), repeat=n - 1):
                    nd = [node[j] for j in range(n)]
                    o = 0
                    for j in range(n):
                        if j != i:
                            nd[j] += ev[o]
                            o += 1
                    gvals.add(float(P.gval[tuple(nd)]))
                assert len(gvals) == 1, "exact check needs a constant datum on each face"
                flux += sgn * area * s_val(cal.C[i][tuple(node)], F(gvals.pop()))
    W, W2 = P.W, P.W2
    total = Fr(0)
    for c in itertools.product(*[range(Ni) for Ni in P.N]):
        fc = []                                               # coefficient vectors of the 2n faces of the cell
        for i in range(n):
            lo = cal.C[i][c]
            hi_idx = list(c)
            hi_idx[i] += 1
            fc.append((lo, cal.C[i][tuple(hi_idx)]))
        Mv = {v: None for v in P.verts}

        def push(v, val):
            if Mv[v] is None or val > Mv[v]:
                Mv[v] = val

        for kap in range(K):
            d2 = sum(((F(fc[i][1][kap]) - 2 * F(fc[i][1][kap + 1]) + F(fc[i][1][kap + 2]))
                      - (F(fc[i][0][kap]) - 2 * F(fc[i][0][kap + 1]) + F(fc[i][0][kap + 2]))) / Dl ** 2 / h[i] for i in range(n))
            prev = None
            for j in range(Qs + 1):
                th = Fr(j, Qs)
                t = tlo + (kap + th) * Dl
                sv, dv = [], []
                for i in range(n):
                    row_s, row_d = [], []
                    for side in (0, 1):
                        cf = fc[i][side]
                        c0, c1, c2 = F(cf[kap]), F(cf[kap + 1]), F(cf[kap + 2])
                        row_s.append(c0 * (1 - th) ** 2 / 2 + c1 * (-2 * th * th + 2 * th + 1) / 2 + c2 * th * th / 2)
                        row_d.append(((c1 - c0) * (1 - th) + (c2 - c1) * th) / Dl)
                    sv.append(row_s)
                    dv.append(row_d)
                D = sum((sv[i][1] - sv[i][0]) / h[i] for i in range(n))
                cur = {v: D + sum(dv[i][v[i]] ** 2 for i in range(n)) / (4 * k) - poly(W, t) for v in P.verts}
                if prev is not None:
                    ta = tlo + (kap + Fr(j - 1, Qs)) * Dl
                    corr = max(ipoly(W2, ta, t)[1] - d2, Fr(0)) * (Dl / Qs) ** 2 / 8
                    for v in P.verts:
                        push(v, max(prev[v], cur[v]) + corr)
                prev = cur
        for e, m0 in ((0, 0), (1, K)):
            Dend = sum(((F(fc[i][1][m0]) + F(fc[i][1][m0 + 1])) / 2 - (F(fc[i][0][m0]) + F(fc[i][0][m0 + 1])) / 2) / h[i] for i in range(n))
            for v in P.verts:
                push(v, Dend - tails[e])
        total += P.cellvol / 2 ** n * sum(Mv[v] for v in P.verts)
    return flux - total


def exact_lower_convex_quadratic(P, cal):
    """K = 0, W = c0 + c1 t + c2 t^2 (c2 > 0):  flux - sum_c [ int_c |sigma0|^2/(4k) + |c| W*(d_c) ],  W*(d) = (d-c1)^2/(4 c2) - c0"""
    n = P.n
    c0, c1, c2 = P.W
    k, h = P.k, P.h
    flux = Fr(0)
    for i in range(n):
        for side in (0, 1):
            if P.bc[i][side] != 'D':
                continue
            sgn = -1 if side == 0 else 1
            area = P.cellvol / h[i]
            for f in itertools.product(*[range(P.N[j]) for j in range(n) if j != i]):
                node = list(f)
                node.insert(i, 0 if side == 0 else P.N[i])
                gs = []
                for ev in itertools.product((0, 1), repeat=n - 1):
                    nd = list(node)
                    o = 0
                    for j in range(n):
                        if j != i:
                            nd[j] += ev[o]
                            o += 1
                    gs.append(F(P.gval[tuple(nd)]))
                flux += sgn * area * F(cal.S0[i][tuple(node)]) * sum(gs) / len(gs)
    total = Fr(0)
    for c in itertools.product(*[range(Ni) for Ni in P.N]):
        a0, d = Fr(0), Fr(0)
        for i in range(n):
            hi = list(c)
            hi[i] += 1
            lo_, hi_ = F(cal.S0[i][c]), F(cal.S0[i][tuple(hi)])
            a0 += (lo_ * lo_ + lo_ * hi_ + hi_ * hi_) / 3 / (4 * k)
            d += (hi_ - lo_) / h[i]
        total += P.cellvol * (a0 + (d - c1) ** 2 / (4 * c2) - c0)
    return flux - total


def exact_energy_quartic(P, U):
    """exact energy of the Q1 function for polynomial W: tensor Gauss-free, by expanding in monomials of the local
    coordinates and integrating term by term (n <= 3, small grids)."""
    n = P.n
    k, h = P.k, P.h
    verts = P.verts
    tot = Fr(0)
    for c in itertools.product(*[range(Ni) for Ni in P.N]):
        uv = {v: F(U[tuple(ci + vi for ci, vi in zip(c, v))]) for v in verts}
        # multilinear u = sum over subsets A of coefficient a_A prod_{i in A} xi_i  (local coordinates in [0,1])
        coef = {}
        for A in itertools.product((0, 1), repeat=n):
            s = Fr(0)
            for v in verts:
                if all(v[i] <= A[i] for i in range(n)):
                    s += uv[v] * (-1) ** (sum(A) - sum(v))
            coef[A] = s

        def mul(p, q):
            r = {}
            for ea, ca in p.items():
                for eb, cb in q.items():
                    e = tuple(x + y for x, y in zip(ea, eb))
                    r[e] = r.get(e, Fr(0)) + ca * cb
            return r

        def integ(p):
            return sum(cf / Fr(np.prod([e + 1 for e in ex])) if False else cf * Fr(1, int(np.prod([e + 1 for e in ex]))) for ex, cf in p.items())

        # potential
        powu = {tuple([0] * n): Fr(1)}
        acc = Fr(0)
        for q, wq in enumerate(P.W):
            if q > 0:
                powu = mul(powu, coef)
            acc += wq * integ(powu)
        # gradient
        for i in range(n):
            der = {}
            for A, cf in coef.items():
                if A[i] == 1:
                    e = list(A)
                    e[i] = 0
                    der[tuple(e)] = der.get(tuple(e), Fr(0)) + cf
            acc += k * integ(mul(der, der)) / h[i] ** 2
        tot += P.cellvol * acc
    return tot


if __name__ == '__main__':
    rng = np.random.default_rng(7)
    W4 = [Fr(1, 4), 0, Fr(-1, 2), 0, Fr(1, 4)]
    ok = True
    print("calibration mode: interval value vs exact rational value of the same formula")
    for (n, N, bc, box, K, Qs) in [(1, (9,), [('D', 'D')], [(-3, 3)], 5, 4),
                                   (2, (6, 3), [('D', 'D'), ('N', 'N')], [(-3, 3), (0, 1)], 4, 4),
                                   (2, (4, 4), [('D', 'D'), ('D', 'D')], [(-1, 1), (0, 2)], 3, 3),
                                   (3, (3, 2, 2), [('D', 'D'), ('N', 'N'), ('N', 'N')], [(-3, 3), (0, 1), (0, 1)], 3, 2),
                                   (3, (2, 2, 2), [('D', 'D'), ('D', 'N'), ('D', 'D')], [(0, 1), (0, 1), (0, 2)], 2, 2)]:
        # data that are constant on every Dirichlet face: a wall (-1 | +1) if only the x_1 faces are Dirichlet,
        # otherwise the constant 1/2
        if all(b == ('N', 'N') for b in bc[1:]):
            gg = lambda *X: np.where(X[0] < 0, -1.0, 1.0)
        else:
            gg = lambda *X: 0.5 + 0 * X[0]
        P = S.Problem(box=box, N=N, k=Fr(1, 2), W=W4, bc=bc, g=gg)
        cal = S.Calibration(P, -1, 1, K, affine=False)
        for trial in range(3):
            if trial == 0:
                S.optimise_calibration(P, cal, Qo=2, taus=(1e-2, 1e-3), iters=200)
            else:
                cal.unpack(cal.pack() + 0.3 * rng.standard_normal(cal.pack().size))
            lb_ia = S.lower_bound(P, cal, Qs=Qs)
            lb_ex = exact_lower_calibration(P, cal, Qs, (Fr(0), Fr(0)))       # inf of W beyond the wells is W(+-1) = 0
            good = Fr(lb_ia) <= lb_ex and float(lb_ex) - lb_ia < 1e-6 * (1 + abs(lb_ia))
            ok &= good
            print(f"  n={n} N={N} K={K} trial {trial}: interval {lb_ia:+.12f}   exact {float(lb_ex):+.12f}   exact - interval = {float(lb_ex - Fr(lb_ia)):.2e}  {'ok' if good else 'MISMATCH'}")
    print("convex mode (quadratic W): interval value vs exact rational value")
    for n, N in ((1, 7), (2, 5), (3, 3)):
        P = S.Problem(box=[(0, 1)] * n, N=(N,) * n, k=Fr(1, 2), W=[Fr(1, 8), -1, 1], bc=[('D', 'D')] * n, g=lambda *X: 0.25 * X[0])
        out, U, cal = S.certify_convex(P)
        for trial in range(3):
            if trial:
                cal.S0 = [s + 0.05 * rng.standard_normal(s.shape) for s in cal.S0]
            lb_ia = S.lower_bound(P, cal)
            lb_ex = exact_lower_convex_quadratic(P, cal)
            good = Fr(lb_ia) <= lb_ex and float(lb_ex) - lb_ia < 1e-6 * (1 + abs(lb_ia))
            ok &= good
            print(f"  n={n} N={N} trial {trial}: interval {lb_ia:+.12f}   exact {float(lb_ex):+.12f}   exact - interval = {float(lb_ex - Fr(lb_ia)):.2e}  {'ok' if good else 'MISMATCH'}")
    print("energy of the Q1 competitor: exact rational value inside the enclosure")
    for n, N in ((1, 6), (2, 4), (3, 2)):
        P = S.Problem(box=[(-1, 1)] * n, N=(N,) * n, k=Fr(1, 2), W=W4, bc=[('D', 'D')] * n, g=lambda *X: X[0])
        U = rng.standard_normal(P.fixed.shape)
        U[P.fixed] = P.gval[P.fixed]
        lo, hi = S.upper_bound(P, U)
        ex = exact_energy_quartic(P, U)
        good = Fr(lo) <= ex <= Fr(hi)
        ok &= good
        print(f"  n={n} N={N}: enclosure [{lo:.10f}, {hi:.10f}]   exact {float(ex):.10f}   {'ok' if good else 'MISMATCH'}")
    print('ALL EXACT CHECKS PASSED' if ok else 'EXACT CHECK FAILED')
    sys.exit(0 if ok else 1)
