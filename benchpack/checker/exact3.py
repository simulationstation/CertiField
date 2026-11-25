"""Exact trivariate tools for the 3D (P3) checker: polynomials on a cell in local coordinates (s, t, u) in [0,1]^3 as
NumPy object arrays of Fractions (P[a, b, c] = coefficient of s^a t^b u^c), products, derivatives, restriction to a
face, exact integrals of products over the cell or a dyadic sub-box, the Beltrami operator inc, and an exact test of
divisibility by a polynomial that is quadratic in one variable with a constant leading coefficient.
Nothing here is shared with the builder."""
import itertools
from fractions import Fraction as Fr
from functools import lru_cache
from math import lcm

import numpy as np

import exact as X

ZERO = Fr(0)


def zeros3(shape):
    out = np.empty(shape, dtype=object)
    out.fill(ZERO)
    return out


def t_add(*Ps):
    shp = tuple(max(P.shape[d] for P in Ps) for d in range(3))
    out = zeros3(shp)
    for P in Ps:
        out[:P.shape[0], :P.shape[1], :P.shape[2]] += P
    return out


def t_diff(P, ax):
    """d/d(local axis ax)"""
    n = P.shape[ax]
    if n == 1:
        shp = list(P.shape)
        return zeros3(tuple(shp))
    sl = [slice(None)] * 3
    sl[ax] = slice(1, None)
    w = np.arange(1, n, dtype=object)
    shape = [1, 1, 1]
    shape[ax] = n - 1
    return P[tuple(sl)] * w.reshape(shape)


def t_face(P, ax, which):
    """P restricted to the face local ax = 0 or 1: a 3D array with that axis of length 1"""
    if which == 0:
        sl = [slice(None)] * 3
        sl[ax] = slice(0, 1)
        return P[tuple(sl)].copy()
    return P.sum(axis=ax, keepdims=True)


def t_is_zero(P):
    return all(v == 0 for v in P.flat)


def t_equal(P, Q):
    return t_is_zero(t_add(P, -Q))


def to_int(*Ps):
    D = 1
    for P in Ps:
        for v in P.flat:
            D = lcm(D, v.denominator)
    out = []
    for P in Ps:
        Xi = np.empty(P.shape, dtype=object)
        for idx, v in np.ndenumerate(P):
            Xi[idx] = v.numerator * (D // v.denominator)
        out.append(Xi)
    return out, D


def t_mul(P, Q):
    """exact product (integers over a common denominator inside)"""
    (A, B), D = to_int(P, Q)
    out = np.zeros(tuple(P.shape[d] + Q.shape[d] - 1 for d in range(3)), dtype=object)
    for (a, b, c), v in np.ndenumerate(A):
        if v:
            out[a:a + Q.shape[0], b:b + Q.shape[1], c:c + Q.shape[2]] += v * B
    den = D * D
    res = np.empty(out.shape, dtype=object)
    for idx, v in np.ndenumerate(out):
        res[idx] = Fr(v, den)
    return res


@lru_cache(maxsize=None)
def moments(lo, hi, n):
    """exact int_lo^hi s^k ds, k < n (lo, hi Fractions)"""
    return tuple((hi ** (k + 1) - lo ** (k + 1)) / (k + 1) for k in range(n))


def t_integral(P, box=((ZERO, Fr(1)),) * 3):
    """exact integral of P over the local box [s0,s1] x [t0,t1] x [u0,u1]"""
    ms = [moments(box[d][0], box[d][1], P.shape[d]) for d in range(3)]
    tot = ZERO
    for (a, b, c), v in np.ndenumerate(P):
        if v:
            tot += v * ms[0][a] * ms[1][b] * ms[2][c]
    return tot


def quad_form3(polys, weights, box=((ZERO, Fr(1)),) * 3):
    """exact sum_{(i,j)} w_ij int_box P_i P_j  (integers inside: mode products with moment matrices)"""
    Xs, D = to_int(*polys)
    total = ZERO
    for (i, j), w in weights.items():
        if not w:
            continue
        A, B = Xs[i], Xs[j]
        Z = B
        dens = 1
        # contract axis by axis: Z[a, ...] = sum_c M[a, c] B[c, ...] with M[a, c] = int s^(a+c) over the box side
        for d in range(3):
            lo, hi = box[d]
            n = A.shape[d] + B.shape[d] - 1
            m = moments(lo, hi, n)
            L = 1
            for v in m:
                L = lcm(L, v.denominator)
            Mi = np.empty((A.shape[d], B.shape[d]), dtype=object)
            for a in range(A.shape[d]):
                for c in range(B.shape[d]):
                    v = m[a + c]
                    Mi[a, c] = v.numerator * (L // v.denominator)
            Z = np.moveaxis(np.tensordot(Mi, np.moveaxis(Z, d, 0), axes=([1], [0])), 0, d)
            dens *= L
        total += w * Fr(int((A * Z).sum()), dens)
    return total / (D * D)


def levi(i, j, k):
    return (i - j) * (j - k) * (k - i) // 2


def inc(Phi, scale):
    """tau_ij = sum eps_ikl eps_jmn d_k d_m Phi_ln for a symmetric potential Phi (dict (l, n) -> local array, both
    orders present); scale = (1/h_x, 1/h_y, 1/h_z) converts local to physical derivatives."""
    D2 = {}
    for (l, n), P in Phi.items():
        for k in range(3):
            dk = t_diff(P, k) * scale[k]
            for m in range(3):
                D2[(l, n, k, m)] = t_diff(dk, m) * scale[m]
    tau = {}
    for i in range(3):
        for j in range(i, 3):
            terms = []
            for k, l, m, n in itertools.product(range(3), repeat=4):
                s = levi(i, k, l) * levi(j, m, n)
                if s:
                    terms.append(D2[(l, n, k, m)] * s)
            tau[(i, j)] = t_add(*terms)
            tau[(j, i)] = tau[(i, j)]
    return tau


def divisible_quadratic(P, Q, ax):
    """exact: does Q divide P?  Q quadratic in local axis ax with a constant (nonzero) coefficient of ax^2 and no other
    terms of degree 2 in ax.  Division in ax over the other two variables; True iff the remainder vanishes."""
    Pm = np.moveaxis(P, ax, 0)
    Qm = np.moveaxis(Q, ax, 0)
    assert Qm.shape[0] == 3
    lead = Qm[2]
    q2 = lead[0, 0]
    assert q2 != 0 and all(v == 0 for v in lead.flat[1:])
    rows = [Pm[b].copy() for b in range(Pm.shape[0])]
    for b in range(len(rows) - 1, 1, -1):
        f = rows[b] / q2
        if all(v == 0 for v in f.flat):
            continue
        rows[b] = zeros3((1,) + f.shape)[0]
        rows[b - 1] = _add2(rows[b - 1], -_mul2(f, Qm[1]))
        rows[b - 2] = _add2(rows[b - 2], -_mul2(f, Qm[0]))
    return all(v == 0 for v in rows[0].flat) and (len(rows) < 2 or all(v == 0 for v in rows[1].flat))


def _add2(A, B):
    shp = (max(A.shape[0], B.shape[0]), max(A.shape[1], B.shape[1]))
    out = np.empty(shp, dtype=object)
    out.fill(ZERO)
    out[:A.shape[0], :A.shape[1]] += A
    out[:B.shape[0], :B.shape[1]] += B
    return out


def _mul2(A, B):
    out = np.empty((A.shape[0] + B.shape[0] - 1, A.shape[1] + B.shape[1] - 1), dtype=object)
    out.fill(ZERO)
    for (a, b), v in np.ndenumerate(A):
        if v:
            out[a:a + B.shape[0], b:b + B.shape[1]] += B * v
    return out


def local_from_dict(coefs, x0, h):
    """sum c_ijk x^i y^j z^k in local coordinates (x_d = x0[d] + h[d] s_d)"""
    deg = [max(key[d] for key in coefs) + 1 for d in range(3)]
    pw = [[X.ppow([x0[d], h[d]], e) for e in range(deg[d])] for d in range(3)]
    P = zeros3(tuple(deg))
    for (i, j, k), c in coefs.items():
        if not c:
            continue
        for a, u in enumerate(pw[0][i]):
            for b, v in enumerate(pw[1][j]):
                for cc, w in enumerate(pw[2][k]):
                    P[a, b, cc] += c * u * v * w
    return P


def tensor3(Mx, My, Mz, C):
    """sum C[i,j,k] Bx_i(s) By_j(t) Bz_k(u) in local power coefficients"""
    T = np.tensordot(Mx.T, C, axes=([1], [0]))            # (a, j, k)
    T = np.tensordot(T, My, axes=([1], [0]))              # (a, k, b)
    T = np.tensordot(T, Mz, axes=([1], [0]))              # (a, b, c)
    return T
