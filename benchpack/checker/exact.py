"""Exact rational tools of the independent checker: univariate and bivariate polynomials, B-spline pieces, and exact
integrals of products of polynomials over the unit square.

Everything here is exact (fractions.Fraction, or Python integers with an explicit common denominator).  Nothing in this
file is shared with the certificate builder.

Conventions:
- a univariate polynomial is a list of Fractions, entry a = coefficient of s^a;
- a bivariate polynomial on a cell is a NumPy object array P of Fractions, P[a, b] = coefficient of s^a t^b, where
  s = (x - x0)/h and t = (y - y0)/k are the local coordinates of the cell [x0, x0 + h] x [y0, y0 + k].
"""
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb, lcm

import numpy as np

ZERO = Fr(0)
ONE = Fr(1)


def F(v):
    """exact Fraction from a float (its binary value), an int, a Fraction or a 'p/q' string"""
    if isinstance(v, Fr):
        return v
    if isinstance(v, str):
        return Fr(v)
    if isinstance(v, (int, np.integer)):
        return Fr(int(v))
    return Fr(float(v))     # float64 -> exact binary rational


def farr(a):
    """object array of exact Fractions from an array of float64 (or anything F accepts)"""
    a = np.asarray(a)
    out = np.empty(a.shape, dtype=object)
    for idx, v in np.ndenumerate(a):
        out[idx] = F(v)
    return out


def zeros(shape):
    out = np.empty(shape, dtype=object)
    out.fill(ZERO)
    return out


# ------------------------------------------------------------------------------------------------- univariate
def ptrim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else ZERO) + (b[i] if i < len(b) else ZERO) for i in range(n)]


def pscale(a, c):
    return [c * v for v in a]


def pmul(a, b):
    if not a or not b:
        return []
    out = [ZERO] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        if u:
            for j, v in enumerate(b):
                out[i + j] += u * v
    return out


def ppow(a, n):
    out = [ONE]
    for _ in range(n):
        out = pmul(out, a)
    return out


def pder(a):
    return [i * a[i] for i in range(1, len(a))]


def pint(a):
    """antiderivative vanishing at 0"""
    return [ZERO] + [a[i] / (i + 1) for i in range(len(a))]


def peval(a, x):
    r = ZERO * x if not isinstance(x, Fr) else ZERO
    for c in reversed(a):
        r = r * x + c
    return r


def pcompose_lin(a, c0, c1):
    """a(c0 + c1 x) as a polynomial in x"""
    out = []
    lin = [c0, c1]
    for c in reversed(a):
        out = padd(pmul(out, lin), [c])
    return out


# ------------------------------------------------------------------------------------------------- B-splines
def bspline_pieces(knots, p):
    """Local polynomial pieces of the B-splines of degree p on the knot vector `knots` (list of Fractions).

    Returns a list of spans (k, t_k, t_{k+1}, M), one per non-empty knot span [t_k, t_{k+1}), where M is a
    (p+1) x (p+1) object array with M[i, a] = coefficient of s^a of B_{k-p+i}(t_k + (t_{k+1} - t_k) s).
    Cox-de Boor recursion in exact arithmetic (terms with a zero denominator omitted)."""
    t = knots
    n = len(t) - p - 1
    spans = []
    for k in range(p, n):
        if not t[k] < t[k + 1]:
            continue
        h = t[k + 1] - t[k]
        B = {k: [ONE]}
        for d in range(1, p + 1):
            NB = {}
            for i in range(k - d, k + 1):
                acc = [ZERO] * (d + 1)
                if i in B and t[i + d] != t[i]:
                    den = t[i + d] - t[i]
                    acc = padd(acc, pmul([(t[k] - t[i]) / den, h / den], B[i]))
                if (i + 1) in B and t[i + d + 1] != t[i + 1]:
                    den = t[i + d + 1] - t[i + 1]
                    acc = padd(acc, pmul([(t[i + d + 1] - t[k]) / den, -h / den], B[i + 1]))
                NB[i] = acc
            B = NB
        M = zeros((p + 1, p + 1))
        for i in range(p + 1):
            row = B[k - p + i]
            for a in range(p + 1):
                M[i, a] = row[a] if a < len(row) else ZERO
        # partition of unity on the span, exactly
        tot = M.sum(axis=0)
        assert tot[0] == 1 and all(v == 0 for v in tot[1:]), 'partition of unity fails'
        spans.append((k, t[k], t[k + 1], M))
    return spans


def knot_report(knots, p, lo, hi):
    """checks of an open knot vector on [lo, hi]; returns (ok, max interior multiplicity, message)"""
    t = knots
    msg = []
    ok = True
    if any(t[i] > t[i + 1] for i in range(len(t) - 1)):
        ok = False
        msg.append('not non-decreasing')
    if t[0] != lo or t[-1] != hi:
        ok = False
        msg.append('does not span [%s, %s]' % (lo, hi))
    if not (all(v == t[0] for v in t[:p + 1]) and all(v == t[-1] for v in t[-p - 1:])):
        ok = False
        msg.append('not open')
    mult = {}
    for v in t:
        mult[v] = mult.get(v, 0) + 1
    inner = [m for v, m in mult.items() if v != t[0] and v != t[-1]]
    return ok, (max(inner) if inner else 0), '; '.join(msg)


# ------------------------------------------------------------------------------------------------- bivariate
def tensor(Mx, C, My):
    """sum_{i,j} C[i,j] Bx_i(s) By_j(t) as local coefficients: Mx^T C My"""
    return Mx.T.dot(C).dot(My)


def b_add(*Ps):
    n1 = max(P.shape[0] for P in Ps)
    n2 = max(P.shape[1] for P in Ps)
    out = zeros((n1, n2))
    for P in Ps:
        out[:P.shape[0], :P.shape[1]] += P
    return out


def b_scale(P, c):
    return P * c


def b_conv_s(P, a):
    """P(s,t) * a(s)"""
    a = ptrim(a)
    if not a:
        return zeros((1, P.shape[1]))
    out = zeros((P.shape[0] + len(a) - 1, P.shape[1]))
    for j, c in enumerate(a):
        if c:
            out[j:j + P.shape[0], :] += P * c
    return out


def b_conv_t(P, b):
    return b_conv_s(P.T, b).T


def b_ds(P):
    if P.shape[0] == 1:
        return zeros((1, P.shape[1]))
    return P[1:, :] * np.arange(1, P.shape[0], dtype=object)[:, None]


def b_dt(P):
    return b_ds(P.T).T


def b_at_s(P, which):
    """P at s = 0 or s = 1, as a univariate polynomial in t (list)"""
    return list(P[0, :]) if which == 0 else list(P.sum(axis=0))


def b_at_t(P, which):
    return list(P[:, 0]) if which == 0 else list(P.sum(axis=1))


def b_mul(P, Q):
    """exact product of two bivariate polynomials (integer arithmetic inside)"""
    (X, Y), D = to_int(P, Q)
    out = np.zeros((P.shape[0] + Q.shape[0] - 1, P.shape[1] + Q.shape[1] - 1), dtype=object)
    for (a, b), v in np.ndenumerate(X):
        if v:
            out[a:a + Q.shape[0], b:b + Q.shape[1]] += v * Y
    den = D * D
    res = np.empty(out.shape, dtype=object)
    for idx, v in np.ndenumerate(out):
        res[idx] = Fr(v, den)
    return res


def same_poly(a, b):
    return ptrim(a) == ptrim(b)


def is_zero_poly(a):
    return not ptrim(a)


# ------------------------------------------------------------------------------------------------- exact integrals
def to_int(*Ps):
    """common denominator D and integer arrays X_i = D * P_i"""
    D = 1
    for P in Ps:
        for v in P.flat:
            D = lcm(D, v.denominator)
    out = []
    for P in Ps:
        X = np.empty(P.shape, dtype=object)
        for idx, v in np.ndenumerate(P):
            X[idx] = v.numerator * (D // v.denominator)
        out.append(X)
    return out, D


@lru_cache(maxsize=None)
def hankel(n1, n2):
    """integer matrix L/(a+c+1), a < n1, c < n2, and L = lcm(1..n1+n2-1)"""
    L = 1
    for m in range(1, n1 + n2):
        L = lcm(L, m)
    H = np.empty((n1, n2), dtype=object)
    for a in range(n1):
        for c in range(n2):
            H[a, c] = L // (a + c + 1)
    return H, L


def quad_form(polys, weights):
    """exact value of  sum_{(i,j) in weights} w_ij * int_0^1 int_0^1 P_i P_j ds dt  (P_i bivariate Fraction arrays)

    Uses int s^a t^b = 1/((a+1)(b+1)) in integer arithmetic: int P Q = sum(X * (Hs Y Ht)) / (D^2 Ls Lt)."""
    Xs, D = to_int(*polys)
    total = Fr(0)
    cache = {}
    for (i, j), w in weights.items():
        if not w:
            continue
        X, Y = Xs[i], Xs[j]
        key = (j, X.shape)
        if key not in cache:
            Hs, Ls = hankel(X.shape[0], Y.shape[0])
            Ht, Lt = hankel(Y.shape[1], X.shape[1])
            cache[key] = (Hs.dot(Y).dot(Ht), Ls * Lt)
        Z, L = cache[key]
        val = int((X * Z).sum())
        total += w * Fr(val, L)
    return total / (D * D)


def integral(P):
    """exact int_0^1 int_0^1 P ds dt"""
    tot = Fr(0)
    for (a, b), v in np.ndenumerate(P):
        if v:
            tot += v / ((a + 1) * (b + 1))
    return tot


def integral_rect(P, s0, s1, t0, t1):
    """exact int_{s0}^{s1} int_{t0}^{t1} P ds dt (rational limits)"""
    tot = Fr(0)
    for (a, b), v in np.ndenumerate(P):
        if v:
            tot += v * (s1 ** (a + 1) - s0 ** (a + 1)) / (a + 1) * (t1 ** (b + 1) - t0 ** (b + 1)) / (b + 1)
    return tot


def binom(n, k):
    return comb(n, k)
