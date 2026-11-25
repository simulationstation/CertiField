"""Rigorous integrals of exact polynomials over curved parts of a grid cell (P1: the cells cut by the hole, and the
patch P_h), by slicing in x.

Region inside the cell [x0, x0+h] x [y0, y0+k] (with x >= 0):
    { x : max(lower functions) <= y <= min(upper functions) },
where every function y = f(x) is one of
    ('const', c)       y = c (c rational),
    ('circle', R2)     y = sqrt(R2 - x^2) on 0 <= x <= sqrt(R2) (R2 rational; undefined beyond: a lower function
                       is then inactive, an upper one makes the slice empty),
    ('line', name)     y = m x with m = cot(1/10) (the ray theta = pi/2 - 1/10 of the patch), an mpmath.iv enclosure.

Method.  All crossing abscissae of the functions (and the cell edges and the circles' end points) are listed; they
are rational, square roots of rationals (compared exactly by squaring), or transcendental (mpmath.iv enclosures,
which never coincide with algebraic numbers).  Between consecutive distinct break points no two functions cross, so
the active lower and upper functions are found at one rational test point.  On each piece
    int dx int_{L(x)}^{U(x)} p dy = k int [G(s, t_U(x)) - G(s, t_L(x))] dx,   G = int_0^t p dt,
and the x-integrals are done in closed form:
    const:  a polynomial in x with rational coefficients;
    circle: A(x) + B(x) sqrt(R2 - x^2), and int x^m sqrt(R2 - x^2) dx = P_m(x) sqrt(R2 - x^2) + c_m asin(x/R)
            (exact rational recursion; the identity is re-verified exactly for every use);
    line:   a polynomial in x with interval coefficients.
Only the final evaluation at the break points uses intervals: mpmath.iv (outward rounding) at PREC bits; asin is
enclosed by bisection-free verification with the interval sine (sin is increasing on [0, pi/2]).
"""
from fractions import Fraction as Fr
from functools import lru_cache

from mpmath import iv, mp
from mpmath.libmp import to_rational

import exact as X

PREC = 1200


def set_prec(bits=PREC):
    iv.prec = bits
    mp.prec = bits + 64


set_prec()


def ivq(q):
    """interval enclosure of a rational"""
    q = Fr(q)
    return iv.mpf(q.numerator) / iv.mpf(q.denominator)


def lo_q(I):
    return Fr(*to_rational(I._mpi_[0]))


def hi_q(I):
    return Fr(*to_rational(I._mpi_[1]))


def iv_lt(A, B):
    """True if A < B surely, False if A >= B surely, None if undecided"""
    if A.b < B.a:
        return True
    if A.a >= B.b:
        return False
    return None


def iv_asin(Z):
    """enclosure of asin(z) for z in Z subset [0, 1], verified with the interval sine"""
    zl, zh = lo_q(Z), hi_q(Z)
    assert zl >= 0 and zh <= 1, 'asin argument outside [0, 1]'
    P = iv.prec
    halfpi = mp.make_mpf((iv.pi / 2)._mpi_[1])          # >= pi/2, representable at P bits
    with mp.workprec(P + 64):                            # candidates only; every end point is verified below
        eps = mp.mpf(2) ** (-(P - 30))
        # lower end a (P bits): a <= 0, or sin(a) <= zl with a in [0, pi/2]; either way a <= asin(zl) <= asin(z)
        a = mp.asin(mp.mpf(zl.numerator) / zl.denominator) - eps if zl else mp.mpf(0)
        while True:
            a_r = mp.make_mpf(iv.mpf(a)._mpi_[0])
            if a_r <= 0 or iv.sin(iv.mpf(a_r)).b <= ivq(zl).a:
                break
            a -= eps
            eps *= 2
        # upper end b (P bits): b >= pi/2, or sin(b) >= zh with b in [0, pi/2]; either way b >= asin(zh)
        if zh == 1:
            b_r = halfpi
        else:
            b = mp.asin(mp.mpf(zh.numerator) / zh.denominator) + eps
            while True:
                b_r = mp.make_mpf(iv.mpf(b)._mpi_[1])
                if b_r >= halfpi or iv.sin(iv.mpf(b_r)).a >= ivq(zh).b:
                    break
                b += eps
                eps *= 2
            b_r = min(b_r, halfpi)
    return iv.mpf([a_r, b_r])


# ------------------------------------------------------------------------------------------------- numbers
class Num:
    """a non-negative break point: ('rat', q) = q, ('sqrt', q) = sqrt(q), ('tr', f) = transcendental (f() -> interval)"""
    __slots__ = ('kind', 'q', 'f', '_iv')

    def __init__(self, kind, q=None, f=None):
        self.kind, self.q, self.f, self._iv = kind, q, f, None

    def iv(self):
        if self._iv is None:
            if self.kind == 'rat':
                self._iv = ivq(self.q)
            elif self.kind == 'sqrt':
                self._iv = iv.sqrt(ivq(self.q))
            else:
                self._iv = self.f()
        return self._iv

    def sq(self):
        """x^2 as an exact rational (algebraic kinds) or None"""
        if self.kind == 'rat':
            return self.q * self.q
        if self.kind == 'sqrt':
            return self.q
        return None


def num_cmp(a, b):
    """-1, 0, 1 exactly for algebraic numbers; by enclosures for transcendental ones (never equal to algebraic)"""
    sa, sb = a.sq(), b.sq()
    if sa is not None and sb is not None:          # both >= 0: compare squares
        return (sa > sb) - (sa < sb)
    A, B = a.iv(), b.iv()
    r = iv_lt(A, B)
    if r is True:
        return -1
    if iv_lt(B, A) is True:
        return 1
    raise ArithmeticError('break points not separated at %d bits' % iv.prec)


@lru_cache(maxsize=None)
def cot_tenth():
    t = ivq(Fr(1, 10))
    return iv.cos(t) / iv.sin(t)


def line_slope(name):
    assert name == 'cot(1/10)'
    return cot_tenth()


# ------------------------------------------------------------------------------------------------- closed forms
@lru_cache(maxsize=None)
def sqrt_moments(R2, m):
    """(P_m, c_m): int x^m sqrt(R2 - x^2) dx = P_m(x) sqrt(R2 - x^2) + c_m asin(x/R) + const"""
    if m == 0:
        return [Fr(0), Fr(1, 2)], R2 / 2
    if m == 1:
        return [-R2 / 3, Fr(0), Fr(1, 3)], Fr(0)
    Pm2, cm2 = sqrt_moments(R2, m - 2)
    head = [Fr(0)] * (m - 1) + [-R2 / (m + 2), Fr(0), Fr(1, m + 2)]      # -x^{m-1}(R2 - x^2)/(m+2)
    f = (m - 1) * R2 / (m + 2)
    return X.padd(head, X.pscale(Pm2, f)), f * cm2


def sqrt_antider(O, R2):
    """P, c with d/dx [P(x) sqrt(R2-x^2) + c asin(x/R)] = O(x) sqrt(R2-x^2), verified exactly"""
    P, c = [], Fr(0)
    for m, om in enumerate(O):
        if om:
            Pm, cm = sqrt_moments(R2, m)
            P = X.padd(P, X.pscale(Pm, om))
            c += om * cm
    # exact check:  P'(R2 - x^2) - x P + c == O (R2 - x^2)
    w = [R2, Fr(0), Fr(-1)]
    lhs = X.padd(X.padd(X.pmul(X.pder(P), w), X.pscale([Fr(0)] + P, -1)), [c])
    assert X.same_poly(lhs, X.pmul(O, w)), 'closed form of int x^m sqrt(R2 - x^2) failed'
    return P, c


def iv_poly(a, Xv):
    """Horner with interval arithmetic; coefficients exact rationals or intervals"""
    r = iv.mpf(0)
    for c in reversed(a):
        r = r * Xv + (ivq(c) if isinstance(c, (Fr, int)) else c)
    return r


# ------------------------------------------------------------------------------------------------- region integral
def _fn_value_at(fn, xq):
    """interval value of the function at the rational point xq, or None if undefined there"""
    kind, par = fn
    if kind == 'const':
        return ivq(par)
    if kind == 'circle':
        if xq * xq > par:
            return None
        return iv.sqrt(ivq(par - xq * xq))
    return line_slope(par) * ivq(xq)


def _crossings(f, g):
    """abscissae x >= 0 where f = g (list of Num)"""
    (kf, pf), (kg, pg) = f, g
    if kf == kg:
        return []
    if kf > kg:
        (kf, pf), (kg, pg) = (kg, pg), (kf, pf)
    # kinds sorted: 'circle' < 'const' < 'line'
    if kf == 'circle' and kg == 'const':
        if pg >= 0 and pg * pg <= pf:
            return [Num('sqrt', pf - pg * pg)]
        return []
    if kf == 'const' and kg == 'line':
        if pf >= 0:
            return [Num('tr', f=lambda c=pf, n=pg: ivq(c) / line_slope(n))]
        return []
    if kf == 'circle' and kg == 'line':
        return [Num('tr', f=lambda R2=pf, n=pg: iv.sqrt(ivq(R2)) / iv.sqrt(1 + line_slope(n) ** 2))]
    return []


def _subs_G(Gx, fn, y0, k):
    """the univariate integrand x -> G(s(x), t_f(x)), as ('poly', coeffs) or ('circ', E, O, R2) or ('ivpoly', coeffs);
    Gx[b] = coefficient polynomial (in x) of t^b."""
    kind, par = fn
    if kind == 'const':
        tc = (par - y0) / k
        out = []
        p = Fr(1)
        for b in range(len(Gx)):
            out = X.padd(out, X.pscale(Gx[b], p))
            p *= tc
        return ('poly', out)
    if kind == 'circle':
        R2 = par
        nb = len(Gx)
        # t = (sigma - y0)/k ; t^b = sum_j C(b,j) sigma^j (-y0)^(b-j) / k^b
        W = [[] for _ in range(nb)]
        for b in range(nb):
            if not X.ptrim(Gx[b]):
                continue
            kb = k ** b
            for j in range(b + 1):
                coef = X.binom(b, j) * (-y0) ** (b - j) / kb
                if coef:
                    W[j] = X.padd(W[j], X.pscale(Gx[b], coef))
        w = [R2, Fr(0), Fr(-1)]
        E, O = [], []
        wp = [Fr(1)]
        for j in range(nb):
            if j >= 2 and j % 2 == 0:
                wp = X.pmul(wp, w)
            if W[j]:
                if j % 2 == 0:
                    E = X.padd(E, X.pmul(W[j], wp))
                else:
                    O = X.padd(O, X.pmul(W[j], wp))
        return ('circ', E, O, R2)
    # line y = m x
    m = line_slope(par)
    lin = [ivq(-y0 / k), m / ivq(k)]            # t = (m x - y0)/k
    out = []
    tp = [iv.mpf(1)]
    for b in range(len(Gx)):
        if X.ptrim(Gx[b]):
            term = _ivpmul([ivq(c) for c in Gx[b]], tp)
            out = _ivpadd(out, term)
        tp = _ivpmul(tp, lin)
    return ('ivpoly', out)


def _ivpmul(a, b):
    if not a or not b:
        return []
    out = [iv.mpf(0)] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            out[i + j] += u * v
    return out


def _ivpadd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def _antider_eval(item, xa, xb):
    """interval value of int_{xa}^{xb} of the univariate integrand `item` (xa, xb are Num)"""
    if item[0] == 'poly':
        A = X.pint(item[1])
        return iv_poly(A, xb.iv()) - iv_poly(A, xa.iv())
    if item[0] == 'ivpoly':
        A = [iv.mpf(0)] + [c / (i + 1) for i, c in enumerate(item[1])]
        return iv_poly(A, xb.iv()) - iv_poly(A, xa.iv())
    _, E, O, R2 = item
    val = iv.mpf(0)
    if E:
        A = X.pint(E)
        val += iv_poly(A, xb.iv()) - iv_poly(A, xa.iv())
    if O:
        P, c = sqrt_antider(O, R2)
        R = iv.sqrt(ivq(R2))

        def at(xn):
            sq = xn.sq()
            if sq is not None:
                sig = iv.sqrt(ivq(R2 - sq)) if R2 - sq > 0 else iv.mpf(0)
                assert R2 - sq >= 0
            else:
                Xv = xn.iv()
                sig = iv.sqrt(ivq(R2) - Xv * Xv)
            v = iv_poly(P, xn.iv()) * sig
            if c:
                v += ivq(c) * iv_asin(_clip01(xn.iv() / R))
            return v
        val += at(xb) - at(xa)
    return val


def _clip01(Z):
    lo, hi = Z.a, Z.b
    if lo < 0:
        lo = iv.mpf(0)
    if hi > 1:
        hi = iv.mpf(1)
    return iv.mpf([lo.a, hi.b])


def integrate_region(P, x0, h, y0, k, lowers, uppers):
    """interval enclosure of  int over {x in [x0, x0+h], max(y0, lowers) <= y <= min(y0+k, uppers)} of p dx dy,
    p given in the cell's local coordinates by the Fraction array P (P[a, b] s^a t^b).  Requires x0 >= 0."""
    assert x0 >= 0
    x1, y1 = x0 + h, y0 + k
    lows = [('const', y0)] + list(lowers)
    ups = [('const', y1)] + list(uppers)
    fns = lows + ups
    # break points
    bps = [Num('rat', x0), Num('rat', x1)]
    for i in range(len(fns)):
        if fns[i][0] == 'circle':
            bps.append(Num('sqrt', fns[i][1]))
        for j in range(i + 1, len(fns)):
            bps += _crossings(fns[i], fns[j])
    lo, hi = bps[0], bps[1]
    inside = [b for b in bps if num_cmp(b, lo) >= 0 and num_cmp(b, hi) <= 0]
    # sort and merge equal ones
    import functools
    inside.sort(key=functools.cmp_to_key(num_cmp))
    pts = []
    for b in inside:
        if not pts or num_cmp(pts[-1], b) != 0:
            pts.append(b)
    # G = int_0^t P dt, as coefficient polynomials in global x of each t^b
    n1, n2 = P.shape
    G = X.zeros((n1, n2 + 1))
    for b in range(n2):
        G[:, b + 1] = P[:, b] / (b + 1)
    Gx = [X.pcompose_lin(list(G[:, b]), -x0 / h, 1 / h) for b in range(n2 + 1)]
    total = iv.mpf(0)
    subs_cache = {}
    for a, b in zip(pts[:-1], pts[1:]):
        A, B = a.iv(), b.iv()
        assert A.b < B.a, 'pieces not separated'
        xq = (hi_q(A) + lo_q(B)) / 2
        # active lower / upper at the test point
        best_lo, best_lo_v = None, None
        for fn in lows:
            v = _fn_value_at(fn, xq)
            if v is None:
                continue
            if best_lo is None or iv_lt(best_lo_v, v):
                best_lo, best_lo_v = fn, v
            elif not iv_lt(v, best_lo_v):
                raise ArithmeticError('lower functions not separated at the test point')
        best_up, best_up_v, empty = None, None, False
        for fn in ups:
            v = _fn_value_at(fn, xq)
            if v is None:
                empty = True
                break
            if best_up is None or iv_lt(v, best_up_v):
                best_up, best_up_v = fn, v
            elif not iv_lt(best_up_v, v):
                raise ArithmeticError('upper functions not separated at the test point')
        if empty:
            continue
        r = iv_lt(best_lo_v, best_up_v)
        if r is None:
            raise ArithmeticError('emptiness undecided at the test point')
        if not r:
            continue
        for fn in (best_up, best_lo):
            if fn not in subs_cache:
                subs_cache[fn] = _subs_G(Gx, fn, y0, k)
        total += ivq(k) * (_antider_eval(subs_cache[best_up], a, b) - _antider_eval(subs_cache[best_lo], a, b))
    return total
