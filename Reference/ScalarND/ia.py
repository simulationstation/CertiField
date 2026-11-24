#!/usr/bin/env python3
"""Minimal rigorous interval arithmetic on numpy float64 arrays.

Every + - * result is widened by one unit in the last place in each direction (numpy.nextafter).  Under IEEE-754
round-to-nearest the exact result of one operation lies within half an ulp of the computed one, so the widened
interval contains it.  Exact rational constants (fractions.Fraction) are converted by exact comparison.
No library elementary function is used.
"""
import math
from fractions import Fraction
import numpy as np

INF = np.inf


def dn(x):
    return np.nextafter(x, -INF)


def up(x):
    return np.nextafter(x, INF)


class I:
    """closed interval [lo, hi]; lo and hi are float64 arrays that broadcast against each other"""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = np.asarray(lo, dtype=np.float64)
        self.hi = self.lo if hi is None else np.asarray(hi, dtype=np.float64)

    def __add__(self, o):
        o = _c(o)
        return I(dn(self.lo + o.lo), up(self.hi + o.hi))

    __radd__ = __add__

    def __sub__(self, o):
        o = _c(o)
        return I(dn(self.lo - o.hi), up(self.hi - o.lo))

    def __rsub__(self, o):
        return _c(o) - self

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __mul__(self, o):
        o = _c(o)
        p1 = self.lo * o.lo
        p2 = self.lo * o.hi
        p3 = self.hi * o.lo
        p4 = self.hi * o.hi
        return I(dn(np.minimum(np.minimum(p1, p2), np.minimum(p3, p4))),
                 up(np.maximum(np.maximum(p1, p2), np.maximum(p3, p4))))

    __rmul__ = __mul__

    def sqr(self):
        l2 = self.lo * self.lo
        h2 = self.hi * self.hi
        mx = np.maximum(l2, h2)
        mn = np.where((self.lo <= 0) & (self.hi >= 0), 0.0, np.minimum(l2, h2))
        return I(np.maximum(dn(mn), 0.0), up(mx))

    def inv(self):
        if not np.all(self.lo > 0):
            raise ZeroDivisionError("inv: interval not strictly positive")
        return I(dn(1.0 / self.hi), up(1.0 / self.lo))

    def hull(self, o):
        o = _c(o)
        return I(np.minimum(self.lo, o.lo), np.maximum(self.hi, o.hi))

    def pos(self):
        """enclosure of max(x, 0)"""
        return I(np.maximum(self.lo, 0.0), np.maximum(self.hi, 0.0))

    def __getitem__(self, k):
        return I(self.lo[k], self.hi[k])

    @property
    def shape(self):
        return np.broadcast(self.lo, self.hi).shape


def _c(o):
    """coerce: I stays; Fraction/int become exact enclosures; float is taken as the exact binary number it is"""
    if isinstance(o, I):
        return o
    if isinstance(o, Fraction):
        return fr(o)
    if isinstance(o, (int, np.integer)):
        return fr(Fraction(int(o)))
    return I(o)


def frac_enclosure(q):
    """floats lo <= q <= hi for an exact Fraction, verified by exact comparison"""
    x = float(q)
    lo = hi = x
    while Fraction(lo) > q:
        lo = math.nextafter(lo, -math.inf)
    while Fraction(hi) < q:
        hi = math.nextafter(hi, math.inf)
    return lo, hi


def fr(q):
    lo, hi = frac_enclosure(Fraction(q))
    return I(lo, hi)


def horner(coefs, t):
    """enclosure of sum_k coefs[k] t^k; coefs are I / Fraction / float (ascending), t is I"""
    p = _c(coefs[-1])
    for c in reversed(coefs[:-1]):
        p = p * t + _c(c)
    return p


def poly_deriv(coefs):
    """derivative of an ascending coefficient list of Fractions"""
    return [k * coefs[k] for k in range(1, len(coefs))] or [Fraction(0)]


# ================================================================================================================
# additions for docs/SCALAR_ND_2.md: logarithm, verified square root, a lower bound for pi.  Still only + - * / of
# doubles with outward rounding, comparisons and exact rational arithmetic; no library elementary function is trusted.
# ================================================================================================================
def _ln2_enclosure(terms=40):
    """ln 2 = 2 artanh(1/3) = 2 sum_j 3^-(2j+1)/(2j+1): the partial sum is a lower bound, and the remainder after the
    term j = terms is below 2 * 3^-(2 terms + 3) / ((2 terms + 3)(1 - 1/9))  (geometric tail)"""
    z = Fraction(1, 3)
    s = sum(z ** (2 * j + 1) / (2 * j + 1) for j in range(terms + 1))
    rem = 2 * z ** (2 * terms + 3) / ((2 * terms + 3) * (1 - z * z))
    return frac_enclosure(2 * s)[0], frac_enclosure(2 * s + rem)[1]


_LN2 = _ln2_enclosure()
_LN_N = 12                                                   # terms j = 0 .. _LN_N of the artanh series
_LN_Z = Fraction(9, 50)                                      # |z| <= 9/50 after the reduction below
_LN_REM = frac_enclosure(2 * _LN_Z ** (2 * _LN_N + 3) / ((2 * _LN_N + 3) * (1 - _LN_Z ** 2)))[1]


def _ln_point(x):
    """enclosure (I) of ln at the positive doubles x (Lemma 12 of docs/SCALAR_ND_2.md):
    x = mu 2^e exactly with mu in [0.70, 1.42], ln x = e ln 2 + 2 artanh((mu - 1)/(mu + 1))"""
    x = np.asarray(x, dtype=np.float64)
    mu, e = np.frexp(x)                                       # x = mu 2^e, 1/2 <= mu < 1, exact
    small = mu < 0.7071067811865476
    mu = np.where(small, 2.0 * mu, mu)                        # exact
    e = np.where(small, e - 1, e)
    assert np.all(np.ldexp(mu, e) == x) and np.all(mu >= 0.70) and np.all(mu <= 1.42)
    m = I(mu)
    z = (m - I(1.0)) * (m + I(1.0)).inv()
    assert np.all(z.lo >= -0.1799) and np.all(z.hi <= 0.1799)   # so |z| <= 9/50, as the remainder bound assumes
    y = z.sqr()
    acc = fr(Fraction(1, 2 * _LN_N + 1))
    for j in range(_LN_N - 1, -1, -1):                        # Horner in y = z^2 for sum_j y^j / (2j + 1)
        acc = acc * y + fr(Fraction(1, 2 * j + 1))
    s = z * acc * I(2.0)
    s = I(dn(s.lo - _LN_REM), up(s.hi + _LN_REM))
    return I(e.astype(np.float64)) * I(_LN2[0], _LN2[1]) + s


def ln(x):
    """RIGOROUS enclosure of the natural logarithm of a positive interval (ln is increasing: end by end)"""
    x = _c(x)
    if not (np.all(x.lo > 0) and np.all(np.isfinite(x.hi))):
        raise ValueError("ln: the interval must be positive and finite")
    return I(_ln_point(x.lo).lo, _ln_point(x.hi).hi)


def xlnx(x):
    """RIGOROUS enclosure of s ln s for an interval s >= 0 (value 0 at s = 0).  Where the lower end is positive this is
    the product of the two enclosures.  Where it is 0:  s ln s <= max(0, b ln b) on [0, b], and s ln s >= b ln b if
    b <= 1/4 (the function decreases on [0, 1/e]), >= -1 otherwise (ln s >= 1 - 1/s gives s ln s >= s - 1)."""
    x = _c(x)
    lo = np.asarray(x.lo, dtype=np.float64) * np.ones(x.shape)
    hi = np.asarray(x.hi, dtype=np.float64) * np.ones(x.shape)
    if not (np.all(lo >= 0) and np.all(np.isfinite(hi))):
        raise ValueError("xlnx: the interval must be non-negative and finite")
    pos = lo > 0
    lo1 = np.where(pos, lo, 1.0)
    hi1 = np.where(hi > 0, hi, 1.0)
    full = I(lo1, np.maximum(hi1, lo1)) * ln(I(lo1, np.maximum(hi1, lo1)))        # used where pos
    b = I(hi1) * ln(I(hi1))                                                       # b ln b at the upper end
    zlo = np.where(hi <= 0.25, np.minimum(b.lo, 0.0), -1.0)
    zhi = np.maximum(b.hi, 0.0)
    out_lo = np.where(pos, full.lo, np.where(hi > 0, zlo, 0.0))
    out_hi = np.where(pos, full.hi, np.where(hi > 0, zhi, 0.0))
    return I(out_lo, out_hi)


def sqrt_up(x):
    """a double s >= sqrt(x) for a double x >= 0, verified with one multiplication: the product s*s rounded to nearest
    and moved one unit down is a lower bound for the exact s^2, and it is checked to be >= x.  The library square root
    only supplies the first guess."""
    x = float(x)
    if not (x >= 0.0 and math.isfinite(x)):
        raise ValueError("sqrt_up: needs a finite non-negative number")
    if x == 0.0:
        return 0.0
    if x < 1e-300:                                   # near the subnormal range the verification below cannot succeed
        return 1e-150                                # in a reasonable number of steps; (1e-150)^2 = 1e-300 > x
    s = float(up(math.sqrt(x)))
    while not float(dn(s * s)) >= x:
        s = float(up(s))
    return s


def pi_lower(terms=12):
    """an exact rational lower bound for pi from the formula of Machin, pi = 16 arctan(1/5) - 4 arctan(1/239).
    arctan(1/q) = sum_j (-1)^j q^-(2j+1)/(2j+1) is alternating with decreasing terms, so a partial sum that ends with a
    negative term is a lower bound and one that ends with a positive term is an upper bound."""
    def atan_partial(q, last):
        return sum(Fraction((-1) ** j, (2 * j + 1) * q ** (2 * j + 1)) for j in range(last + 1))
    lo5 = atan_partial(5, 2 * terms + 1)                     # ends with a negative term
    hi239 = atan_partial(239, 2 * terms)                     # ends with a positive term
    return 16 * lo5 - 4 * hi239
