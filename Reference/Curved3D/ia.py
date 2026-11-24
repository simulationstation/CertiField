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
    __array_ufunc__ = None

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


def sqr(a):
    """square of an interval or of a plain array"""
    return a.sqr() if isinstance(a, I) else a * a


def hi_of(a):
    return a.hi if isinstance(a, I) else a


def lo_of(a):
    return a.lo if isinstance(a, I) else a


def sum_up(a, axis=None):
    """upper bound of the exact sum of NON-NEGATIVE floats: any summation order errs by at most (n-1)u relative"""
    a = np.asarray(a, dtype=np.float64)
    if a.size and np.min(a) < 0:
        raise ValueError("sum_up needs non-negative terms")
    n = a.size if axis is None else a.shape[axis]
    s = np.sum(a, axis=axis)
    return np.nextafter(s * (1.0 + 4.0 * n * 2.0 ** -53), np.inf)


def sum_dn(a, axis=None):
    """lower bound of the exact sum of NON-NEGATIVE floats"""
    a = np.asarray(a, dtype=np.float64)
    if a.size and np.min(a) < 0:
        raise ValueError("sum_dn needs non-negative terms")
    n = a.size if axis is None else a.shape[axis]
    s = np.sum(a, axis=axis)
    return np.maximum(np.nextafter(s * (1.0 - 4.0 * n * 2.0 ** -53), -np.inf), 0.0)
