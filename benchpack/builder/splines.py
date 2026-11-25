"""One-dimensional B-splines with exact rational polynomial pieces, float evaluation and rigorous Gauss tables.

A knot vector t (Fractions, open) and a degree p define n = len(t) - p - 1 B-splines (Cox-de Boor).  On every
non-empty span [t_k, t_{k+1}] (a "cell" e) the p+1 non-zero B-splines B_{k-p}, ..., B_k are polynomials; they are
stored exactly as coefficient lists in the local coordinate s in [0, 1], x = a_e + h_e s.
"""
import os
import sys
from fractions import Fraction as Fr

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', '..', 'Reference', 'Curved3D'))
from ia import I, frac_enclosure, fr          # noqa: E402  (outward-rounded interval arithmetic of curved3d)
from curved3d import gauss_enclosures, RI, ri_horner     # noqa: E402  (rational Gauss enclosures of curved3d)


def poly_mul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r


def poly_add(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else Fr(0)) + (b[i] if i < len(b) else Fr(0)) for i in range(n)]


def poly_der(a):
    return [k * a[k] for k in range(1, len(a))] or [Fr(0)]


def open_knots(breaks, p, mult=None):
    """open knot vector on the sorted breakpoints; mult: {interior breakpoint: multiplicity} (default 1)"""
    breaks = [Fr(b) for b in breaks]
    mult = {Fr(k): v for k, v in (mult or {}).items()}
    t = [breaks[0]] * (p + 1)
    for b in breaks[1:-1]:
        t += [b] * mult.get(b, 1)
    t += [breaks[-1]] * (p + 1)
    return t


class Spline1D:
    def __init__(self, knots, p):
        self.t = [Fr(k) for k in knots]
        self.p = p
        t = self.t
        assert all(t[i] <= t[i + 1] for i in range(len(t) - 1))
        assert t[0] == t[p] and t[-1] == t[-p - 1], "knot vector must be open"
        self.n = len(t) - p - 1
        self.spans = [k for k in range(p, len(t) - p - 1) if t[k] < t[k + 1]]
        self.a = [t[k] for k in self.spans]
        self.b = [t[k + 1] for k in self.spans]
        self.h = [t[k + 1] - t[k] for k in self.spans]
        self.first = np.array([k - p for k in self.spans], dtype=np.int64)
        self.ne = len(self.spans)
        self.pieces = [self._pieces(k) for k in self.spans]      # [e][m] -> coefficient list in s (ascending)
        self.af = np.array([float(x) for x in self.a])
        self.hf = np.array([float(x) for x in self.h])

    def _pieces(self, k):
        t, p = self.t, self.p
        a, h = t[k], t[k + 1] - t[k]
        N = {k: [Fr(1)]}
        for d in range(1, p + 1):
            M = {}
            for j in range(k - d, k + 1):
                acc = [Fr(0)]
                if j in N:                                   # omega_{j,d} B_j^{d-1}
                    den = t[j + d] - t[j]
                    if den:
                        acc = poly_add(acc, poly_mul([(a - t[j]) / den, h / den], N[j]))
                if j + 1 in N:                               # (1 - omega_{j+1,d}) B_{j+1}^{d-1}
                    den = t[j + d + 1] - t[j + 1]
                    if den:
                        acc = poly_add(acc, poly_mul([(t[j + d + 1] - a) / den, -h / den], N[j + 1]))
                M[j] = acc
            N = M
        out = []
        for j in range(k - p, k + 1):
            c = N[j] + [Fr(0)] * (p + 1 - len(N[j]))
            out.append(c[:p + 1])
        return out

    def piece_x(self, e, m, d=0):
        """exact coefficients (in s) of the d-th PHYSICAL derivative of local function m on cell e"""
        c = self.pieces[e][m]
        for _ in range(d):
            c = poly_der(c)
        f = Fr(1) / self.h[e] ** d
        return [x * f for x in c]

    def eval_float(self, e, s, nder):
        """float values: e scalar or array (cell per point), s array -> [npt, p+1, nder+1] physical derivatives"""
        s = np.asarray(s, dtype=float)
        e = np.broadcast_to(np.asarray(e), s.shape)
        out = np.zeros(s.shape + (self.p + 1, nder + 1))
        for ee in np.unique(e):
            msk = e == ee
            ss = s[msk]
            for m in range(self.p + 1):
                c = self.pieces[ee][m]
                for d in range(nder + 1):
                    cf = [float(x) for x in c]
                    v = np.zeros_like(ss)
                    for x in reversed(cf):
                        v = v * ss + x
                    out[msk, m, d] = v / self.hf[ee] ** d
                    c = poly_der(c)
        return out

    def cell_of(self, x):
        """cell index for float points (right end belongs to the last cell)"""
        x = np.asarray(x, dtype=float)
        e = np.searchsorted(np.array([float(b) for b in self.b]), x, side='right')
        return np.minimum(e, self.ne - 1)

    def rig_table(self, e, nodes, nder):
        """rigorous enclosures at rational node intervals (list of RI in [0,1]) -> I [nq, p+1, nder+1]"""
        nq = len(nodes)
        lo = np.zeros((nq, self.p + 1, nder + 1))
        hi = np.zeros((nq, self.p + 1, nder + 1))
        for m in range(self.p + 1):
            for d in range(nder + 1):
                c = self.piece_x(e, m, d)
                for q, T in enumerate(nodes):
                    r = ri_horner(c, T)
                    lo[q, m, d] = frac_enclosure(r.lo)[0]
                    hi[q, m, d] = frac_enclosure(r.hi)[1]
        return I(lo, hi)


_GAUSS = {}


def gauss01(n):
    """float Gauss-Legendre rule on [0, 1]"""
    x, w = np.polynomial.legendre.leggauss(n)
    return (x + 1) / 2, w / 2


def gauss01_rig(n):
    """rational enclosures (RI) of the Gauss-Legendre nodes and weights on [0, 1] (curved3d.gauss_enclosures)"""
    if n not in _GAUSS:
        _GAUSS[n] = gauss_enclosures(n)
    return _GAUSS[n]
