"""Exact and ball arithmetic shared by the rigorous evaluators.

- Polynomials in (x, y) with rational coefficients: python-flint fmpq_mpoly (exact).
- Transcendental and algebraic numbers (square roots, angles, sin(1/10)): python-flint arb balls at PREC bits;
  every arb operation returns a ball that contains the exact result (Arb's guarantee).
- Integrals of polynomials over rectangles: exact rationals.
- Integrals along circular arcs centred at the origin: exact reduction to the integrals of cos^n, whose values
  need only the endpoint angles and their sines and cosines (balls).
"""
from fractions import Fraction as Fr

from flint import fmpq_mpoly_ctx, fmpq, arb, ctx as _fctx

PREC = 320
_fctx.prec = PREC
CTX = fmpq_mpoly_ctx.get(('x', 'y'), 'lex')
X, Y = CTX.gens()
ZERO = CTX.from_dict({})


def q(v):
    """exact fmpq of a Fraction / int / float (a float is the binary rational it represents)"""
    if isinstance(v, fmpq):
        return v
    if isinstance(v, float):
        n, d = v.as_integer_ratio()
        return fmpq(n, d)
    v = Fr(v)
    return fmpq(v.numerator, v.denominator)


def fr_of_fmpq(v):
    return Fr(int(v.p), int(v.q))


def arb_bounds(b):
    """exact rational (lo, hi) with lo <= every point of the ball b <= hi"""
    out = []
    for e in (b.lower(), b.upper()):
        m, x = e.man_exp()
        m, x = int(m), int(x)
        out.append(Fr(m) * Fr(2) ** x)
    return out[0], out[1]


def uni(coefs, var):
    """polynomial sum coefs[k] var^k (var an mpoly generator)"""
    p = ZERO
    for c in reversed(coefs):
        p = p * var + q(c)
    return p


def local_to_x(coefs_s, a, h, var):
    """polynomial in var of sum c_k ((var - a)/h)^k"""
    t = (var - q(a)) * q(1 / Fr(h))
    p = ZERO
    for c in reversed(coefs_s):
        p = p * t + q(c)
    return p


def rect_integral(p, x0, x1, y0, y1):
    """exact integral of the mpoly p over [x0,x1] x [y0,y1] (fmpq)"""
    px = p.integral('x')
    px = px.subs({'x': q(x1)}) - px.subs({'x': q(x0)})
    py = px.integral('y')
    v = py.subs({'y': q(y1)}) - py.subs({'y': q(y0)})
    return v.coefficient(0) if not v.is_zero() else fmpq(0)


def eval_arb(p, xv, yv):
    """ball value of the mpoly p at ball (or exact) arguments"""
    xv, yv = arb(xv), arb(yv)
    tot = arb(0)
    cache_x, cache_y = {}, {}
    for (i, j), c in p.terms():
        if i not in cache_x:
            cache_x[i] = xv ** i
        if j not in cache_y:
            cache_y[j] = yv ** j
        tot += arb(c) * cache_x[i] * cache_y[j]
    return tot


class APoint:
    """a point on the circle of radius R centred at the origin: angle t, cos t, sin t as balls"""
    __slots__ = ('t', 'c', 's')

    def __init__(self, t, c, s):
        self.t, self.c, self.s = t, c, s

    @staticmethod
    def from_x(R, x):
        """the point with abscissa x (rational, 0 <= x <= R) in the first quadrant"""
        c = arb(q(Fr(x) / Fr(R)))
        s = (1 - c * c).sqrt()
        return APoint(c.acos(), c, s)

    @staticmethod
    def from_y(R, y):
        s = arb(q(Fr(y) / Fr(R)))
        c = (1 - s * s).sqrt()
        return APoint(s.asin(), c, s)


def cos_power_integrals(nmax, lo, hi):
    """I_n = int_{lo.t}^{hi.t} cos^n t dt, n = 0..nmax (balls)"""
    I = [hi.t - lo.t, hi.s - lo.s]
    for n in range(2, nmax + 1):
        I.append((hi.c ** (n - 1) * hi.s - lo.c ** (n - 1) * lo.s) / n + I[n - 2] * fmpq(n - 1, n))
    return I[:nmax + 1]


def arc_integral(F, R, wc, ws, lo, hi):
    """int_{lo.t}^{hi.t} F(R cos t, R sin t) (wc cos t + ws sin t) dt  for an mpoly F with rational coefficients and
    rational R, wc, ws: reduced exactly to A(cos t) + sin t B(cos t)."""
    R, wc, ws = q(R), q(wc), q(ws)
    A, B = {}, {}

    def add(i, j, c):
        # c * cos^i * sin^j  with sin^2 = 1 - cos^2
        target = A if j % 2 == 0 else B
        m = j // 2
        binom = 1
        for k in range(m + 1):              # (1 - c^2)^m = sum_k (-1)^k binom(m,k) c^{2k}
            target[i + 2 * k] = target.get(i + 2 * k, fmpq(0)) + c * ((-1) ** k * binom)
            binom = binom * (m - k) // (k + 1)

    for (i, j), c in F.terms():
        cc = c * R ** (i + j)
        if wc != 0:
            add(i + 1, j, cc * wc)
        if ws != 0:
            add(i, j + 1, cc * ws)
    tot = arb(0)
    if A:
        I = cos_power_integrals(max(A), lo, hi)
        for n, a in A.items():
            if a != 0:
                tot += arb(a) * I[n]
    # int sin t B(cos t) dt = [-Btilde(cos t)]
    for n, b in B.items():
        if b != 0:
            tot += arb(b) * (lo.c ** (n + 1) - hi.c ** (n + 1)) / (n + 1)
    return tot


def ray_integral(F, c0, s0, w, r0, r1):
    """int_{r0}^{r1} F(r c0, r s0) w dr for an mpoly F, balls c0, s0, w, r0, r1"""
    tot = arb(0)
    for (i, j), c in F.terms():
        n = i + j
        tot += arb(c) * c0 ** i * s0 ** j * (r1 ** (n + 1) - r0 ** (n + 1)) / (n + 1)
    return tot * w


def disc_part_integral(p, x0, x1, y0, y1, a=Fr(1, 2)):
    """ball enclosure of the integral of the mpoly p over the cell [x0,x1]x[y0,y1] (x0, y0 >= 0) intersected with
    the open disc x^2 + y^2 < a^2.  The cell must be cut by the circle."""
    x0, x1, y0, y1, a = Fr(x0), Fr(x1), Fr(y0), Fr(y1), Fr(a)
    a2 = a * a
    assert 0 <= x0 < x1 and 0 <= y0 < y1 and x0 * x0 + y0 * y0 < a2 < x1 * x1 + y1 * y1
    Qy = p.integral('y')
    Qy = Qy - Qy.subs({'y': q(y0)})
    tot = arb(0)
    has_xs = y1 < a and x0 * x0 + y1 * y1 < a2           # the circle meets y = y1 at xs > x0
    if has_xs:                                            # x in [x0, xs]: full height
        H = Qy.subs({'y': q(y1)}).integral('x')
        xs = arb(q(a2 - y1 * y1)).sqrt()
        tot += eval_arb(H, xs, 0) - eval_arb(H, q(x0), 0)
        pa = APoint.from_y(a, y1)
    else:
        pa = APoint.from_x(a, x0)
    if x1 * x1 + y0 * y0 < a2:                            # right end of the arc part: x = x1
        pb = APoint.from_x(a, x1)
    else:                                                 # the circle meets y = y0
        pb = APoint.from_y(a, y0)
    # x from xa to xb along y = sqrt(a^2 - x^2):  x = a cos t, dx = -a sin t dt, t from pa.t down to pb.t
    tot += arc_integral(Qy, a, 0, a, pb, pa)
    return tot
