"""P3: integrals of a trivariate polynomial over (cell minus ball), exactly up to Arb balls.

For a cell C = [x0,x1] x [y0,y1] x [z0,z1] in the first octant that is cut by the sphere r = a,

    int_{C minus ball} p = int_C p (exact rational) - int_{C cap ball} p,

and with Q(x,y,z) = int_{z0}^{z} p dz (exact), Z(x,y) = sqrt(a^2 - x^2 - y^2), R = [x0,x1] x [y0,y1]:

    int_{C cap ball} p = int_{R cap D(r1)} Q(x,y,z1) + int_{R cap (D(r0) minus D(r1))} Q(x,y,Z),

D(r) the disc x^2 + y^2 <= r^2 (r0^2 = a^2 - z0^2, r1^2 = a^2 - z1^2, D(r1) empty if z1 >= a).  Q(x,y,Z) = P + Z W
with polynomials P, W (Z^2 = a^2 - x^2 - y^2).  The integrals of P over R cap D(r) are reduced exactly to integrals
along circle arcs (exact.arc_integral with an Arb radius); the integrals of Z W are integrated exactly in y (closed
forms with sqrt and asin) and then in x by Arb's rigorous integration (acb.integral, Petras' algorithm), which encloses
the integral of the analytic integrand; branch cuts are checked (analytic flag).
"""
from fractions import Fraction as Fr

from flint import fmpq_mpoly_ctx, fmpq, arb, acb

from exact import q, APoint, cos_power_integrals

CTX3 = fmpq_mpoly_ctx.get(('x', 'y', 'z'), 'lex')
X3, Y3, Z3 = CTX3.gens()
ZERO3 = CTX3.from_dict({})


def uni3(coefs, var):
    p = ZERO3
    for c in reversed(coefs):
        p = p * var + q(c)
    return p


def local3(coefs_s, a, h, var):
    t = (var - q(a)) * q(1 / Fr(h))
    p = ZERO3
    for c in reversed(coefs_s):
        p = p * t + q(c)
    return p


def box_integral(p, box):
    x0, x1, y0, y1, z0, z1 = [q(v) for v in box]
    for v, lo, hi in (('x', x0, x1), ('y', y0, y1), ('z', z0, z1)):
        p = p.integral(v)
        p = p.subs({v: hi}) - p.subs({v: lo})
    return p.coefficient(0) if not p.is_zero() else fmpq(0)


# ---------------------------------------------------------------------------------------------------------------------
# 2D pieces over R cap D(r), r^2 = r2 rational, R in the first quadrant
def _by_y(poly2):
    """dict j -> list of (i, coefficient) for a polynomial in x, y (z absent)"""
    out = {}
    for (i, j, k), c in poly2.terms():
        assert k == 0
        out.setdefault(j, []).append((i, c))
    return out


def _eval_xy(poly, xv, yv):
    tot = arb(0)
    for (i, j, k), c in poly.terms():
        tot += arb(c) * xv ** i * yv ** j
    return tot


def _arc_integral_arb(F, R, wc, ws, lo, hi):
    """int_{lo.t}^{hi.t} F(R cos t, R sin t)(wc cos t + ws sin t) dt for a polynomial F(x, y) and an Arb radius R"""
    A, B = {}, {}

    def add(i, j, c):
        target = A if j % 2 == 0 else B
        m = j // 2
        binom = 1
        for k in range(m + 1):
            target[i + 2 * k] = target.get(i + 2 * k, arb(0)) + c * ((-1) ** k * binom)
            binom = binom * (m - k) // (k + 1)

    for (i, j, k), c in F.terms():
        cc = arb(c) * R ** (i + j)
        if wc != 0:
            add(i + 1, j, cc * wc)
        if ws != 0:
            add(i, j + 1, cc * ws)
    tot = arb(0)
    if A:
        I = cos_power_integrals(max(A), lo, hi)
        for n, a in A.items():
            tot += a * I[n]
    for n, b in B.items():
        tot += b * (lo.c ** (n + 1) - hi.c ** (n + 1)) / (n + 1)
    return tot


def _apt_x(R, x):
    c = arb(q(x)) / R
    s = (1 - c * c).sqrt()
    return APoint(c.acos(), c, s)


def _apt_y(R, y):
    s = arb(q(y)) / R
    c = (1 - s * s).sqrt()
    return APoint(s.asin(), c, s)


def _layout(rect, r2):
    """pieces of R cap D(r): ('rect', xa, xb) full height in x in [xa, xb], ('arc', pa, pb) top on the circle.
    Returns None if R cap D(r) is empty, 'all' if R lies in the closed disc"""
    x0, x1, y0, y1 = rect
    if x0 * x0 + y0 * y0 >= r2:
        return None
    if x1 * x1 + y1 * y1 <= r2:
        return 'all'
    pieces = []
    has_xs = y1 * y1 < r2 and x0 * x0 + y1 * y1 < r2            # the circle meets y = y1 at xs > x0
    if has_xs:
        pieces.append(('rect', Fr(x0), ('sqrt', r2 - y1 * y1)))  # x in [x0, xs]: y up to y1
    right = ('x', x1) if x1 * x1 + y0 * y0 < r2 else ('y', y0)  # where the arc part ends
    left = ('y', y1) if has_xs else ('x', x0)
    pieces.append(('arc', left, right))
    return pieces


def poly_disc(P, rect, r2):
    """ball of int_{R cap D(r)} P(x, y) dx dy, P a polynomial in x, y; r2 rational"""
    lay = _layout(rect, r2)
    x0, x1, y0, y1 = [q(v) for v in rect]
    if lay is None:
        return arb(0)
    if lay == 'all':
        p = P.integral('x')
        p = p.subs({'x': x1}) - p.subs({'x': x0})
        p = p.integral('y')
        p = p.subs({'y': y1}) - p.subs({'y': y0})
        return arb(p.coefficient(0)) if not p.is_zero() else arb(0)
    R = arb(q(r2)).sqrt()
    Py = P.integral('y')
    Py = Py - Py.subs({'y': y0})
    tot = arb(0)
    for piece in lay:
        if piece[0] == 'rect':
            H = Py.subs({'y': y1}).integral('x')
            xs = arb(q(piece[2][1])).sqrt()
            tot += _eval_xy(H, xs, arb(0)) - _eval_xy(H, arb(x0), arb(0))
        else:
            left, right = piece[1], piece[2]
            pa = _apt_y(R, left[1]) if left[0] == 'y' else _apt_x(R, left[1])
            pb = _apt_x(R, right[1]) if right[0] == 'x' else _apt_y(R, right[1])
            # x from xa to xb along y = sqrt(r^2 - x^2): x = R cos t, dx = -R sin t dt, t from pa.t down to pb.t
            tot += _arc_integral_arb(Py, R, 0, R, pb, pa)
    return tot


# ---------------------------------------------------------------------------------------------------------------------
def _sqrt_a(z, analytic):
    return z.sqrt(analytic=analytic)


def _asin_a(u, analytic):
    if analytic and (u.real.upper() >= 1 or u.real.lower() <= -1):
        return acb('nan')
    return u.asin()


def _K(jmax, y, b2, analytic, s=None, asn=None):
    """K_j(y) = int_0^y t^j sqrt(b^2 - t^2) dt, j = 0..jmax (complex balls); s = sqrt(b^2 - y^2) and
    asn = asin(y / b) may be given when they are known in closed form"""
    b = _sqrt_a(b2, analytic)
    if s is None:
        s = _sqrt_a(b2 - y * y, analytic)
    if asn is None:
        asn = _asin_a(y / b, analytic)
    K = [(y * s + b2 * asn) / 2, (b2 * b - s * s * s) / 3]
    for j in range(2, jmax + 1):
        K.append(-y ** (j - 1) * s * s * s / (j + 2) + K[j - 2] * b2 * fmpq(j - 1, j + 2))
    return K[:jmax + 1]


def hemi_disc(W, rect, r2, a2, rel_bits=60):
    """ball of int_{R cap D(r)} sqrt(a^2 - x^2 - y^2) W(x, y) dx dy, with r^2 = r2 <= a2 (rationals)"""
    lay = _layout(rect, r2)
    if lay is None or W.is_zero():
        return arb(0)
    x0, x1, y0, y1 = [Fr(v) for v in rect]
    byj = _by_y(W)
    jmax = max(byj)
    a2c = acb(q(a2))
    y0c, y1c = acb(q(y0)), acb(q(y1))
    r2c = acb(q(r2))

    s_circ = acb(q(a2 - r2)).sqrt()                 # on the circle x^2 + y^2 = r^2: sqrt(b^2 - y^2) = sqrt(a^2 - r^2)
    equator = (r2 == a2)                              # then y / b = 1 on the whole circle: asin = pi/2

    def make_f(top):
        def f(x, analytic):
            b2 = a2c - x * x
            if top == 'y1':
                Ku = _K(jmax, y1c, b2, analytic)
            else:
                yu = _sqrt_a(r2c - x * x, analytic)
                Ku = _K(jmax, yu, b2, analytic, s=s_circ, asn=(acb(arb.pi() / 2) if equator else None))
            Kl = _K(jmax, y0c, b2, analytic) if y0 != 0 else [acb(0)] * (jmax + 1)
            tot = acb(0)
            for j, lst in byj.items():
                poly = acb(0)
                for (i, c) in lst:
                    poly += acb(c) * x ** i
                tot += poly * (Ku[j] - Kl[j])
            return tot
        return f

    if lay == 'all':
        segs = [(Fr(x0), Fr(x1), 'y1')]
    else:
        segs = []
        for piece in lay:
            if piece[0] == 'rect':
                segs.append((Fr(x0), piece[2], 'y1'))
            else:
                left, right = piece[1], piece[2]
                xa = ('sqrt', r2 - left[1] ** 2) if left[0] == 'y' else Fr(left[1])
                xb = Fr(right[1]) if right[0] == 'x' else ('sqrt', r2 - right[1] ** 2)
                segs.append((xa, xb, 'circle'))
    tot = arb(0)
    for xa, xb, top in segs:
        tot += _integrate(make_f(top), xa, xb, rel_bits)
    return tot


def _endpoint(v):
    """an exact float-like arb midpoint m and a radius bound e with |v - m| <= e"""
    if isinstance(v, tuple):
        b = arb(q(v[1])).sqrt()
        return arb(b.mid()), b.rad()
    return arb(q(v)), arb(0)


def _integrate(f, xa, xb, rel_bits):
    """rigorous enclosure of int_{xa}^{xb} f(x) dx (endpoints rational or ('sqrt', r) for sqrt(r))"""
    ma, ea = _endpoint(xa)
    mb, eb = _endpoint(xb)
    val = acb.integral(f, acb(ma), acb(mb), rel_tol=arb(2) ** (-rel_bits))
    out = val.real
    # endpoint uncertainty: |int over [m, v]| <= |v - m| * sup |f| on a ball around the endpoint
    for m, e in ((ma, ea), (mb, eb)):
        if e != 0:
            ball = acb(arb(m, 2 * e))
            fb = f(ball, False)
            out += arb(0, 2 * e * fb.abs_upper())
    return out


# ---------------------------------------------------------------------------------------------------------------------
def ball_part(p, box, a=Fr(1, 2)):
    """ball of int_{box cap {r < a}} p (box in the first octant, cut by the sphere)"""
    x0, x1, y0, y1, z0, z1 = [Fr(v) for v in box]
    a2 = Fr(a) ** 2
    rect = (x0, x1, y0, y1)
    Qz = p.integral('z')
    Qz = Qz - Qz.subs({'z': q(z0)})
    tot = arb(0)
    r0sq = a2 - z0 * z0
    r1sq = a2 - z1 * z1
    if r1sq > 0:
        tot += poly_disc(Qz.subs({'z': q(z1)}), rect, r1sq)
    # Q(x, y, Z) = P + Z W with Z^2 = a^2 - x^2 - y^2
    S2 = q(a2) - X3 * X3 - Y3 * Y3
    P, W = ZERO3, ZERO3
    pw = {0: CTX3.from_dict({(0, 0, 0): fmpq(1)})}
    for (i, j, k), c in Qz.terms():
        m = k // 2
        if m not in pw:
            pw[m] = S2 ** m
        term = CTX3.from_dict({(i, j, 0): c}) * pw[m]
        if k % 2 == 0:
            P += term
        else:
            W += term
    tot += poly_disc(P, rect, r0sq) + hemi_disc(W, rect, r0sq, a2)
    if r1sq > 0:
        tot -= poly_disc(P, rect, r1sq) + hemi_disc(W, rect, r1sq, a2)
    return tot


def outside_integral(p, box, a=Fr(1, 2)):
    """ball of int_{box minus ball} p"""
    return arb(box_integral(p, box)) - ball_part(p, box, a)
