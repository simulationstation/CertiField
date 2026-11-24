#!/usr/bin/env python3
"""CertiField magnetisation -- reference implementation (Python), thin films.

    E[m] = int_Omega  A |grad m|^2 + kappa m.(curl m) + K (1 - m3^2) + h (1 - m3) dx ,     m : Omega -> S^2 ,
    m = e3 on the boundary, degree Q.

docs/MAGNETISATION.md has the statements and proofs.  This file has

  the lower bound (Proposition 5):  a certificate is a pair of piecewise polynomials (phi, u) on [-1, 1] with rational
      coefficients and a rational number Gamma.  check_certificate verifies the polynomial inequalities of
      Proposition 5 in EXACT rational arithmetic (Bernstein coefficients with bisection); then
          E[m] >= 2 pi Q (int u dt - 4A) - Gamma |Omega|      for every admissible m of degree Q.
  the wall bound (Proposition 7) for the strip, checked in the same way.

Nothing here uses floating point for a certified statement except the final conversion of an exact rational bound to
a decimal, which is rounded in the safe direction.
"""
import math
from fractions import Fraction as Fr

__version__ = "0.1.0-reference"

# pi between two rationals (35 digits)
PI_LO = Fr(314159265358979323846264338327950288, 10 ** 35)
PI_HI = Fr(314159265358979323846264338327950289, 10 ** 35)


# ================================================================================================================
# exact polynomials in one variable: ascending lists of Fractions
# ================================================================================================================
def P_(*c):
    return [Fr(x) for x in c]


def padd(p, q, s=1):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) + s * (q[i] if i < len(q) else 0) for i in range(n)]


def pmul(p, q):
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, c in enumerate(q):
                r[i + j] += a * c
    return r


def pscale(p, c):
    return [Fr(c) * x for x in p]


def pder(p):
    return [i * c for i, c in enumerate(p)][1:] or [Fr(0)]


def pval(p, t):
    r = Fr(0)
    for c in reversed(p):
        r = r * t + c
    return r


def pint(p, a, b):
    return sum(c * (Fr(b) ** (i + 1) - Fr(a) ** (i + 1)) / (i + 1) for i, c in enumerate(p))


def pzero(p):
    return all(c == 0 for c in p)


def bernstein(p, a, b):
    """exact Bernstein coefficients of p on [a, b]"""
    a, b = Fr(a), Fr(b)
    n = len(p) - 1
    # q(s) = p(a + (b - a) s) in the monomial basis
    q = [Fr(0)]
    lin = [a, b - a]
    for c in reversed(p):
        q = padd(pmul(q, lin), [c])
    q = (q + [Fr(0)] * (n + 1))[:n + 1]
    return [sum(Fr(math.comb(k, i), math.comb(n, i)) * q[i] for i in range(k + 1)) for k in range(n + 1)]


def _split(bz):
    """de Casteljau at 1/2: the Bernstein coefficients on the two halves"""
    left, right = [bz[0]], [bz[-1]]
    cur = list(bz)
    while len(cur) > 1:
        cur = [(cur[i] + cur[i + 1]) / 2 for i in range(len(cur) - 1)]
        left.append(cur[0])
        right.append(cur[-1])
    return left, right[::-1]


def nonneg(p, a, b, depth=60):
    """True only if p >= 0 on [a, b]; exact.  (False means 'not proved', or a point with p < 0 was found.)"""
    stack = [(bernstein(p, a, b), 0)]
    while stack:
        bz, d = stack.pop()
        if min(bz) >= 0:
            continue
        if bz[0] < 0 or bz[-1] < 0 or d >= depth:
            return False
        l, r = _split(bz)
        stack.append((l, d + 1))
        stack.append((r, d + 1))
    return True


# ================================================================================================================
# the certificate of Proposition 5
# ================================================================================================================
class Model:
    """A > 0, kappa, K, h rational.  V(t) = K (1 - t^2) + h (1 - t)."""

    def __init__(self, A=Fr(1, 2), kappa=0, K=0, h=0):
        self.A, self.kappa, self.K, self.h = Fr(A), Fr(kappa), Fr(K), Fr(h)
        assert self.A > 0

    def V(self):
        return [self.K + self.h, -self.h, -self.K]

    def spec(self):
        return dict(A=str(self.A), kappa=str(self.kappa), K=str(self.K), h=str(self.h))

    @staticmethod
    def from_spec(s):
        return Model(Fr(s['A']), Fr(s['kappa']), Fr(s['K']), Fr(s['h']))


def hermite_cell(t0, t1, f0, d0, f1, d1):
    """the cubic on [t0, t1] with the given end values and slopes, as a polynomial in t (exact)"""
    t0, t1 = Fr(t0), Fr(t1)
    hh = t1 - t0
    s = [-t0 / hh, 1 / hh]                                   # the local variable (t - t0)/h
    basis = [P_(1, 0, -3, 2), P_(0, 1, -2, 1), P_(0, 0, 3, -2), P_(0, 0, -1, 1)]
    coef = [Fr(f0), Fr(d0) * hh, Fr(f1), Fr(d1) * hh]
    out = [Fr(0)]
    for bpoly, c in zip(basis, coef):
        if c == 0:
            continue
        q = [Fr(0)]
        for a in reversed(bpoly):
            q = padd(pmul(q, s), [a])
        out = padd(out, pscale(q, c))
    return out


class Certificate:
    """mesh: knots -1 = t_0 < ... < t_n = 1 (rationals);
    phi: per cell a polynomial in t (the cells must join in C^1: use from_hermite), u: per cell a polynomial in t;
    Gamma: rational."""

    def __init__(self, model, mesh, phi, u, Gamma=0):
        self.model = model
        self.mesh = [Fr(x) for x in mesh]
        self.phi = [[Fr(c) for c in p] for p in phi]
        self.u = [[Fr(c) for c in p] for p in u]
        self.Gamma = Fr(Gamma)
        assert len(self.mesh) == len(self.phi) + 1 == len(self.u) + 1

    @staticmethod
    def from_hermite(model, mesh, values, slopes, u, Gamma=0):
        phi = [hermite_cell(mesh[i], mesh[i + 1], values[i], slopes[i], values[i + 1], slopes[i + 1]) for i in range(len(mesh) - 1)]
        return Certificate(model, mesh, phi, u, Gamma)

    def polynomials(self, i):
        """P, R, V + Gamma and Pi_Gamma on cell i (exact)"""
        m = self.model
        phi, u = self.phi[i], self.u[i]
        a = pmul([Fr(1), Fr(1)], padd(phi, [m.kappa]))            # (1 + t)(phi + kappa)
        b = pmul([Fr(1), Fr(-1)], padd(phi, [m.kappa], -1))       # (1 - t)(phi - kappa)
        wp = pmul([Fr(1), Fr(-1)], pder(a))                        # (1 - t) a'
        wm = pmul([Fr(1), Fr(1)], pder(b))                         # (1 + t) b'
        Pp = pscale(pmul(wp, wp), Fr(1, 4))
        Rr = pscale(pmul(wm, wm), Fr(1, 4))
        VG = padd(m.V(), [self.Gamma])
        fourA_u = padd([4 * m.A], u, -1)
        Pi = padd(padd(pmul(Pp, fourA_u), pmul(Rr, u)), pmul(VG, pmul(u, fourA_u)), -1)
        return Pp, Rr, VG, Pi, fourA_u

    def smooth(self):
        """phi is C^1 across the knots"""
        for i in range(len(self.phi) - 1):
            t = self.mesh[i + 1]
            if pval(self.phi[i], t) != pval(self.phi[i + 1], t) or pval(pder(self.phi[i]), t) != pval(pder(self.phi[i + 1]), t):
                return False
        return True

    def check(self, verbose=False):
        """EXACT check of the conditions of Proposition 5.  Returns True only if they hold."""
        if not (self.mesh[0] == -1 and self.mesh[-1] == 1 and all(x < y for x, y in zip(self.mesh, self.mesh[1:]))):
            return False
        if not self.smooth():
            if verbose:
                print("phi is not C^1")
            return False
        A4 = 4 * self.model.A
        for i in range(len(self.u)):
            t0, t1 = self.mesh[i], self.mesh[i + 1]
            Pp, Rr, VG, Pi, fourA_u = self.polynomials(i)
            u = self.u[i]
            if pzero(u):
                ok = pzero(Pp) and nonneg(padd(pscale(VG, A4), Rr, -1), t0, t1)
            elif pzero(fourA_u):
                ok = pzero(Rr) and nonneg(padd(pscale(VG, A4), Pp, -1), t0, t1)
            else:
                ok = nonneg(u, t0, t1) and nonneg(fourA_u, t0, t1) and nonneg(pscale(Pi, -1), t0, t1)
            if not ok:
                if verbose:
                    print(f"condition fails on cell {i}: [{float(t0)}, {float(t1)}]")
                return False
        return True

    def integral_u(self):
        return sum(pint(self.u[i], self.mesh[i], self.mesh[i + 1]) for i in range(len(self.u)))

    def bound(self, Q, area=0):
        """the exact lower bound as (c, d):  E >= c pi + d,  c = 2 Q (int u - 4A),  d = -Gamma |Omega|.
        Valid only if check() is True."""
        c = 2 * Q * (self.integral_u() - 4 * self.model.A)
        d = -self.Gamma * Fr(area)
        return c, d

    def lower(self, Q, area=0):
        """the bound as an exact rational that is <= c pi + d"""
        c, d = self.bound(Q, area)
        return c * (PI_LO if c >= 0 else PI_HI) + d

    def to_json(self):
        return dict(format='CertiField-magnet-dual-1', model=self.model.spec(), Gamma=str(self.Gamma), mesh=[str(x) for x in self.mesh],
                    phi=[[str(c) for c in p] for p in self.phi], u=[[str(c) for c in p] for p in self.u])

    @staticmethod
    def from_json(d):
        return Certificate(Model.from_spec(d['model']), [Fr(x) for x in d['mesh']], [[Fr(c) for c in p] for p in d['phi']],
                           [[Fr(c) for c in p] for p in d['u']], Fr(d['Gamma']))


def constant_certificate(model, phi, u, Gamma=0):
    """phi a constant, u a polynomial, on the single cell [-1, 1]"""
    return Certificate(model, [-1, 1], [[Fr(phi)]], [[Fr(c) for c in u]], Gamma)


def floor_float(q):
    """the largest float that is <= the rational q"""
    x = float(q)
    while Fr(x) > q:
        x = math.nextafter(x, -math.inf)
    return x


def ceil_float(q):
    x = float(q)
    while Fr(x) < q:
        x = math.nextafter(x, math.inf)
    return x


# ================================================================================================================
# the wall bound of Proposition 7 (strip, kappa = 0)
# ================================================================================================================
def wall_bound(model, mesh, s, width):
    """s: per cell a polynomial in t (C^1 across knots is not needed for a Lipschitz s; continuity is checked).
    Returns the exact rational  w (s(1) - s(-1))  if  (1 - t^2) s'^2 <= 4 A V  holds on every cell, else None."""
    assert model.kappa == 0
    mesh = [Fr(x) for x in mesh]
    for i in range(len(s) - 1):
        if pval(s[i], mesh[i + 1]) != pval(s[i + 1], mesh[i + 1]) or pval(pder(s[i]), mesh[i + 1]) != pval(pder(s[i + 1]), mesh[i + 1]):
            return None
    for i, si in enumerate(s):
        d = pder([Fr(c) for c in si])
        g = padd(pscale(model.V(), 4 * model.A), pmul([Fr(1), Fr(0), Fr(-1)], pmul(d, d)), -1)
        if not nonneg(g, mesh[i], mesh[i + 1]):
            return None
    return Fr(width) * (pval(s[-1], Fr(1)) - pval(s[0], Fr(-1)))
