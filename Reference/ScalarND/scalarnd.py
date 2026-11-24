#!/usr/bin/env python3
"""CertiField ND -- reference implementation (Python/NumPy).

Certified two-sided bounds for scalar variational problems on boxes in n = 1, 2, 3 dimensions,

    E[u] = int_Omega  k |grad u|^2 + W(u) dx ,      u = g on the Dirichlet faces, free on the natural faces,

with a user-supplied polynomial potential W (it may be non-convex).  NOTE.md has the theorem and its proof.

  upper bound : the energy of the piecewise multilinear (Q1) function with the given nodal values, enclosed rigorously.
  lower bound : for ANY field   S(x, t) = t sigma0(x) + sum_m a_m(x) beta_m(clamp(t, tlo, thi))
                (sigma0, a_m lowest-order Raviart-Thomas on the box grid, beta_m uniform C^1 quadratic B-splines),
                every admissible u satisfies

                    E[u]  >=  int_{Gamma_D} S(x, g) . nu   -   int_Omega sup_t Phi(x, t) dx ,
                    Phi   =   div_x S + |d_t S|^2 / (4k) - W(t) .

                Two ways of choosing S are provided (both heuristic: they affect tightness, never validity):
                  convex mode      (K = 0)  S = t sigma0, sigma0 from a cell-centred finite-volume solve.
                                            This is the classical complementary-energy bound; tight for convex W.
                  calibration mode (K > 0)  the spline part alone, optimised by L-BFGS on a smoothed copy of the
                                            certified functional, coarse-to-fine.  Needed when W is not convex.

Everything in upper_bound, flux_lower, lower_bound, conj_upper, conj_cells is rigorous: outward-rounded interval
arithmetic (ia.py) on exact inputs (binary floats and rationals).  Everything else is plain floating point.
"""
import itertools
import json
import math
import time
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import minimize

from ia import I, fr, horner, poly_deriv, dn, up, frac_enclosure, xlnx, sqrt_up, pi_lower

INF = np.inf
__version__ = "0.3.0-reference"


# ================================================================================================================
# potentials (docs/SCALAR_ND_2.md, Section 1).  A list of ascending rational coefficients, as in version 0.3.0, is a
# Poly on the whole line.
# ================================================================================================================
class Poly:
    """W(t) = sum_j W[j] t^j with rational coefficients on the value range `domain` = (lo, hi), +infinity outside;
    domain None is the whole line"""
    kind = 'poly'

    def __init__(self, W, domain=None):
        self.W = [Fr(c) for c in W]
        self.W1 = poly_deriv(self.W)
        self.W2 = poly_deriv(self.W1)
        self.Wf = np.array([float(c) for c in self.W])
        self.W1f = np.array([float(c) for c in self.W1])
        self.W2f = np.array([float(c) for c in self.W2])
        self.domain = None if domain is None else (Fr(domain[0]), Fr(domain[1]))
        self.c1 = self.c2 = Fr(0)
        assert self.domain is None or self.domain[1] > self.domain[0]

    # ---- floating point, for the heuristics only
    def f(self, t):
        return np.polynomial.polynomial.polyval(t, self.Wf)

    def df(self, t):
        return np.polynomial.polynomial.polyval(t, self.W1f)

    def d2f(self, t):
        return np.polynomial.polynomial.polyval(t, self.W2f)

    # ---- rigorous
    def enc(self, t):
        """enclosure of W at the exact rational t"""
        return fr(sum(c * t ** q for q, c in enumerate(self.W)))

    def enc_iv(self, x):
        """enclosure of W over the interval array x (inside the value range)"""
        return horner(self.W, x)

    def w2_iv(self, x):
        """enclosure of W'' over the interval array x; the upper end is +inf where x touches a singular end"""
        return horner(self.W2, x)

    def sub_bound(self, a, b):
        """For the subinterval [a, b] (exact rationals) of the value range: (h2, ca, cb) such that, for C convex and
        D a quadratic polynomial,
            sup_[a,b] (C + D - W) <= max(Phi(a) + ca, Phi(b) + cb) + (b - a)^2/8 (h2 - D'')^+
        (Lemma 3 of docs/SCALAR_ND.md; Lemma 11 of docs/SCALAR_ND_2.md at a singular end).  h2 is an upper bound for
        the second derivative of the part of W that is regular on [a, b]."""
        return float(horner(self.W2, I(frac_enclosure(a)[0], frac_enclosure(b)[1])).hi), 0.0, 0.0

    def spec(self):
        return dict(kind=self.kind, W=[str(c) for c in self.W],
                    domain=None if self.domain is None else [str(x) for x in self.domain])


class PolyLog(Poly):
    """W(t) = P(t) + c1 (t - lo) ln(t - lo) + c2 (hi - t) ln(hi - t) on [lo, hi], +infinity outside; c1, c2 >= 0.
    The entropy terms are 0 at the ends.  (t ln t + (1 - t) ln(1 - t) + chi t (1 - t): PolyLog([0, chi, -chi], 1, 1,
    (0, 1)).)"""
    kind = 'polylog'

    def __init__(self, W, c1, c2, domain):
        super().__init__(W, domain)
        assert self.domain is not None
        self.c1, self.c2 = Fr(c1), Fr(c2)
        assert self.c1 >= 0 and self.c2 >= 0
        self._lo, self._hi, self._c1, self._c2 = (float(x) for x in (self.domain[0], self.domain[1], self.c1, self.c2))

    @staticmethod
    def _ent(s):
        s = np.asarray(s, dtype=float)
        return np.where(s > 0, s * np.log(np.where(s > 0, s, 1.0)), 0.0)

    def f(self, t):
        return super().f(t) + self._c1 * self._ent(t - self._lo) + self._c2 * self._ent(self._hi - t)

    def df(self, t):
        with np.errstate(divide='ignore', invalid='ignore'):
            return super().df(t) + self._c1 * (np.log(t - self._lo) + 1.0) - self._c2 * (np.log(self._hi - t) + 1.0)

    def d2f(self, t):
        with np.errstate(divide='ignore', invalid='ignore'):
            return super().d2f(t) + self._c1 / (t - self._lo) + self._c2 / (self._hi - t)

    def enc(self, t):
        lo, hi = self.domain
        assert lo <= t <= hi
        v = super().enc(t)
        if self.c1:
            v = v + xlnx(fr(t - lo)) * fr(self.c1)
        if self.c2:
            v = v + xlnx(fr(hi - t)) * fr(self.c2)
        return v

    def _dist(self, x):
        """enclosures of x - lo and hi - x for x inside the value range (so both are >= 0)"""
        lo, hi = self.domain
        s1, s2 = x - fr(lo), fr(hi) - x
        return (I(np.maximum(s1.lo, 0.0), np.maximum(s1.hi, 0.0)), I(np.maximum(s2.lo, 0.0), np.maximum(s2.hi, 0.0)))

    def enc_iv(self, x):
        v = horner(self.W, x)
        s1, s2 = self._dist(x)
        if self.c1:
            v = v + xlnx(s1) * fr(self.c1)
        if self.c2:
            v = v + xlnx(s2) * fr(self.c2)
        return v

    def w2_iv(self, x):
        v = horner(self.W2, x)
        for c, s in zip((self.c1, self.c2), self._dist(x)):
            if c:
                pos = s.lo > 0                                       # c / s on [s.lo, s.hi]; unbounded above if s.lo = 0
                t = I(np.where(pos, s.lo, 1.0), np.where(pos, s.hi, 1.0)).inv() * fr(c)
                v = v + I(np.where(pos, t.lo, 0.0), np.where(pos, t.hi, INF))
        return v

    def sub_bound(self, a, b):
        lo, hi = self.domain
        assert lo <= a < b <= hi
        w = b - a
        ca = float(frac_enclosure(self.c1 * w)[1]) if a == lo else 0.0          # Lemma 11: the entropy term costs c w
        cb = float(frac_enclosure(self.c2 * w)[1]) if b == hi else 0.0
        x = I(frac_enclosure(a)[0], frac_enclosure(b)[1])
        h2 = horner(self.W2, x)
        if self.c1 and a != lo:
            h2 = h2 + (x - fr(lo)).inv() * fr(self.c1)
        if self.c2 and b != hi:
            h2 = h2 + (fr(hi) - x).inv() * fr(self.c2)
        return float(h2.hi), ca, cb

    def spec(self):
        return dict(super().spec(), c=[str(self.c1), str(self.c2)])


def potential_from_spec(s):
    dom = None if s.get('domain') is None else (Fr(s['domain'][0]), Fr(s['domain'][1]))
    if s['kind'] == 'poly':
        return Poly([Fr(c) for c in s['W']], dom)
    if s['kind'] == 'polylog':
        return PolyLog([Fr(c) for c in s['W']], Fr(s['c'][0]), Fr(s['c'][1]), dom)
    raise ValueError("unknown kind of potential: %r" % (s['kind'],))


# ================================================================================================================
# the problem
# ================================================================================================================
class Problem:
    """box: [(a_1, b_1), ...] rationals;  N: cells per direction;  k: stiffness (rational > 0);
    W: ascending polynomial coefficients (rationals);  bc[i] = (type at x_i = a_i, type at x_i = b_i), 'D' or 'N';
    g: function of the node coordinate arrays giving the Dirichlet nodal values (each float is taken as the exact
    binary number it is).  THE DIRICHLET DATUM OF THE CERTIFIED PROBLEM IS THE PIECEWISE MULTILINEAR INTERPOLANT OF
    THESE NODAL VALUES on the boundary grid (equal to g itself when g is multilinear on each boundary face).

    Extensions (docs/SCALAR_ND_2.md).  W may also be a Poly or PolyLog object, or a list of them together with
    mat, an integer array over the cells that says which one acts in each cell (they share one value range).
    lin: the coefficient a of a term a(x) u in the energy density: a function of the node coordinate arrays or a node
    array (a is then the continuous piecewise multilinear function with these nodal values) or an array over the cells
    (a constant on each cell); every float is taken as the exact binary number it is.
    mean: the prescribed mean of u (a rational), or None.
    bc[i] = ('P', 'P'): direction i is periodic."""

    def __init__(self, box, N, k, W, bc, g=None, name='', gval=None, lin=None, mean=None, mat=None):
        self.name = name
        self.box = [(Fr(a), Fr(b)) for a, b in box]
        self.N = tuple(int(x) for x in N)
        self.n = len(self.N)
        self.k = Fr(k)
        if isinstance(W, Poly):
            self.pots = [W]
        elif len(W) > 0 and isinstance(W[0], Poly):
            self.pots = list(W)
        else:
            self.pots = [Poly(W)]
        self.W = self.pots[0].W                    # as in version 0.3.0; it is the potential when there is only one
        self.domain = self.pots[0].domain
        assert all(p.domain == self.domain for p in self.pots), "the potentials must share one value range"
        self.mat = None if mat is None else np.array(mat, dtype=int)
        assert (self.mat is None and len(self.pots) == 1) or (self.mat is not None and self.mat.shape == self.N
                                                                and self.mat.min() >= 0 and self.mat.max() < len(self.pots))
        self.plain = len(self.pots) == 1 and type(self.pots[0]) is Poly and self.domain is None
        self.bc = [tuple(x) for x in bc]
        assert all(set(b) <= {'D', 'N'} or b == ('P', 'P') for b in self.bc), "bc entries: 'D', 'N', or ('P', 'P')"
        self.periodic = [b == ('P', 'P') for b in self.bc]
        self.mean = None if mean is None else Fr(mean)
        assert len(self.box) == self.n == len(self.bc) and self.k > 0 and all(b > a for a, b in self.box)
        self.h = [(b - a) / Ni for (a, b), Ni in zip(self.box, self.N)]
        self.hf = [float(x) for x in self.h]
        self.cellvol = math.prod(self.h)                                   # exact
        self.vol = self.cellvol * math.prod(self.N)                        # exact
        self.axes = [np.array([float(a + j * h) for j in range(Ni + 1)])
                     for (a, b), Ni, h in zip(self.box, self.N, self.h)]
        fixed = np.zeros([Ni + 1 for Ni in self.N], dtype=bool)
        for i in range(self.n):
            for side in (0, 1):
                if self.bc[i][side] == 'D':
                    sl = [slice(None)] * self.n
                    sl[i] = 0 if side == 0 else -1
                    fixed[tuple(sl)] = True
        self.fixed = fixed
        if gval is not None:
            self.gval = np.array(gval, dtype=float)
            assert self.gval.shape == fixed.shape
        else:
            self.gval = np.zeros(fixed.shape)
            if g is not None:
                X = np.meshgrid(*self.axes, indexing='ij')
                gv = np.asarray(g(*X), dtype=float) * np.ones(fixed.shape)
                self.gval[fixed] = gv[fixed]
        self.gval[~fixed] = 0.0
        self.gval = self.tie(self.gval)
        assert np.all(np.isfinite(self.gval))
        self.verts = list(itertools.product((0, 1), repeat=self.n))
        self.Wf = np.array([float(c) for c in self.W])
        self.W1 = poly_deriv(self.W)
        self.W2 = poly_deriv(self.W1)
        self.W1f = np.array([float(c) for c in self.W1])
        self.W2f = np.array([float(c) for c in self.W2])
        # the term a(x) u: av[v] = the value of a at the vertex v of every cell (an array over the cells)
        self.lin_nodes = self.lin_cells = self.av = None
        if lin is not None:
            if callable(lin):
                X = np.meshgrid(*self.axes, indexing='ij')
                L = np.asarray(lin(*X), dtype=float) * np.ones(fixed.shape)
            else:
                L = np.array(lin, dtype=float)
            assert np.all(np.isfinite(L))
            if L.shape == fixed.shape:
                self.lin_nodes = self.tie(L)
                self.av = {v: self.lin_nodes[self.vslice(v)] for v in self.verts}
            elif L.shape == self.N:
                self.lin_cells = L
                self.av = {v: L for v in self.verts}
            else:
                raise ValueError("lin must be an array over the nodes or over the cells")
        Wtrim = list(self.W)
        while len(Wtrim) > 1 and Wtrim[-1] == 0:
            Wtrim.pop()
        self.linear = self.plain and len(Wtrim) <= 2          # W = w0 + w1 t: Proposition 16 of docs/SCALAR_ND_2.md
        self.affine_ok = self.linear or all(p.domain is None and is_superlinear(p.W) for p in self.pots)

    def tie(self, A):
        """in every periodic direction the last node (or face) is the first one: copy it there"""
        for i in range(self.n):
            if self.periodic[i]:
                first = [slice(None)] * self.n
                first[i] = 0
                last = [slice(None)] * self.n
                last[i] = -1
                A[tuple(last)] = A[tuple(first)]
        return A

    def tied(self, A):
        return all(not self.periodic[i] or np.array_equal(np.take(A, 0, axis=i), np.take(A, -1, axis=i))
                   for i in range(self.n))

    # ---- the potential of each cell.  One polynomial: scalars, exactly as in version 0.3.0.  Several: arrays over cells.
    def _sel(self, vals):
        if self.mat is None:
            return vals[0]
        if isinstance(vals[0], I):
            return I(self._sel([v.lo for v in vals]), self._sel([v.hi for v in vals]))
        out = np.broadcast_to(np.asarray(vals[0], dtype=float), self.N).copy()
        for j in range(1, len(vals)):
            out = np.where(self.mat == j, vals[j], out)
        return out

    def Wenc(self, t):
        """enclosure of W_c at the exact rational t"""
        return self._sel([p.enc(t) for p in self.pots])

    def sub_bound(self, a, b):
        r = [p.sub_bound(a, b) for p in self.pots]
        return tuple(self._sel([x[q] for x in r]) for q in range(3))

    def Wenc_cells(self, x):
        """enclosure of W_c over the interval array x (one interval per cell)"""
        return self._sel([p.enc_iv(x) for p in self.pots])

    def W2rng_cells(self, x):
        return self._sel([p.w2_iv(x) for p in self.pots])

    def vslice(self, v):
        """for every cell, its vertex v in {0,1}^n, from a node array"""
        return tuple(slice(vi, Ni + vi) for vi, Ni in zip(v, self.N))

    def fslice(self, i, side):
        """for every cell, its lower (side 0) / upper (side 1) face in direction i, from a face array of direction i"""
        sl = [slice(None)] * self.n
        sl[i] = slice(0, self.N[i]) if side == 0 else slice(1, self.N[i] + 1)
        return tuple(sl)

    def face_shape(self, i):
        return tuple(Ni + (1 if j == i else 0) for j, Ni in enumerate(self.N))

    def Wfloat(self, t):
        if self.plain:
            return np.polynomial.polynomial.polyval(t, self.Wf)
        return self._sel([p.f(t) for p in self.pots])

    def dWfloat(self, t):
        if self.plain:
            return np.polynomial.polynomial.polyval(t, self.W1f)
        return self._sel([p.df(t) for p in self.pots])

    def spec(self):
        s = dict(name=self.name, box=[[str(a), str(b)] for a, b in self.box], N=list(self.N), k=str(self.k),
                 W=[str(c) for c in self.W], bc=[list(x) for x in self.bc])
        if not self.plain:                                   # the arrays (mat, lin) are stored beside the spec
            s['potentials'] = [p.spec() for p in self.pots]
        if self.mean is not None:
            s['mean'] = str(self.mean)
        return s

    @staticmethod
    def from_spec(s, gval, lin=None, mat=None):
        W = [potential_from_spec(q) for q in s['potentials']] if 'potentials' in s else [Fr(c) for c in s['W']]
        return Problem(box=[(Fr(a), Fr(b)) for a, b in s['box']], N=s['N'], k=Fr(s['k']), W=W, bc=s['bc'],
                       name=s.get('name', ''), gval=gval, lin=lin, mat=mat,
                       mean=(Fr(s['mean']) if s.get('mean') is not None else None))


def is_superlinear(W):
    """even degree >= 2 with positive leading coefficient"""
    W = list(W)
    while len(W) > 1 and W[-1] == 0:
        W.pop()
    d = len(W) - 1
    return d >= 2 and d % 2 == 0 and W[-1] > 0


# ================================================================================================================
# heuristic: a candidate (minimiser of the Q1 energy with 2-point Gauss quadrature for the potential)
# ================================================================================================================
def _shape_matrices(n):
    gp = [0.5 - 0.5 / math.sqrt(3.0), 0.5 + 0.5 / math.sqrt(3.0)]
    verts = list(itertools.product((0, 1), repeat=n))
    pts = list(itertools.product(gp, repeat=n))
    Nm = np.zeros((len(pts), len(verts)))
    Gm = np.zeros((n, len(pts), len(verts)))
    for p, xi in enumerate(pts):
        for a, v in enumerate(verts):
            w = [xi[j] if v[j] else 1.0 - xi[j] for j in range(n)]
            Nm[p, a] = math.prod(w)
            for i in range(n):
                Gm[i, p, a] = (1.0 if v[i] else -1.0) * math.prod(w[j] for j in range(n) if j != i)
    return Nm, Gm


def discrete_energy(P, U):
    n = P.n
    Nm, Gm = _shape_matrices(n)
    Uc = np.stack([U[P.vslice(v)] for v in P.verts], axis=-1)
    vol, k = float(P.cellvol), float(P.k)
    wq = vol / Nm.shape[0]
    E = 0.0
    gUc = np.zeros_like(Uc)
    for p in range(Nm.shape[0]):
        u = Uc @ Nm[p]
        E += wq * float(P.Wfloat(u).sum())
        gUc += wq * P.dWfloat(u)[..., None] * Nm[p]
        for i in range(n):
            d = (Uc @ Gm[i, p]) / P.hf[i]
            E += wq * k * float((d * d).sum())
            gUc += (wq * 2.0 * k / P.hf[i]) * d[..., None] * Gm[i, p]
    g = np.zeros_like(U)
    for a, v in enumerate(P.verts):
        g[P.vslice(v)] += gUc[..., a]
    return E, g


def solve_candidate(P, U0=None, maxiter=6000):
    """any nodal array is a valid candidate; this one is a local minimiser of the discrete energy started at U0"""
    U = np.zeros(P.fixed.shape) if U0 is None else np.array(U0, dtype=float) * np.ones(P.fixed.shape)
    U[P.fixed] = P.gval[P.fixed]
    free = ~P.fixed

    def fun(x):
        U[free] = x
        E, g = discrete_energy(P, U)
        return E, g[free]

    res = minimize(fun, U[free].copy(), jac=True, method='L-BFGS-B',
                   options=dict(maxiter=maxiter, maxfun=4 * maxiter, ftol=1e-15, gtol=1e-10, maxcor=30))
    U[free] = res.x
    return U


# ================================================================================================================
# RIGOROUS: energy of the Q1 competitor
# ================================================================================================================
def _fsum_up(a):
    """upper bound for the exact sum of the stored doubles: math.fsum is the correctly rounded sum under IEEE-754
    double arithmetic; it is widened by two units in the last place (one would do; the second covers builds on which
    fsum can be off in its last bit)"""
    return float(up(up(math.fsum(np.asarray(a, dtype=np.float64).ravel().tolist()))))


def _fsum_dn(a):
    return float(dn(dn(math.fsum(np.asarray(a, dtype=np.float64).ravel().tolist()))))


def upper_bound(P, U):
    """RIGOROUS enclosure (lo, hi) of E[u_h], u_h the Q1 function with nodal values U.
    Gradient term: the exact quadratic form of the nodal values, outward rounded.
    Potential term: midpoint rule on each cell with its remainder,
        | int_c F - |c| F(centre) | <= |c|/24 sum_i h_i^2 sup_c |d_ii F| ,   F = W(u_h) ,   d_ii F = W''(u_h) (d_i u_h)^2
    (u_h is multilinear, so d_ii u_h = 0; u_h(centre) is the mean of the vertex values; the ranges of u_h and of d_i u_h
    over a cell are spanned by the vertex values and by the edge difference quotients)."""
    n = P.n
    assert U.shape == P.fixed.shape and np.all(np.isfinite(U))
    assert np.array_equal(U[P.fixed], P.gval[P.fixed]), "Dirichlet nodal values are not imposed exactly"
    assert P.tied(U), "the nodal values on identified periodic nodes differ"
    Ulo, Uhi = nodal_enclosure(P, U)
    if P.domain is not None and not (Fr(float(Ulo.min())) >= P.domain[0] and Fr(float(Uhi.max())) <= P.domain[1]):
        raise ValueError("the candidate leaves the value range of the potential: its energy is infinite")
    Uc = {v: I(Ulo[P.vslice(v)], Uhi[P.vslice(v)]) for v in P.verts}
    invh = [fr(1 / h) for h in P.h]
    h2 = [fr(h * h) for h in P.h]
    umin = np.minimum.reduce([Uc[v].lo for v in P.verts])
    umax = np.maximum.reduce([Uc[v].hi for v in P.verts])
    W2rng = P.W2rng_cells(I(umin, umax))
    W2pm = I(np.minimum(W2rng.lo, 0.0), np.maximum(W2rng.hi, 0.0))
    if not P.plain:                               # W'' may be unbounded on a cell that touches a singular end
        fin = np.isfinite(W2rng.lo) & np.isfinite(W2rng.hi) & np.ones(P.N, dtype=bool)
        W2pm = I(np.where(fin, W2pm.lo, 0.0), np.where(fin, W2pm.hi, 0.0))
    kin = None
    rem_hi = np.zeros(P.N)
    rem_lo = np.zeros(P.N)
    corners = list(itertools.product((0, 1), repeat=n - 1))
    for i in range(n):
        d = {}
        for e in corners:
            d[e] = (Uc[e[:i] + (1,) + e[i:]] - Uc[e[:i] + (0,) + e[i:]]) * invh[i]
        acc = None
        for e in corners:                 # (1/|c|) int_c (d_i u_h)^2 = sum_{e,e'} d_e d_e' prod_j M[e_j, e'_j]
            for e2 in corners:
                w = math.prod((Fr(1, 3) if a == b else Fr(1, 6)) for a, b in zip(e, e2))
                term = (d[e].sqr() if e == e2 else d[e] * d[e2]) * fr(w)
                acc = term if acc is None else acc + term
        kin = acc if kin is None else kin + acc
        dmax2 = np.maximum.reduce([d[e].sqr().hi for e in corners])
        ri = W2pm * I(np.zeros(P.N), dmax2) * h2[i] * fr(Fr(1, 24))
        rem_hi = up(rem_hi + ri.hi)
        rem_lo = dn(rem_lo + ri.lo)
    kin = kin * fr(P.k)
    uc = None
    for v in P.verts:
        uc = Uc[v] if uc is None else uc + Uc[v]
    if P.plain:
        Wc = horner(P.W, uc * fr(Fr(1, 2 ** n)))
        cell_hi = up(up(kin.hi + Wc.hi) + rem_hi)
        cell_lo = dn(dn(kin.lo + Wc.lo) + rem_lo)
    else:                                          # Proposition 13.2 of docs/SCALAR_ND_2.md
        Wc = P.Wenc_cells(uc * fr(Fr(1, 2 ** n)))
        pot_hi = np.where(fin, up(Wc.hi + rem_hi), INF)                # midpoint rule where W'' is bounded on the cell
        pot_lo = np.where(fin, dn(Wc.lo + rem_lo), -INF)
        convex = (W2rng.lo >= 0) & np.ones(P.N, dtype=bool)            # vertex rule where W_c is convex on the range
        if np.any(convex):
            Wv = None
            for v in P.verts:
                t = P.Wenc_cells(Uc[v])
                Wv = t if Wv is None else Wv + t
            Wv = Wv * fr(Fr(1, 2 ** n))
            pot_hi = np.where(convex, np.minimum(pot_hi, Wv.hi), pot_hi)
            pot_lo = np.where(convex, np.maximum(pot_lo, Wc.lo), pot_lo)
        cell_hi = up(kin.hi + pot_hi)
        cell_lo = dn(kin.lo + pot_lo)
    if P.av is not None:                           # the term a(x) u: Proposition 13.1
        lin = None
        for v in P.verts:
            for v2 in P.verts:
                w = math.prod((Fr(1, 3) if a == b else Fr(1, 6)) for a, b in zip(v, v2))
                term = I(P.av[v]) * Uc[v2] * fr(w)
                lin = term if lin is None else lin + term
        cell_hi = up(cell_hi + lin.hi)
        cell_lo = dn(cell_lo + lin.lo)
    tot = I(_fsum_dn(cell_lo), _fsum_up(cell_hi)) * fr(P.cellvol)
    return float(tot.lo), float(tot.hi)


def nodal_enclosure(P, U):
    """The nodal values of the competitor as intervals (lo, hi).  Without a prescribed mean they are U itself.  With
    one, the competitor is u_h + delta phi_h, phi_h the Q1 function that is 1 at the nodes off the Dirichlet faces and
    0 on them, and delta the number that makes the mean exact (Proposition 13.4 of docs/SCALAR_ND_2.md); delta is
    enclosed, not computed exactly."""
    if P.mean is None:
        return U, U
    free = ~P.fixed
    phi = free.astype(float)

    def integral(X):                               # of the Q1 function with nodal values X: |c| sum over cells of vertex means
        acc = None
        for v in P.verts:
            acc = I(X[P.vslice(v)]) if acc is None else acc + I(X[P.vslice(v)])
        return I(_fsum_dn(acc.lo), _fsum_up(acc.hi)) * fr(P.cellvol / 2 ** P.n)

    Iphi = integral(phi)
    if not float(Iphi.lo) > 0.0:
        raise ValueError("the mean cannot be prescribed: every node lies on a Dirichlet face")
    delta = (fr(P.mean * P.vol) - integral(U)) * Iphi.inv()
    return (np.where(free, dn(U + float(delta.lo)), U), np.where(free, up(U + float(delta.hi)), U))


# ================================================================================================================
# the calibration field
# ================================================================================================================
class Calibration:
    """S(x, t) = t sigma0(x) + sum_m a_m(x) beta_m(clamp(t, tlo, thi)).
    S0[i]: normal component of sigma0 on the faces of direction i (shape face_shape(i)).
    C[i] : B-spline coefficients of the spline part on those faces (shape face_shape(i) + (K + 2,)).
    K = 0 (tlo = thi): the affine part alone.  affine=False: the spline part alone."""
    freeze_affine = False

    def __init__(self, P, tlo, thi, K, affine=None):
        self.tlo, self.thi, self.K = Fr(tlo), Fr(thi), int(K)
        sup = is_superlinear(P.W)
        self.affine = sup if affine is None else (bool(affine) and sup)
        if self.K == 0:
            assert self.tlo == self.thi and self.affine, "K = 0 needs tlo = thi and a superlinear potential"
        else:
            assert self.thi > self.tlo
        nb = self.K + 2 if self.K > 0 else 0
        self.C = [np.zeros(P.face_shape(i) + (nb,)) for i in range(P.n)]
        self.S0 = [np.zeros(P.face_shape(i)) for i in range(P.n)]
        self.mask = [np.ones(P.face_shape(i)) for i in range(P.n)]          # natural faces carry no flux
        for i in range(P.n):
            for side in (0, 1):
                if P.bc[i][side] == 'N':
                    sl = [slice(None)] * P.n
                    sl[i] = 0 if side == 0 else -1
                    self.mask[i][tuple(sl)] = 0.0

    def pack(self):
        return np.concatenate([c.ravel() for c in self.C] + [s.ravel() for s in self.S0])

    def unpack(self, x):
        o = 0
        for i, c in enumerate(self.C):
            self.C[i] = x[o:o + c.size].reshape(c.shape) * self.mask[i][..., None]
            o += c.size
        for i, s in enumerate(self.S0):
            if not self.freeze_affine:
                self.S0[i] = (x[o:o + s.size].reshape(s.shape) * self.mask[i]) if self.affine else np.zeros(s.shape)
            o += s.size

    def size(self):
        return sum(c.size for c in self.C) + (sum(s.size for s in self.S0) if self.affine else 0)

    def copy(self, P):
        new = Calibration(P, self.tlo, self.thi, self.K, affine=self.affine)
        new.freeze_affine = self.freeze_affine
        new.C = [c.copy() for c in self.C]
        new.S0 = [s.copy() for s in self.S0]
        return new

    def respects_natural_faces(self):
        return all(np.all(c[m == 0.0] == 0.0) and np.all(s[m == 0.0] == 0.0)
                   for c, s, m in zip(self.C, self.S0, self.mask))


def spline_tables(K, tlo, thi, Qo):
    """per knot interval, Qo + 1 samples (both ends): sample points, value / first-derivative matrices of the K + 2
    B-splines, second-derivative matrix per knot interval, Greville abscissae (sum_m xi_m beta_m(t) = t on [tlo, thi])"""
    D = (thi - tlo) / K
    Q = K * (Qo + 1)
    tq = np.zeros(Q)
    B = np.zeros((Q, K + 2))
    dB = np.zeros((Q, K + 2))
    B2 = np.zeros((K, K + 2))
    for kap in range(K):
        B2[kap, kap:kap + 3] = np.array([1.0, -2.0, 1.0]) / D ** 2
        for j in range(Qo + 1):
            q = kap * (Qo + 1) + j
            th = j / Qo
            tq[q] = tlo + (kap + th) * D
            B[q, kap:kap + 3] = [(1 - th) ** 2 / 2, (-2 * th * th + 2 * th + 1) / 2, th * th / 2]
            dB[q, kap:kap + 3] = np.array([th - 1, 1 - 2 * th, th]) / D
    xi = tlo + (np.arange(K + 2) - 0.5) * D
    return tq, B, dB, B2, xi


def _basis_row_float(K, tlo, thi, t):
    D = (thi - tlo) / K
    z = np.clip((np.asarray(t, dtype=float) - tlo) / D, 0.0, K)
    kap = np.minimum(np.floor(z).astype(int), K - 1)
    th = z - kap
    row = np.zeros(z.shape + (K + 2,))
    idx = np.indices(z.shape)
    for off, val in ((0, (1 - th) ** 2 / 2), (1, (-2 * th * th + 2 * th + 1) / 2), (2, th * th / 2)):
        row[tuple(idx) + (kap + off,)] = val
    return row


def _dirichlet_faces(P):
    """for every Dirichlet side: (i, side, sign of the outward normal, nodal g values at the 2^(n-1) vertices of each
    boundary face, as arrays over the faces)"""
    out = []
    for i in range(P.n):
        for side in (0, 1):
            if P.bc[i][side] != 'D':
                continue
            sl = [slice(None)] * P.n
            sl[i] = 0 if side == 0 else -1
            G = P.gval[tuple(sl)]
            shp = tuple(Nj for j, Nj in enumerate(P.N) if j != i)
            vv = [G[tuple(slice(e, Nj + e) for e, Nj in zip(ev, shp))]
                  for ev in itertools.product((0, 1), repeat=P.n - 1)]
            out.append((i, side, -1.0 if side == 0 else 1.0, vv))
    return out


# ---- restricted convex conjugates  T(d) = sup { t d - W(t) : t >= t0 } (side +1)  or  t <= t0 (side -1) -------------
def _reflect(W):
    return [c * (-1) ** j for j, c in enumerate(W)]


def conj_float(Wf, d, t0, side):
    """float values and maximisers of the restricted conjugate (heuristic use)"""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    W1 = np.array([j * Wf[j] for j in range(1, len(Wf))])
    T = np.empty(d.shape)
    A = np.empty(d.shape)
    pv = np.polynomial.polynomial.polyval
    for ix in np.ndindex(d.shape):
        p = W1.copy()
        p[0] -= d[ix]
        r = np.roots(p[::-1])
        r = r[np.abs(r.imag) < 1e-9 * (1 + np.abs(r.real))].real
        r = r[(r - t0) * side > 0]
        cand = np.concatenate([[t0], r])
        val = cand * d[ix] - pv(cand, Wf)
        j = int(np.argmax(val))
        T[ix], A[ix] = val[j], cand[j]
    return T, A


class ConjFloat:
    """float table of T(d) with C^1 cubic Hermite interpolation (heuristic use: smooth for the optimiser)"""

    def __init__(self, Wf, t0, side, dmax=4.0, nd=2049):
        self.Wf, self.t0, self.side, self.nd = Wf, float(t0), side, nd
        self._build(dmax)

    def _build(self, dmax):
        self.dmax = dmax
        self.g = np.linspace(-dmax, dmax, self.nd)
        self.T, self.A = conj_float(self.Wf, self.g, self.t0, self.side)
        self.hh = self.g[1] - self.g[0]

    def __call__(self, d):
        m = float(np.max(np.abs(d))) if np.size(d) else 0.0
        while m > 0.98 * self.dmax:
            self._build(2.0 * self.dmax)
        j = np.clip(np.floor((d - self.g[0]) / self.hh).astype(int), 0, self.nd - 2)
        s = (d - self.g[j]) / self.hh
        y0, y1, m0, m1 = self.T[j], self.T[j + 1], self.A[j] * self.hh, self.A[j + 1] * self.hh
        val = (1 + 2 * s) * (1 - s) ** 2 * y0 + s * (1 - s) ** 2 * m0 + s * s * (3 - 2 * s) * y1 + s * s * (s - 1) * m1
        der = ((6 * s * s - 6 * s) * (y0 - y1) + (3 * s * s - 4 * s + 1) * m0 + (3 * s * s - 2 * s) * m1) / self.hh
        return val, der


def conj_upper(P, dvals, t0, side, m_fine=6000, m_coarse=1500):
    """RIGOROUS upper bounds for T(d) at the exact floats dvals.  W must be superlinear.
    On a subinterval [a, b] of width w:   t d - W(t) <= max of the two end values + w^2/8 sup_[a,b] (W'')^+ .
    Beyond R the function is non-increasing, because there  W'(t) >= t^(deg-2) (deg a_deg t - sum_{j<deg} j |a_j|) >= |d|."""
    W = [Fr(c) for c in P.W]
    while len(W) > 1 and W[-1] == 0:
        W.pop()
    deg = len(W) - 1
    assert deg >= 2 and deg % 2 == 0 and W[-1] > 0
    dvals = np.asarray(dvals, dtype=float)
    t0 = Fr(t0)
    if side < 0:
        W, dv, t0 = _reflect(W), -dvals, -t0
    else:
        dv = dvals
    dmax = Fr(float(np.max(np.abs(dv)))) if dv.size else Fr(0)
    R = max(Fr(1), t0, (sum(j * abs(W[j]) for j in range(1, deg)) + dmax) / (deg * W[deg])) + 1
    # where the maxima are (float; only used to place the fine part of the partition)
    Wf = np.array([float(c) for c in W])
    _, arg = conj_float(Wf, np.array([float(np.min(dv)), float(np.max(dv)), float(np.median(dv))]), float(t0), +1)
    span = max(float(arg.max() - arg.min()), 1e-3 * float(R - t0), 1e-6)
    a = max(t0, Fr(float(arg.min()) - 0.25 * span).limit_denominator(10 ** 9))
    b = min(R, Fr(float(arg.max()) + 0.25 * span).limit_denominator(10 ** 9))
    if not (t0 <= a < b <= R):
        a, b = t0, R
    W2 = poly_deriv(poly_deriv(W))
    best = np.full(dv.shape, -INF)
    for (l, r, m) in ((t0, a, m_coarse), (a, b, m_fine), (b, R, m_coarse)):
        if r <= l:
            continue
        w = (r - l) / m
        ts = [l + j * w for j in range(m + 1)]
        enc = [frac_enclosure(t) for t in ts]
        tl = np.array([e[0] for e in enc])
        th = np.array([e[1] for e in enc])
        Wlo = np.array([frac_enclosure(sum(c * t ** q for q, c in enumerate(W)))[0] for t in ts])
        W2hi = np.maximum(horner(W2, I(tl[:-1], th[1:])).hi, 0.0)
        corr = up(frac_enclosure(w * w / 8)[1] * W2hi)
        f = up((I(tl, th) * I(dv[..., None])).hi - Wlo)
        seg = up(np.maximum(f[..., :-1], f[..., 1:]) + corr)
        best = np.maximum(best, seg.max(axis=-1))
    return best


def conj_cells(P, d0, t0, side, nd=513):
    """RIGOROUS upper bound, cell by cell, for T over the interval d0 (an I array).  T is convex: on an interval it is at
    most the larger of its two end values, and between two table nodes it is at most the chord through (upper bounds
    of) the node values."""
    lo, hi = float(np.min(d0.lo)), float(np.max(d0.hi))
    pad = max(1e-9, 1e-6 * max(abs(lo), abs(hi), hi - lo))
    g = np.linspace(lo - pad, hi + pad, nd)
    assert g[0] < lo and g[-1] > hi and np.all(np.diff(g) > 0)
    Tg = conj_upper(P, g, t0, side)

    def at(p):
        j = np.clip(np.searchsorted(g, p, side='right') - 1, 0, nd - 2)
        assert np.all(g[j] <= p) and np.all(p <= g[j + 1])
        lam = (I(p) - I(g[j])) * (I(g[j + 1]) - I(g[j])).inv()
        lam = I(np.clip(lam.lo, 0.0, 1.0), np.clip(lam.hi, 0.0, 1.0))
        return (I(Tg[j]) + lam * (I(Tg[j + 1]) - I(Tg[j]))).hi

    return np.maximum(at(d0.lo), at(d0.hi))


def tail_inf(P, t0, side):
    """RIGOROUS lower bound for inf { W(t) : t >= t0 } (side +1) or t <= t0 (side -1); -inf if W is unbounded below"""
    key = (Fr(t0), side)
    cache = P.__dict__.setdefault('_tail_cache', {})
    if key not in cache:
        W = [Fr(c) for c in P.W]
        while len(W) > 1 and W[-1] == 0:
            W.pop()
        d = len(W) - 1
        if d == 0:
            cache[key] = frac_enclosure(W[0])[0]
        elif d % 2 == 1 or W[-1] <= 0:
            cache[key] = -INF
        else:
            cache[key] = -float(conj_upper(P, np.array([0.0]), t0, side)[0])
    return cache[key]


# ================================================================================================================
# heuristic optimiser for the calibration field
# ================================================================================================================
def optimise_calibration(P, cal, Qo=4, taus=(1e-2, 1e-3, 1e-4), iters=300, scale=1.0, verbose=False, curvature=True):
    """Maximise a smoothed copy of the certified functional
           J(S) = flux(S) - sum_cells |c| 2^-n sum_vertices sup_t Phi_c(v, t)
    with L-BFGS.  J is concave in S.  sup_t is replaced by tau * log-mean-exp over the samples (with the curvature
    allowance the verifier uses, so nothing can hide between samples) and the two tails."""
    n, K = P.n, cal.K
    tlo, thi = float(cal.tlo), float(cal.thi)
    k = float(P.k)
    vol = float(P.cellvol)
    wv = vol / 2 ** n
    hf = P.hf
    cells = P.N
    if K > 0:
        tq, B, dB, B2, xi = spline_tables(K, tlo, thi, Qo)
        Q = len(tq)
        kq = np.repeat(np.arange(K), Qo + 1)
        Vq = P.Wfloat(tq)
        w = (thi - tlo) / K / Qo
        cw = w * w / 8 if curvature else 0.0      # curvature=False is for stress tests only: it lets peaks hide
        fine = np.polynomial.polynomial.polyval(np.linspace(tlo, thi, K * Qo * 8 + 1), P.W2f)
        segmax = np.array([fine[8 * s:8 * s + 9].max() for s in range(K * Qo)])          # sup W'' per subinterval
        W2max = np.empty(Q)
        for kap in range(K):
            for j in range(Qo + 1):
                W2max[kap * (Qo + 1) + j] = max(segmax[kap * Qo + jj] for jj in (j - 1, j) if 0 <= jj < Qo)
    else:
        Q = 0
    if cal.affine:
        cl, ch = ConjFloat(P.Wf, tlo, -1), ConjFloat(P.Wf, thi, +1)
    else:
        tl_c, th_c = -tail_inf(P, cal.tlo, -1), -tail_inf(P, cal.thi, +1)
        if not (np.isfinite(tl_c) and np.isfinite(th_c)):
            raise ValueError("the potential is not bounded below outside the t-range: no finite lower bound")
    # the flux is linear in the coefficients: 2-point Gauss rule on every Dirichlet face
    FWC = [np.zeros_like(c) for c in cal.C]
    FW0 = [np.zeros_like(s) for s in cal.S0]
    gp = [0.5 - 0.5 / math.sqrt(3.0), 0.5 + 0.5 / math.sqrt(3.0)]
    corners = list(itertools.product((0, 1), repeat=n - 1))
    pts = list(itertools.product(gp, repeat=n - 1))
    for i, side, sign, vv in _dirichlet_faces(P):
        area = vol / hf[i]
        sl = [slice(None)] * n
        sl[i] = 0 if side == 0 else -1
        if K > 0:
            for xg in pts:
                gq = sum(math.prod(x if e else 1 - x for x, e in zip(xg, ev)) * vv[a] for a, ev in enumerate(corners))
                FWC[i][tuple(sl)] += sign * area / len(pts) * _basis_row_float(K, tlo, thi, gq)
        FW0[i][tuple(sl)] += sign * area * sum(vv) / len(vv)
    hist = []
    opt_affine = cal.affine and not cal.freeze_affine

    def fun(x, tau):
        cal.unpack(x)
        C, S0 = cal.C, cal.S0
        obj = 0.0
        d0 = np.zeros(cells)
        A0 = []
        for i in range(n):
            lo, hi = P.fslice(i, 0), P.fslice(i, 1)
            d0 += (S0[i][hi] - S0[i][lo]) / hf[i]
            A0.append((S0[i][lo] ** 2 / (4 * k), S0[i][hi] ** 2 / (4 * k)))
        if cal.affine:
            Tl, dTl = cl(d0)
            Th, dTh = ch(d0)
        else:
            Tl, dTl, Th, dTh = tl_c, 0.0, th_c, 0.0
        if K > 0:
            Ct = [C[i] + S0[i][..., None] * xi for i in range(n)]
            ds = [c @ dB.T for c in Ct]
            D = np.zeros(cells + (Q,))
            D2 = np.zeros(cells + (K,))
            A = []
            Dsl = np.zeros(cells)
            Dsh = np.zeros(cells)
            for i in range(n):
                lo, hi = P.fslice(i, 0), P.fslice(i, 1)
                s = Ct[i] @ B.T
                s2 = Ct[i] @ B2.T
                D += (s[hi] - s[lo]) / hf[i]
                D2 += (s2[hi] - s2[lo]) / hf[i]
                A.append((ds[i][lo] ** 2 / (4 * k), ds[i][hi] ** 2 / (4 * k)))
                el = (C[i][..., 0] + C[i][..., 1]) / 2
                eh = (C[i][..., K] + C[i][..., K + 1]) / 2
                Dsl += (el[hi] - el[lo]) / hf[i]
                Dsh += (eh[hi] - eh[lo]) / hf[i]
            carg = W2max - D2[..., kq]
            ind = carg > 0
            base = D + cw * np.where(ind, carg, 0.0) - Vq
            GD = np.zeros(cells + (Q,))
            GA = [[np.zeros(cells + (Q,)), np.zeros(cells + (Q,))] for _ in range(n)]
        else:
            Dsl = Dsh = 0.0
        Gtl = np.zeros(cells)
        Gth = np.zeros(cells)
        GA0 = [[np.zeros(cells), np.zeros(cells)] for _ in range(n)]
        Phi = np.empty(cells + (Q + 2,))
        for v in P.verts:
            a0 = sum(A0[i][v[i]] for i in range(n))
            if K > 0:
                Phi[..., :Q] = base
                for i in range(n):
                    Phi[..., :Q] += A[i][v[i]]
            Phi[..., Q] = Dsl + a0 + Tl
            Phi[..., Q + 1] = Dsh + a0 + Th
            mx = Phi.max(axis=-1, keepdims=True)
            ex = np.exp((Phi - mx) / tau)
            sm = ex.sum(axis=-1, keepdims=True)
            obj -= wv * float((mx[..., 0] + tau * np.log(sm[..., 0] / (Q + 2))).sum())
            ex *= (-wv) / sm
            Gtl += ex[..., Q]
            Gth += ex[..., Q + 1]
            for i in range(n):
                GA0[i][v[i]] += ex[..., Q] + ex[..., Q + 1]
            if K > 0:
                GD += ex[..., :Q]
                for i in range(n):
                    GA[i][v[i]] += ex[..., :Q]
        gC, gS = [], []
        Gd0 = Gtl * dTl + Gth * dTh
        if K > 0:
            G2 = (GD * np.where(ind, -cw, 0.0)).reshape(cells + (K, Qo + 1)).sum(axis=-1)
        for i in range(n):
            lo, hi = P.fslice(i, 0), P.fslice(i, 1)
            g0 = np.zeros(P.face_shape(i))
            if opt_affine:
                g0[hi] += Gd0 / hf[i]
                g0[lo] -= Gd0 / hf[i]
                g0[lo] += GA0[i][0] * S0[i][lo] / (2 * k)
                g0[hi] += GA0[i][1] * S0[i][hi] / (2 * k)
                g0 += FW0[i]
            obj += float((FW0[i] * S0[i]).sum())
            if K > 0:
                gs = np.zeros(P.face_shape(i) + (Q,))
                gs[hi] += GD / hf[i]
                gs[lo] -= GD / hf[i]
                gds = np.zeros_like(gs)
                gds[lo] += GA[i][0] * ds[i][lo] / (2 * k)
                gds[hi] += GA[i][1] * ds[i][hi] / (2 * k)
                gs2 = np.zeros(P.face_shape(i) + (K,))
                gs2[hi] += G2 / hf[i]
                gs2[lo] -= G2 / hf[i]
                gc = gs @ B + gds @ dB + gs2 @ B2
                if opt_affine:
                    g0 += gc @ xi
                gel = np.zeros(P.face_shape(i))
                gel[hi] += Gtl / hf[i]
                gel[lo] -= Gtl / hf[i]
                geh = np.zeros(P.face_shape(i))
                geh[hi] += Gth / hf[i]
                geh[lo] -= Gth / hf[i]
                gc[..., 0] += gel / 2
                gc[..., 1] += gel / 2
                gc[..., K] += geh / 2
                gc[..., K + 1] += geh / 2
                obj += float((FWC[i] * C[i]).sum())
                gc += FWC[i]
                gC.append((gc * cal.mask[i][..., None]).ravel())
            else:
                gC.append(np.zeros(0))
            gS.append((g0 * cal.mask[i]).ravel())
        return -obj, -np.concatenate(gC + gS)

    x = cal.pack()
    t0 = time.time()
    for tau in taus:
        res = minimize(fun, x, args=(tau * scale,), jac=True, method='L-BFGS-B',
                       options=dict(maxiter=iters, maxfun=3 * iters, ftol=1e-15, gtol=1e-12, maxcor=20))
        x = res.x
        hist.append((tau, float(-res.fun), int(res.nit)))
        if verbose:
            print(f"   calibration: tau {tau:.1e}  J_tau = {-res.fun:.8f}  iters {res.nit}  ({time.time() - t0:.1f} s)",
                  flush=True)
    cal.unpack(x)
    optimise_calibration.last_fun = fun
    return hist


# ================================================================================================================
# RIGOROUS: the lower bound for a given calibration field
# ================================================================================================================
def _spline_eval(Cb, K, tlo, thi, t):
    """RIGOROUS enclosure of the clamped spline s_f(t_f) for a stack of faces: Cb has shape faces + (K+2,), t is an I
    array over the faces (thin intervals).  z = (t - tlo)/Delta is clamped to [0, K]; the piece containing z.lo is
    evaluated on the part of the interval inside it and the next piece on the rest, if any."""
    D = (thi - tlo) / K
    z = (t - fr(tlo)) * fr(1 / D)
    zlo = np.clip(z.lo, 0.0, float(K))
    zhi = np.clip(z.hi, 0.0, float(K))
    assert np.all(zhi - zlo < 0.5), "_spline_eval is for thin intervals (at most two pieces are visited)"
    kap = np.minimum(np.floor(zlo).astype(int), K - 1)
    idx = np.indices(kap.shape)

    def piece(kk, th):
        c0, c1, c2 = (I(Cb[tuple(idx) + (kk + r,)]) for r in range(3))
        one = I(1.0)
        n0 = (one - th).sqr() * fr(Fr(1, 2))
        n2 = th.sqr() * fr(Fr(1, 2))
        n1 = one - n0 - n2                                           # the three weights sum to 1
        return c0 * n0 + c1 * n1 + c2 * n2

    th1 = I(np.clip(dn(zlo - kap), 0.0, 1.0), np.clip(up(zhi - kap), 0.0, 1.0))
    r = piece(kap, th1)
    over = (up(zhi - kap) > 1.0) & (kap < K - 1)                  # the last piece has no successor: z is clamped to K
    if np.any(over):
        th2 = I(np.zeros(kap.shape), np.clip(up(zhi - kap - 1.0), 0.0, 1.0))
        r2 = piece(np.minimum(kap + 1, K - 1), th2)
        r = I(np.where(over, np.minimum(r.lo, r2.lo), r.lo), np.where(over, np.maximum(r.hi, r2.hi), r.hi))
    return r


def flux_lower(P, cal, r=8):
    """RIGOROUS lower bound for  int_{Gamma_D} S(x, g_h(x)) . nu ,  g_h multilinear on each boundary face.
    Affine part: linear in g_h, and the face mean of g_h is the mean of its vertex values.
    Spline part: each face is cut into r^(n-1) equal sub-faces (g_h is multilinear on each).  On a sub-face with mean
    gbar and dev >= max |g_h - gbar|, the clamped spline s is Lipschitz (constant L1) and, inside [tlo, thi], C^1 with
    |s''| <= L2, so     mean s(g_h) >= s(gbar) - min( L1 dev , L2 dev^2 / 2 )     (Taylor at gbar; the first-order term
    integrates to zero).  If g_h is constant on a face the value is s(g) and the penalty is zero."""
    K, tlo, thi = cal.K, cal.tlo, cal.thi
    n = P.n
    total = []
    corners = list(itertools.product((0, 1), repeat=n - 1))
    for i, side, sign, vv in _dirichlet_faces(P):
        sl = [slice(None)] * n
        sl[i] = 0 if side == 0 else -1
        area = fr(P.cellvol / P.h[i])
        shp = vv[0].shape
        acc = None
        for a in vv:
            acc = I(a) if acc is None else acc + I(a)
        gbar = acc * fr(Fr(1, len(vv)))
        face_lo = (I(cal.S0[i][tuple(sl)]) * gbar * I(sign)).lo               # affine part
        if K > 0:
            D = (thi - tlo) / K
            invD = fr(1 / D)
            Cb = cal.C[i][tuple(sl)]
            c0, c1, c2 = Cb[..., :-2], Cb[..., 1:-1], Cb[..., 2:]
            d01 = (I(c1) - I(c0)) * invD                                      # s' at the left end of each interval
            d12 = (I(c2) - I(c1)) * invD                                      # s' at the right end
            s2 = (d12 - d01) * invD
            s2abs = np.maximum(np.abs(s2.lo), np.abs(s2.hi))                  # |s''| on each knot interval
            s1abs = np.maximum.reduce([np.abs(d01.lo), np.abs(d01.hi), np.abs(d12.lo), np.abs(d12.hi)])
            tl, th, Df = float(tlo), float(thi), float(D)
            kk = np.arange(K)
            const = np.ones(shp, dtype=bool)
            for a in vv[1:]:
                const &= (a == vv[0])
            rr = 1 if bool(np.all(const)) else r
            sub_lo = np.zeros(shp)
            for sub in itertools.product(range(rr), repeat=n - 1):
                sv = []                                                        # g_h at the sub-face vertices
                for ev in corners:
                    loc = [Fr(s_ + e_, rr) for s_, e_ in zip(sub, ev)]
                    val = None
                    for b, cv in enumerate(corners):
                        wgt = math.prod(x if c_ else 1 - x for x, c_ in zip(loc, cv))
                        if wgt != 0:
                            term = I(vv[b]) * fr(wgt)
                            val = term if val is None else val + term
                    sv.append(val)
                gmin = np.minimum.reduce([x.lo for x in sv])
                gmax = np.maximum.reduce([x.hi for x in sv])
                m = None
                for x in sv:
                    m = x if m is None else m + x
                m = m * fr(Fr(1, len(sv)))                                     # sub-face mean of g_h
                val = _spline_eval(Cb, K, tlo, thi, m) * I(sign)
                dev = np.maximum(up(gmax - m.lo), up(m.hi - gmin))
                ka = np.clip(np.floor((gmin - tl) / Df).astype(int) - 1, 0, K - 1)
                kb = np.clip(np.floor((gmax - tl) / Df).astype(int) + 1, 0, K - 1)
                inr = (kk >= ka[..., None]) & (kk <= kb[..., None])
                L2 = np.where(inr, s2abs, 0.0).max(axis=-1)
                L1 = np.where(inr, s1abs, 0.0).max(axis=-1)
                inside = (gmin >= float(up(tl))) & (gmax <= float(dn(th)))
                pen2 = (I(L2) * I(dev).sqr() * fr(Fr(1, 2))).hi
                pen1 = (I(L1) * I(dev)).hi
                pen = np.where(const, 0.0, np.where(inside, np.minimum(pen2, pen1), pen1))
                sub_lo = dn(sub_lo + dn(val.lo - pen))
            face_lo = dn(face_lo + (I(sub_lo) * fr(Fr(1, rr ** (n - 1)))).lo)
        total.append(float((I(_fsum_dn(face_lo), _fsum_up(face_lo)) * area).lo))
    return float(dn(math.fsum(total))) if total else 0.0


def lower_bound(P, cal, Qs=16, return_parts=False):
    """RIGOROUS lower bound  flux - int_Omega sup_t Phi  for the calibration field cal (any coefficients).

    Cell integral: for fixed t, Phi(x, t) is convex in each x_i on a cell (each normal component is affine in its own
    coordinate and constant in the others), so sup_t Phi is at most the multilinear interpolant of its vertex values
    and its integral is at most |c| times their mean.
    At a vertex, on a subinterval [a, b] (width w) of a knot interval, Phi = (D - W) + (convex in t), so
        sup Phi <= max(Phi(a), Phi(b)) + w^2/8 sup (W'' - D'')^+ .
    Outside [tlo, thi] the spline part is constant in t:  Phi = D_spline(end) + |sigma0|^2/(4k) + (t d0 - W(t)), and
    the supremum of the last bracket is a restricted convex conjugate (conj_cells).
    K = 0: Phi = t d0 + |sigma0(x)|^2/(4k) - W(t), so sup_t Phi = |sigma0(x)|^2/(4k) + W*(d0) and the cell integral of
    the first term is computed exactly (the square of an affine function)."""
    n, K = P.n, cal.K
    tlo, thi = cal.tlo, cal.thi
    assert cal.respects_natural_faces(), "the field must have zero normal component on the natural faces"
    invh = [fr(1 / h) for h in P.h]
    inv4k = fr(1 / (4 * P.k))
    S0 = [I(s) for s in cal.S0]
    if not cal.affine:
        assert all(np.all(s == 0.0) for s in cal.S0)
    d0 = None
    A0 = []
    for i in range(n):
        lo, hi = P.fslice(i, 0), P.fslice(i, 1)
        term = (S0[i][hi] - S0[i][lo]) * invh[i]
        d0 = term if d0 is None else d0 + term
        A0.append((S0[i][lo].sqr() * inv4k, S0[i][hi].sqr() * inv4k))
    if cal.affine:
        Tl = conj_cells(P, d0, tlo, -1)
        Th = conj_cells(P, d0, thi, +1)
    else:
        Tl = np.full(P.N, -tail_inf(P, tlo, -1))
        Th = np.full(P.N, -tail_inf(P, thi, +1))
    fl = flux_lower(P, cal)
    if K == 0:
        a0 = None
        third = fr(Fr(1, 3))
        for i in range(n):
            lo, hi = P.fslice(i, 0), P.fslice(i, 1)
            term = (S0[i][lo].sqr() + S0[i][lo] * S0[i][hi] + S0[i][hi].sqr()) * third * inv4k
            a0 = term if a0 is None else a0 + term
        Mc = up(a0.hi + np.maximum(Tl, Th))
        if not np.all(np.isfinite(Mc)):
            return (-INF, None) if return_parts else -INF
        intM = I(_fsum_dn(Mc), _fsum_up(Mc)) * fr(P.cellvol)
        lb = float(dn(fl - float(intM.hi)))
        return (lb, dict(flux_lo=fl, intM_hi=float(intM.hi))) if return_parts else lb
    Mup = {v: np.full(P.N, -INF) for v in P.verts}
    Dl = (thi - tlo) / K
    invD = fr(1 / Dl)
    w = Dl / Qs
    w28 = fr(w * w / 8)
    xi = [tlo + (m - Fr(1, 2)) * Dl for m in range(K + 2)]
    C = cal.C
    half = fr(Fr(1, 2))
    Dend = []
    for m0 in (0, K):                                                # spline part at the two ends of the t-range
        acc = None
        for i in range(n):
            lo, hi = P.fslice(i, 0), P.fslice(i, 1)
            ev = (I(C[i][..., m0]) + I(C[i][..., m0 + 1])) * half
            term = (ev[hi] - ev[lo]) * invh[i]
            acc = term if acc is None else acc + term
        Dend.append(acc)
    for kap in range(K):
        c = [[I(C[i][..., kap + r]) + S0[i] * fr(xi[kap + r]) for r in range(3)] for i in range(n)]
        d01 = [c[i][1] - c[i][0] for i in range(n)]
        d12 = [c[i][2] - c[i][1] for i in range(n)]
        D2 = None                                                    # D'' on this knot interval
        for i in range(n):
            s2 = (d12[i] - d01[i]) * invD * invD
            lo, hi = P.fslice(i, 0), P.fslice(i, 1)
            term = (s2[hi] - s2[lo]) * invh[i]
            D2 = term if D2 is None else D2 + term
        prev = None
        for j in range(Qs + 1):
            th = Fr(j, Qs)
            t = tlo + (kap + th) * Dl
            n0, n1, n2 = fr((1 - th) ** 2 / 2), fr((-2 * th * th + 2 * th + 1) / 2), fr(th * th / 2)
            a0, a1 = fr(1 - th), fr(th)
            Dv = None
            A = []
            for i in range(n):
                s = c[i][0] * n0 + c[i][1] * n1 + c[i][2] * n2
                ds = (d01[i] * a0 + d12[i] * a1) * invD
                lo, hi = P.fslice(i, 0), P.fslice(i, 1)
                term = (s[hi] - s[lo]) * invh[i]
                Dv = term if Dv is None else Dv + term
                q = ds.sqr() * inv4k
                A.append((q[lo], q[hi]))
            base = Dv - fr(sum(cf * t ** q for q, cf in enumerate(P.W)))     # exact W(t), enclosed
            cur = {}
            for v in P.verts:
                ph = base
                for i in range(n):
                    ph = ph + A[i][v[i]]
                cur[v] = ph.hi
            if prev is not None:
                ta = frac_enclosure(tlo + (kap + Fr(j - 1, Qs)) * Dl)[0]
                tb = frac_enclosure(t)[1]
                W2hi = float(horner(P.W2, I(ta, tb)).hi)
                corr = ((I(W2hi) - D2).pos() * w28).hi
                for v in P.verts:
                    Mup[v] = np.maximum(Mup[v], up(np.maximum(prev[v], cur[v]) + corr))
            prev = cur
    for v in P.verts:
        a0 = None
        for i in range(n):
            a0 = A0[i][v[i]] if a0 is None else a0 + A0[i][v[i]]
        Mup[v] = np.maximum(Mup[v], (Dend[0] + a0 + I(Tl)).hi)
        Mup[v] = np.maximum(Mup[v], (Dend[1] + a0 + I(Th)).hi)
    allM = np.concatenate([Mup[v].ravel() for v in P.verts])
    if not np.all(np.isfinite(allM)):
        return (-INF, None) if return_parts else -INF
    intM = I(_fsum_dn(allM), _fsum_up(allM)) * fr(P.cellvol / 2 ** n)
    lb = float(dn(fl - float(intM.hi)))
    if return_parts:
        return lb, dict(flux_lo=fl, intM_hi=float(intM.hi), Mmax=float(allM.max()), Mmin=float(allM.min()))
    return lb


# ================================================================================================================
# heuristic: the affine part from a cell-centred finite-volume solve
# ================================================================================================================
def fv_energy(P, p, bmeans):
    """sum_faces k |c| (dp/h)^2 + Dirichlet faces 2k |c| ((g - p)/h)^2 + sum_cells |c| W(p), and its gradient.
    With the vertex rule, the optimal affine part is the two-point flux field of a stationary point of this energy."""
    n = P.n
    k = float(P.k)
    vol = float(P.cellvol)
    E = vol * float(P.Wfloat(p).sum())
    g = vol * P.dWfloat(p)
    for i in range(n):
        a = [slice(None)] * n
        a[i] = slice(0, P.N[i] - 1)
        b = [slice(None)] * n
        b[i] = slice(1, P.N[i])
        d = (p[tuple(b)] - p[tuple(a)]) / P.hf[i]
        E += k * vol * float((d * d).sum())
        g[tuple(b)] += 2 * k * vol * d / P.hf[i]
        g[tuple(a)] -= 2 * k * vol * d / P.hf[i]
    for i, side, gb in bmeans:
        c = [slice(None)] * n
        c[i] = 0 if side == 0 else P.N[i] - 1
        d = (p[tuple(c)] - gb) / P.hf[i]
        E += 2 * k * vol * float((d * d).sum())
        g[tuple(c)] += 4 * k * vol * d / P.hf[i]
    return E, g


def fv_solve(P, p0, maxiter=20000):
    bm = [(i, side, sum(vv) / len(vv)) for i, side, sign, vv in _dirichlet_faces(P)]
    res = minimize(lambda x: (lambda E, g: (E, g.ravel()))(*fv_energy(P, x.reshape(P.N), bm)), p0.ravel().copy(),
                   jac=True, method='L-BFGS-B',
                   options=dict(maxiter=maxiter, maxfun=4 * maxiter, ftol=1e-16, gtol=1e-12, maxcor=30))
    return res.x.reshape(P.N)


def affine_from_candidate(P, cal, U, solve=True):
    """sigma0 . e_i = 2k (p_R - p_L)/h_i on interior faces, 2k (g - p)/(h_i/2) on Dirichlet faces, 0 on natural faces,
    p = finite-volume solution started from the cell means of the candidate"""
    assert cal.affine
    n = P.n
    p = sum(U[P.vslice(v)] for v in P.verts) / 2 ** n
    if solve:
        p = fv_solve(P, p)
    k = float(P.k)
    bm = {(i, side): sum(vv) / len(vv) for i, side, sign, vv in _dirichlet_faces(P)}
    for i in range(n):
        s = np.zeros(P.face_shape(i))
        inner = [slice(None)] * n
        inner[i] = slice(1, P.N[i])
        a = [slice(None)] * n
        a[i] = slice(0, P.N[i] - 1)
        b = [slice(None)] * n
        b[i] = slice(1, P.N[i])
        s[tuple(inner)] = 2 * k * (p[tuple(b)] - p[tuple(a)]) / P.hf[i]
        for side in (0, 1):
            if P.bc[i][side] != 'D':
                continue
            f = [slice(None)] * n
            f[i] = 0 if side == 0 else -1
            c = [slice(None)] * n
            c[i] = 0 if side == 0 else P.N[i] - 1
            sgn = -1.0 if side == 0 else 1.0
            s[tuple(f)] = sgn * 2 * k * (bm[(i, side)] - p[tuple(c)]) / (P.hf[i] / 2)
        cal.S0[i] = s * cal.mask[i]
    return p


# ================================================================================================================
# heuristic: exact embeddings of a calibration field into finer spaces (multilevel)
# ================================================================================================================
def refine_t(P, cal):
    """the same field with 2K knot intervals (Chaikin subdivision of the B-spline coefficients)"""
    new = Calibration(P, cal.tlo, cal.thi, 2 * cal.K, affine=cal.affine)
    new.freeze_affine = cal.freeze_affine
    for i in range(P.n):
        c = cal.C[i]
        out = np.zeros(c.shape[:-1] + (2 * cal.K + 2,))
        out[..., 0::2] = (3 * c[..., :-1] + c[..., 1:]) / 4
        out[..., 1::2] = (c[..., :-1] + 3 * c[..., 1:]) / 4
        new.C[i] = out
        new.S0[i] = cal.S0[i].copy()
    return new


def prolong_x(Pc, cal, Pf):
    """the same field on the grid with twice as many cells in every direction (RT0 is nested under bisection)"""
    assert all(b == 2 * a for a, b in zip(Pc.N, Pf.N))
    new = Calibration(Pf, cal.tlo, cal.thi, cal.K, affine=cal.affine)
    new.freeze_affine = cal.freeze_affine

    def up_(a, i):
        for j in range(Pc.n):
            if j == i:
                sh = list(a.shape)
                sh[j] = 2 * Pc.N[j] + 1
                b = np.zeros(sh)
                ev = [slice(None)] * a.ndim
                ev[j] = slice(0, None, 2)
                od = [slice(None)] * a.ndim
                od[j] = slice(1, None, 2)
                lo = [slice(None)] * a.ndim
                lo[j] = slice(0, -1)
                hi = [slice(None)] * a.ndim
                hi[j] = slice(1, None)
                b[tuple(ev)] = a
                b[tuple(od)] = 0.5 * (a[tuple(lo)] + a[tuple(hi)])
                a = b
            else:
                a = np.repeat(a, 2, axis=j)
        return a

    for i in range(Pc.n):
        new.C[i] = up_(cal.C[i], i) * new.mask[i][..., None]
        new.S0[i] = up_(cal.S0[i], i) * new.mask[i]
    return new


# ================================================================================================================
# drivers
# ================================================================================================================
def certify_convex(P, U=None, U0=None):
    """convex mode: candidate -> upper bound; finite-volume solve -> affine field -> lower bound.  No optimiser."""
    T0 = time.time()
    if U is None:
        U = solve_candidate(P, U0)
    ub = upper_bound(P, U)
    cal = Calibration(P, 0, 0, 0)
    affine_from_candidate(P, cal, U)
    lb = lower_bound(P, cal)
    return dict(mode='convex', N=list(P.N), K=0, lower=lb, upper=ub[1], upper_enclosure=list(ub),
                rel_gap=(ub[1] - lb) / abs(ub[1]) if ub[1] != 0 else float('nan'), seconds=time.time() - T0,
                unknowns=cal.size()), U, cal


def certify_calibration(make_problem, levels, trange, U0=None, Qo=4, Qs=16,
                        first=((1e-2, 300), (1e-3, 400), (1e-4, 600)), later=((1e-4, 300), (3e-5, 400)),
                        last=((1e-4, 300), (3e-5, 600), (1e-5, 600)), scale=1.0, verbose=True, log=None):
    """calibration mode, coarse to fine.  levels: list of (N tuple, K); consecutive levels double every N_i and/or
    double K.  U0: function of the node coordinate arrays (start for the candidate).  At every level the field with the
    best certified value so far is kept, so the reported bound never decreases along a level."""
    cal = None
    Pprev = None
    rows = []
    T0 = time.time()
    for li, (N, K) in enumerate(levels):
        P = make_problem(N)
        if cal is None:
            cal = Calibration(P, trange[0], trange[1], K, affine=False)
        else:
            if tuple(N) != tuple(Pprev.N):
                cal = prolong_x(Pprev, cal, P)
            if K != cal.K:
                assert K == 2 * cal.K
                cal = refine_t(P, cal)
        best_lb, best = lower_bound(P, cal, Qs=Qs), cal.copy(P)
        for tau, it in (first if li == 0 else (last if li == len(levels) - 1 else later)):
            optimise_calibration(P, cal, Qo=Qo, taus=(tau,), iters=it, scale=scale)
            lb = lower_bound(P, cal, Qs=Qs)
            if lb > best_lb:
                best_lb, best = lb, cal.copy(P)
        cal = best.copy(P)
        X = np.meshgrid(*P.axes, indexing='ij')
        U = solve_candidate(P, None if U0 is None else U0(*X))
        ub = upper_bound(P, U)
        rows.append(dict(mode='calibration', N=list(N), K=K, lower=best_lb, upper=ub[1], upper_enclosure=list(ub),
                         rel_gap=(ub[1] - best_lb) / abs(ub[1]), seconds=time.time() - T0, unknowns=cal.size()))
        if verbose:
            msg = (f"   level {li}: N {tuple(N)} K {K}  {best_lb:.6f} <= inf E <= {ub[1]:.6f}  "
                   f"gap {100 * (ub[1] - best_lb) / abs(ub[1]):.3f} %  ({time.time() - T0:.0f} s)")
            print(msg, flush=True)
            if log:
                open(log, 'a').write(msg + '\n')
        Pprev = P
    return rows, U, cal, P


# ================================================================================================================
# certificates on disk:  <stem>.json (problem, claimed bounds)  +  <stem>.npz (exact floats: g, U, field coefficients)
# ================================================================================================================
def claim_strings(lower, upper, digits=10):
    """The two numbers a certificate claims, as exact decimal strings: the computed bounds weakened by one part in 10^9
    of the larger modulus and then rounded outward to `digits` significant decimal digits.
    The margin makes verification independent of last-bit differences between platforms: the float root finder that
    places the partition of the conjugate tables (conj_upper) is not bitwise reproducible.  Any partition gives a
    valid bound; only its last digits can move."""
    lo, hi = Fr(float(lower)), Fr(float(upper))                      # exact values of the doubles
    scale = max(abs(lo), abs(hi))
    if scale == 0:
        return '0', '0'
    ex = math.floor(math.log10(float(scale))) - digits + 1          # the claims are integers times 10^ex
    q = Fr(10) ** ex
    mlo = math.floor((lo - scale / 10 ** 9) / q)
    mhi = math.ceil((hi + scale / 10 ** 9) / q)

    def dec(m):                                                      # exact decimal string of m * 10^ex
        if ex >= 0:
            return str(m * 10 ** ex)
        sgn, a = ('-' if m < 0 else ''), str(abs(m)).rjust(-ex + 1, '0')
        return f'{sgn}{a[:ex]}.{a[ex:]}'

    return dec(mlo), dec(mhi)


def save_certificate(stem, P, U, cal, lower, upper, extra=None):
    lo, hi = claim_strings(lower, upper)
    meta = dict(format='CertiField-ND-certificate-1', problem=P.spec(), tlo=str(cal.tlo), thi=str(cal.thi), K=cal.K,
                affine=bool(cal.affine), lower=lo, upper=hi, reference=__version__)
    if extra:
        meta['info'] = extra
    json.dump(meta, open(stem + '.json', 'w'), indent=1)
    arrs = dict(gval=P.gval, U=U)
    for i in range(P.n):
        arrs[f'C{i}'] = cal.C[i]
        arrs[f'S0_{i}'] = cal.S0[i]
    np.savez_compressed(stem + '.npz', **arrs)


def load_certificate(stem):
    meta = json.load(open(stem + '.json'))
    z = np.load(stem + '.npz')
    P = Problem.from_spec(meta['problem'], z['gval'])
    cal = Calibration(P, Fr(meta['tlo']), Fr(meta['thi']), meta['K'], affine=meta['affine'])
    for i in range(P.n):
        assert z[f'C{i}'].shape == cal.C[i].shape and z[f'S0_{i}'].shape == cal.S0[i].shape
        cal.C[i] = np.array(z[f'C{i}'], dtype=float)
        cal.S0[i] = np.array(z[f'S0_{i}'], dtype=float)
    return meta, P, np.array(z['U'], dtype=float), cal


def verify_certificate(stem, verbose=True):
    """recompute both bounds from the stored exact data; True if they are at least as strong as the stored claims
    (the claims are exact decimals, the recomputed bounds exact doubles: the comparison is done in rationals)"""
    meta, P, U, cal = load_certificate(stem)
    lb = lower_bound(P, cal)
    ub = upper_bound(P, U)[1]
    ok = bool(math.isfinite(lb) and math.isfinite(ub) and Fr(lb) >= Fr(meta['lower']) and Fr(ub) <= Fr(meta['upper'])
              and Fr(meta['lower']) <= Fr(meta['upper']))
    if verbose:
        print(f"{stem}: recomputed {lb!r} <= inf E <= {ub!r};  claimed {meta['lower']} <= inf E <= {meta['upper']};  {'OK' if ok else 'FAILED'}")
    return ok, lb, ub


# ================================================================================================================
# raw export for the Wolfram Language port: little-endian float64, C order, plus field.json
# ================================================================================================================
def export_field_dir(d, P, U, cal, expected=None):
    import os
    os.makedirs(d, exist_ok=True)
    meta = dict(format='CertiField-ND-field-1', problem=P.spec(), TRange=[str(cal.tlo), str(cal.thi)], K=cal.K,
                Affine=bool(cal.affine), NodeDims=list(P.fixed.shape),
                FaceDims=[list(P.face_shape(i)) for i in range(P.n)])
    if expected is not None:
        meta['expected'] = expected
    json.dump(meta, open(os.path.join(d, 'field.json'), 'w'), indent=1)
    np.asarray(P.gval, dtype='<f8').tofile(os.path.join(d, 'gval.f64'))
    np.asarray(U, dtype='<f8').tofile(os.path.join(d, 'U.f64'))
    for i in range(P.n):
        np.asarray(cal.S0[i], dtype='<f8').tofile(os.path.join(d, f'S0_{i}.f64'))
        if cal.K > 0:
            np.asarray(cal.C[i], dtype='<f8').tofile(os.path.join(d, f'C{i}.f64'))
