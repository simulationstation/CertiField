#!/usr/bin/env python3
"""curved3d.py -- guaranteed two-sided stiffness bounds for a 3D linear-elastic body with a curved traction-free surface.

Setting (one octant of a body with three mirror planes):
    box (0,1)^3 on a uniform grid of N^3 cells, minus the "void"  kappa.x^2 < a^2  (sphere: kappa=(1,1,1)),
    material: isotropic, Young's modulus 1, Poisson ratio nu,
    faces x_k = 0: mirror planes;  faces y = 1, z = 1: traction-free;  the curved surface: traction-free;
    face x = 1: rigid grip moved by delta/2 in x  (clamped: u = (delta/2, 0, 0);  frictionless: u_x = delta/2 only).
    Full body = 8 mirror images; the stiffness K is defined by  minimal energy = K delta^2 / 2.

Upper bound:   K <= int_{Omega+} eps(v):C:eps(v)   for any admissible displacement v (unit delta), Omega+ a superset of Omega.
Lower bound:   K >= A^2 / int_{Omega+} tau:C^-1:tau  for any symmetric tau with div tau = 0, zero traction on the free
               surfaces and unit mean axial stress (A = full cross-section area):
               tau = tau_p + inc(rho^2 Psi),  rho a polynomial whose zero set is the curved surface, Psi a spline tensor.
Both integrands are polynomials on every box, so Gauss quadrature is exact; rounding is enclosed by interval arithmetic.
See NOTES.md for the argument, tests.py for the checks and verify.py for the independent re-check of stored fields.
"""
import itertools
import math
import time
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ia import I, fr, frac_enclosure, poly_deriv, sqr, hi_of, lo_of, sum_up, sum_dn, up, dn

# ---------------------------------------------------------------------------------------------------------------------
# index conventions
COMPS = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]          # xx yy zz xy yz xz
CIDX = {}
for _n, (_i, _j) in enumerate(COMPS):
    CIDX[(_i, _j)] = _n
    CIDX[(_j, _i)] = _n
VAR = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (0, 1, 1), (1, 0, 1)]
VIDX = {v: i for i, v in enumerate(VAR)}


def var_of(*axes):
    v = [0, 0, 0]
    for a in axes:
        v[a] += 1
    return tuple(v)


FIRST = [VIDX[var_of(k)] for k in range(3)]
SECOND = [VIDX[var_of(k, m)] for (k, m) in COMPS]


def _levi(i, j, k):
    return (i - j) * (j - k) * (k - i) // 2


def beltrami_tensor():
    """G[s, p, c]: tau_s = sum G[s,p,c] d_c Phi_p, from tau_ij = e_ikl e_jmn d_k d_m Phi_ln"""
    G = np.zeros((6, 6, 6))
    for s, (i, j) in enumerate(COMPS):
        for k, l, m, n in itertools.product(range(3), repeat=4):
            c = _levi(i, k, l) * _levi(j, m, n)
            if c:
                G[s, CIDX[(l, n)], CIDX[(k, m)]] += c
    return G


G = beltrami_tensor()
GS = np.zeros((6, 3, 3))            # engineering strain e_s = sum GS[s,a,k] d_k v_a
GS[0, 0, 0] = GS[1, 1, 1] = GS[2, 2, 2] = 1
GS[3, 0, 1] = GS[3, 1, 0] = 1
GS[4, 1, 2] = GS[4, 2, 1] = 1
GS[5, 0, 2] = GS[5, 2, 0] = 1


# ---------------------------------------------------------------------------------------------------------------------
# exact rational intervals (for the Gauss rule and the 1D tables)
class RI:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = Fr(lo)
        self.hi = self.lo if hi is None else Fr(hi)

    def __add__(self, o):
        o = _ri(o)
        return RI(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __sub__(self, o):
        o = _ri(o)
        return RI(self.lo - o.hi, self.hi - o.lo)

    def __rsub__(self, o):
        return _ri(o) - self

    def __mul__(self, o):
        o = _ri(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return RI(min(p), max(p))

    __rmul__ = __mul__


def _ri(o):
    return o if isinstance(o, RI) else RI(o)


def ri_horner(coefs, T):
    p = RI(coefs[-1])
    for c in reversed(coefs[:-1]):
        p = p * T + c
    return p


def legendre_coeffs(n):
    Pm, P = [Fr(1)], [Fr(0), Fr(1)]
    if n == 0:
        return Pm
    for k in range(1, n):
        a = [Fr(0)] + [Fr(2 * k + 1, k + 1) * c for c in P]
        b = [Fr(k, k + 1) * c for c in Pm] + [Fr(0)] * (len(a) - len(Pm))
        Pm, P = P, [u - v for u, v in zip(a, b)]
    return P


def gauss_enclosures(n, bits=120):
    """Gauss-Legendre rule on [0,1] with n points: rational enclosures of the nodes and of the weights"""
    P = legendre_coeffs(n)
    dP = poly_deriv(P)

    def val(c, x):
        r = Fr(0)
        for a in reversed(c):
            r = r * x + a
        return r

    M = 1024
    xs = [Fr(-1) + Fr(2 * k, M) for k in range(M + 1)]
    vs = [val(P, x) for x in xs]
    roots = []
    for k in range(M):
        if vs[k] == 0:
            roots.append((xs[k], xs[k]))
        elif vs[k + 1] != 0 and (vs[k] < 0) != (vs[k + 1] < 0):
            lo, hi, neg = xs[k], xs[k + 1], vs[k] < 0
            for _ in range(bits):
                mid = (lo + hi) / 2
                fm = val(P, mid)
                if fm == 0:
                    lo = hi = mid
                    break
                if (fm < 0) == neg:
                    lo = mid
                else:
                    hi = mid
            roots.append((lo, hi))
    assert len(roots) == n, (n, len(roots))
    nodes, weights = [], []
    for lo, hi in roots:
        X = RI(lo, hi)
        d = ri_horner(dP, X)
        den = (1 - X * X) * (d * d)
        assert den.lo > 0
        nodes.append(RI((1 + lo) / 2, (1 + hi) / 2))
        weights.append(RI(Fr(1) / den.hi, Fr(1) / den.lo))           # (2/den)/2
    return nodes, weights


class Enc:
    """float arrays lo <= exact <= hi and a nearby float mid"""

    def __init__(self, shape):
        self.lo = np.zeros(shape)
        self.hi = np.zeros(shape)
        self.mid = np.zeros(shape)

    def set(self, idx, r):
        self.lo[idx] = frac_enclosure(r.lo)[0]
        self.hi[idx] = frac_enclosure(r.hi)[1]
        self.mid[idx] = float((r.lo + r.hi) / 2)

    def iv(self, idx):
        return I(self.lo[idx], self.hi[idx])


SHAPES = [[Fr(1, 2), Fr(-1), Fr(1, 2)], [Fr(1, 2), Fr(1), Fr(-1)], [Fr(0), Fr(0), Fr(1, 2)]]   # quadratic B-spline on a cell


def poly_mul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r


def make_funcs(N, plateau, kcol=None, kflat=None):
    """1D functions as polynomials in the local coordinate t of cell e (physical x = (e+t)/N)"""
    kb1, kbm, kb2 = [Fr(b) * N for b in plateau]
    assert all(k.denominator == 1 for k in (kb1, kbm, kb2)) and 2 * kbm == kb1 + kb2, "plateau knots must be grid lines"
    kb1, kbm, kb2 = int(kb1), int(kbm), int(kb2)
    D = Fr(kb2 - kb1)

    def f_x(e):
        return [Fr(e, N), Fr(1, N)]

    def f_one(e):
        return [Fr(1)]

    def f_sq(e):
        return [Fr(e * e, N * N), Fr(2 * e, N * N), Fr(1, N * N)]

    def f_s(e):
        if e < kb1:
            return [Fr(1)]
        if e < kbm:
            u0 = Fr(e - kb1) / D
            return [1 - 2 * u0 * u0, -4 * u0 / D, -2 / (D * D)]
        if e < kb2:
            w0 = Fr(kb2 - e) / D
            return [2 * w0 * w0, -4 * w0 / D, 2 / (D * D)]
        return [Fr(0)]

    def f_ssq(e):
        return poly_mul(f_s(e), f_sq(e))

    kc = kb1 if kcol is None else int(kcol)

    def f_ind(e):
        return [Fr(1)] if e < kc else [Fr(0)]

    kf = N if kflat is None else int(kflat)

    def f_sqf(e):
        """eta(x) = x^2 for x <= b = kf/N, and x^2 - (x-b)^2/(1-b) beyond: C^1, increasing, eta'(1) = 0"""
        if e < kf:
            return f_sq(e)
        g = Fr(N, N - kf)
        return [Fr(e * e - g * (e - kf) ** 2, N * N), Fr(2 * e - 2 * g * (e - kf), N * N), Fr(1 - g, N * N)]

    return {'x': f_x, 'one': f_one, 'sq': f_sq, 'sqf': f_sqf, 's': f_s, 'ssq': f_ssq, 'ind': f_ind}


class Tables:
    """1D tables at the Gauss points of every dyadic sub-interval (level <= L) of a cell"""

    def __init__(self, N, L, nq, funcs):
        self.N, self.L, self.nq, self.funcs = N, L, nq, funcs
        nodes, weights = gauss_enclosures(nq)
        nLJ = 2 ** (L + 1) - 1
        self.bas = Enc((nLJ, 3, 3, nq))
        self.w = Enc((nLJ, nq))
        self.fn = {name: Enc((N, nLJ, 3, nq)) for name in funcs}
        for lev in range(L + 1):
            for j in range(2 ** lev):
                lj = 2 ** lev - 1 + j
                for q in range(nq):
                    T = RI((j + nodes[q].lo) / 2 ** lev, (j + nodes[q].hi) / 2 ** lev)
                    self.w.set((lj, q), RI(weights[q].lo / (N * 2 ** lev), weights[q].hi / (N * 2 ** lev)))
                    for m in range(3):
                        c = SHAPES[m]
                        for d in range(3):
                            self.bas.set((lj, d, m, q), ri_horner([x * N ** d for x in c], T))
                            c = poly_deriv(c)
                    for name, f in funcs.items():
                        for e in range(N):
                            c = f(e)
                            for d in range(3):
                                self.fn[name].set((e, lj, d, q), ri_horner([x * N ** d for x in c], T))
                                c = poly_deriv(c)


def float_tables_at(N, funcs, e, t):
    """float tables at arbitrary local points: e [n] ints, t [n] floats -> bas [n,3,3,1], fn[name] [n,3,1]"""
    e = np.asarray(e)
    t = np.asarray(t, dtype=float)
    n = len(t)
    bas = np.zeros((n, 3, 3, 1))
    for m in range(3):
        c = SHAPES[m]
        for d in range(3):
            cc = [float(x * N ** d) for x in c]
            bas[:, d, m, 0] = sum(cc[k] * t ** k for k in range(len(cc)))
            c = poly_deriv(c)
    fn = {name: np.zeros((n, 3, 1)) for name in funcs}
    for name, f in funcs.items():
        for i in range(n):
            c = f(int(e[i]))
            for d in range(3):
                fn[name][i, d, 0] = sum(float(x * N ** d) * t[i] ** k for k, x in enumerate(c))
                c = poly_deriv(c)
    return bas, fn


# ---------------------------------------------------------------------------------------------------------------------
# geometry: boxes of an octree per cell; exact rational classification
def build_boxes(N, L, a2, kappa=(1, 1, 1), solid_inside=False):
    """boxes covering a superset of Omega.  Returns int array [nbox, 8]: ex ey ez lev jx jy jz cut
    void = {kappa.x^2 < a^2} (solid_inside=False) or {kappa.x^2 > a^2} (solid_inside=True)"""
    a2 = Fr(a2)
    out = []

    def classify(lo, n):
        rmin = sum(kappa[k] * lo[k] * lo[k] for k in range(3))
        rmax = sum(kappa[k] * (lo[k] + 1) ** 2 for k in range(3))
        thr = a2 * n * n
        if solid_inside:
            if rmax <= thr:
                return 0
            if rmin >= thr:
                return 2
            return 1
        if rmin >= thr:
            return 0          # inside the solid
        if rmax <= thr:
            return 2          # inside the void
        return 1              # cut

    def rec(e, lev, j):
        lo = tuple(e[k] * 2 ** lev + j[k] for k in range(3))
        c = classify(lo, N * 2 ** lev)
        if c == 2:
            return
        if c == 0 or lev == L:
            out.append((*e, lev, *j, c))
            return
        for d in itertools.product((0, 1), repeat=3):
            rec(e, lev + 1, tuple(2 * j[k] + d[k] for k in range(3)))

    for e in itertools.product(range(N), repeat=3):
        rec(e, 0, (0, 0, 0))
    return np.array(out, dtype=np.int64)


def box_volume_bounds(boxes, N):
    """exact volumes: (superset volume, volume of the boxes proved inside)"""
    v_all = sum(Fr(1, (N * 2 ** int(l)) ** 3) for l in boxes[:, 3])
    v_in = sum(Fr(1, (N * 2 ** int(l)) ** 3) for l, c in zip(boxes[:, 3], boxes[:, 7]) if c == 0)
    return v_all, v_in


# ---------------------------------------------------------------------------------------------------------------------
# ties between spline coefficients (mirror planes, free flat faces, grips)
def expand1d(N, left, right):
    """E [N+2, nred]: full coefficients 0..N+1 from the reduced ones.
    left  = ('par', s): c0 = s c1          (mirror plane at 0;  s=+1 even, s=-1 odd / zero value)
    right = 'free' | 'zero2' (c_N = c_{N+1} = 0: value and normal derivative vanish) | ('par', s): c_{N+1} = s c_N"""
    if right == 'free':
        cols = list(range(1, N + 2))
    elif right == 'zero2':
        cols = list(range(1, N))
    else:
        cols = list(range(1, N + 1))
    pos = {j: i for i, j in enumerate(cols)}
    E = np.zeros((N + 2, len(cols)))
    for j in cols:
        E[j, pos[j]] = 1.0
    if 1 in pos:
        E[0, pos[1]] = float(left[1])
    if isinstance(right, tuple) and N in pos:
        E[N + 1, pos[N]] = float(right[1])
    return sp.csr_matrix(E), np.array(cols, dtype=np.int64)


def tie_matrix(N, specs, owners=False):
    """specs: per component a list of three (left, right) pairs -> sparse P [ncomp (N+2)^3, nred]
    (every row of P has at most one entry, +1 or -1).  With owners=True also the full index that carries each reduced one."""
    blocks, own = [], []
    M = N + 2
    for comp, sp3 in enumerate(specs):
        (Ex, cx), (Ey, cy), (Ez, cz) = [expand1d(N, l, r) for (l, r) in sp3]
        blocks.append(sp.kron(sp.kron(Ex, Ey), Ez, format='csr'))
        own.append((comp * M ** 3 + (cx[:, None, None] * M + cy[None, :, None]) * M + cz[None, None, :]).ravel())
    Pm = sp.block_diag(blocks, format='csr')
    return (Pm, np.concatenate(own)) if owners else Pm


def ties_hold(Pm, owners, cfull):
    """True if the stored coefficient array satisfies every tie EXACTLY (it is reproduced from its reduced part)"""
    flat = np.ascontiguousarray(cfull).reshape(-1)
    return bool(np.array_equal(Pm @ flat[owners], flat))


def psi_ties(N, grip='clamped', owners=False):
    specs = []
    for (a, b) in COMPS:
        sp3 = []
        for k in range(3):
            s = (-1) ** ((a == k) + (b == k))
            right = 'zero2' if k > 0 else ('free' if grip == 'clamped' else ('par', s))
            sp3.append((('par', s), right))
        specs.append(sp3)
    return tie_matrix(N, specs, owners)


def v_ties(N, grip='clamped', owners=False):
    specs = []
    for a in range(3):
        sp3 = []
        for k in range(3):
            s = -1 if a == k else 1
            if k == 0:
                right = ('par', -1) if (a == 0 or grip == 'clamped') else 'free'
            else:
                right = 'free'
            sp3.append((('par', s), right))
        specs.append(sp3)
    return tie_matrix(N, specs, owners)


# ---------------------------------------------------------------------------------------------------------------------
# particular Beltrami potential: sum of separable terms (component, coefficient, fx, fy, fz)
_ZERO = {('one', 1), ('one', 2), ('ind', 1), ('ind', 2)}


def _coef(c, rigorous):
    return fr(Fr(c)) if rigorous else float(c)


def particular_tau(pterms, sterms, FX, FY, FZ, rigorous=False):
    """stress components of inc(Phi_p) + direct stress terms; FX[name][d] etc. broadcast over the point grid.
    pterms: (potential component, coefficient, fx, fy, fz);  sterms: (stress component, coefficient, fx, fy, fz)"""
    Hp = {}
    for (p, coef, nx, ny, nz) in pterms:
        for c, (k, m) in enumerate(COMPS):
            dv = var_of(k, m)
            if (nx, dv[0]) in _ZERO or (ny, dv[1]) in _ZERO or (nz, dv[2]) in _ZERO:
                continue
            val = (FX[nx][dv[0]] * FY[ny][dv[1]]) * FZ[nz][dv[2]] * _coef(coef, rigorous)
            Hp[(p, c)] = val if (p, c) not in Hp else Hp[(p, c)] + val
    tau = [None] * 6
    for s in range(6):
        for (p, c), val in Hp.items():
            g = G[s, p, c]
            if g == 0:
                continue
            term = val * float(g)
            tau[s] = term if tau[s] is None else tau[s] + term
    for (s, coef, nx, ny, nz) in sterms:
        term = (FX[nx][0] * FY[ny][0]) * FZ[nz][0] * _coef(coef, rigorous)
        tau[s] = term if tau[s] is None else tau[s] + term
    return tau


def tension_terms():
    """potential (1 - w) diag(0, z^2/4, y^2/4), w = s(x) s(y) s(z): unit tension where w = 0, zero where w = 1"""
    q = Fr(1, 4)
    return [(1, q, 'one', 'one', 'sq'), (1, -q, 's', 's', 'ssq'), (2, q, 'one', 'sq', 'one'), (2, -q, 's', 'ssq', 's')], []


def uniform_terms():
    """uniform unit axial stress (admissible only if no free surface has a normal with an x-component)"""
    return [], [(0, Fr(1), 'one', 'one', 'one')]


def column_ramp_terms(plateau):
    """direct axial stress c (1 - s(y) s(z)) e_x e_x: zero on the column through the cavity, mean 1 over the section"""
    b1, bm, b2 = [Fr(b) for b in plateau]
    m1 = b1 + (b2 - b1) / 2
    c = 1 / (1 - m1 * m1)
    return [], [(0, c, 'one', 'one', 'one'), (0, -c, 'one', 's', 's')]


def column_step_terms(b):
    """direct axial stress c (1 - 1_{y<b} 1_{z<b}) e_x e_x"""
    b = Fr(b)
    c = 1 / (1 - b * b)
    return [], [(0, c, 'one', 'one', 'one'), (0, -c, 'one', 'ind', 'ind')]


# ---------------------------------------------------------------------------------------------------------------------
class Problem:
    def __init__(self, N, a=Fr(1, 2), nu=Fr(3, 10), kappa=(1, 1, 1), solid_inside=False, grip='clamped',
                 particular='step', column=Fr(1, 2), plateau=(Fr(1, 2), Fr(3, 4), Fr(1)), area=Fr(4), L_max=4, kflat=None):
        """particular: the force-carrying stress field
             'step'      axial stress outside the square column |y|,|z| < column     (cavity inside the column)
             'ramp'      the same with a C^1 ramp given by plateau
             'potential' Beltrami potential cut off by plateau
             'uniform'   uniform axial stress (only if no free surface has a normal with an x-component)"""
        self.N, self.a, self.nu, self.kappa = N, Fr(a), Fr(nu), tuple(kappa)
        self.solid_inside, self.grip, self.particular = solid_inside, grip, particular
        self.a2 = self.a * self.a
        column = Fr(column)
        kcol = column * N
        plateau = tuple(Fr(b) for b in plateau)
        assert kcol.denominator == 1, "the column must end on a grid line"
        self.funcs = make_funcs(N, plateau, int(kcol), kflat)
        self.eta = [('sqf' if (k == 0 and kflat is not None) else 'sq') if self.kappa[k] else None for k in range(3)]
        self.pterms, self.sterms = {'step': lambda: column_step_terms(column), 'ramp': lambda: column_ramp_terms(plateau),
                                    'potential': tension_terms, 'uniform': uniform_terms}[particular]()
        # full cross-section area carrying the unit axial stress: a Fraction or a pair (lower, upper) of Fractions
        self.area = (Fr(area), Fr(area)) if not isinstance(area, tuple) else (Fr(area[0]), Fr(area[1]))
        self.M = N + 2
        self.L_max = L_max
        self._tabs = {}
        self.lam = self.nu / ((1 + self.nu) * (1 - 2 * self.nu))
        self.mu = 1 / (2 * (1 + self.nu))
        self._check_admissible(column, plateau, kflat)

    def _check_admissible(self, column, plateau, kflat):
        """conditions under which the constructed stress field is admissible; the bounds are claimed only if they hold"""
        a, N, kap = self.a, self.N, self.kappa
        assert self.grip in ('clamped', 'frictionless') and 0 < self.nu < Fr(1, 2) and 0 <= a < 1
        assert all(k in (0, 1) for k in kap)
        x_normal = bool(kap[0]) and a > 0          # does the curved surface have points whose normal has an x-component?
        if self.grip == 'frictionless' and kap[0]:
            assert kflat is not None, "frictionless grips: the weight must be flat at the grip (kflat)"
        if kflat is not None:
            assert 0 < kflat < N and Fr(kflat, N) >= a, "the flattening must start beyond the curved surface"
        if not self.solid_inside:
            assert self.area == (Fr(4), Fr(4)), "box cross-section: the transmitted force of the unit field is 4"
        if self.particular == 'uniform':
            assert not x_normal, "uniform axial stress is not traction-free on this surface"
        else:
            assert not self.solid_inside and kap[1] and kap[2], "column fields need a cavity that is bounded in y and z"
            if self.particular == 'step':
                assert a <= column < 1, "the column must contain the cavity"
            else:
                assert a <= plateau[0] < plateau[2] <= 1, "the plateau must contain the cavity"

    def tables(self, nq):
        if nq not in self._tabs:
            self._tabs[nq] = Tables(self.N, self.L_max, nq, self.funcs)
        return self._tabs[nq]

    def boxes(self, L):
        return build_boxes(self.N, L, self.a2, self.kappa, self.solid_inside)

    def compliance(self):
        nu = float(self.nu)
        D = (1 + nu) * np.diag([1.0, 1, 1, 2, 2, 2])
        D[:3, :3] -= nu
        return D

    def stiffness(self):
        lam, mu = float(self.lam), float(self.mu)
        C = np.diag([2 * mu] * 3 + [mu] * 3)
        C[:3, :3] += lam
        return C


def local_to_global(ex, ey, ez, ncomp, M):
    i = np.arange(3)
    base = ((ex[:, None, None, None] + i[None, :, None, None]) * M + (ey[:, None, None, None] + i[None, None, :, None])) * M \
        + (ez[:, None, None, None] + i[None, None, None, :])
    base = base.reshape(len(ex), 27)
    return (np.arange(ncomp)[None, :, None] * M ** 3 + base[:, None, :]).reshape(len(ex), ncomp * 27)


def _lj(bx):
    lev = bx[:, 3]
    off = 2 ** lev - 1
    return off + bx[:, 4], off + bx[:, 5], off + bx[:, 6]


# ---------------------------------------------------------------------------------------------------------------------
# float kernels (used only to FIND good fields; nothing here needs to be rigorous)
def basis_derivs(Bx, By, Bz, variants):
    nb, nq = Bx.shape[0], Bx.shape[-1]
    out = np.empty((nb, len(variants), 27, nq ** 3))
    for v, (dx, dy, dz) in enumerate(variants):
        out[:, v] = np.einsum('bim,bjn,bkp->bijkmnp', Bx[:, dx], By[:, dy], Bz[:, dz], optimize=True).reshape(nb, 27, nq ** 3)
    return out


def grid3(x, y, z):
    nb, nq = x.shape
    X = np.broadcast_to(x[:, :, None, None], (nb, nq, nq, nq)).reshape(nb, -1)
    Y = np.broadcast_to(y[:, None, :, None], (nb, nq, nq, nq)).reshape(nb, -1)
    Z = np.broadcast_to(z[:, None, None, :], (nb, nq, nq, nq)).reshape(nb, -1)
    return [X, Y, Z]


def dual_core_float(P, Bx, By, Bz, fnx, fny, fnz):
    """T [nb, nQ, 6, 162] (stress of each local Psi-coefficient) and taup [nb, nQ, 6]"""
    nb, nq = Bx.shape[0], Bx.shape[-1]
    nQ = nq ** 3
    Bd = basis_derivs(Bx, By, Bz, VAR)
    fns = (fnx, fny, fnz)
    E = [None] * 3
    for k in range(3):
        if P.eta[k] is not None:
            E[k] = []
            for d in range(3):
                a3 = [np.ones((nb, nq))] * 3
                a3[k] = fns[k][P.eta[k]][:, d]
                E[k].append(grid3(*a3)[k])
    rho = -float(P.a2)
    for k in range(3):
        if E[k] is not None:
            rho = rho + E[k][0]
    R = rho * rho
    dR = [2 * rho * E[k][1] if E[k] is not None else None for k in range(3)]
    H = np.empty((nb, 27, 6, nQ))
    for c, (k, m) in enumerate(COMPS):
        h = R[:, None, :] * Bd[:, SECOND[c]]
        if dR[k] is not None:
            h = h + dR[k][:, None, :] * Bd[:, FIRST[m]]
        if dR[m] is not None:
            h = h + dR[m][:, None, :] * Bd[:, FIRST[k]]
        if E[k] is not None and E[m] is not None:
            d2 = 2 * E[k][1] * E[m][1]
            if k == m:
                d2 = d2 + 2 * rho * E[k][2]
            h = h + d2[:, None, :] * Bd[:, 0]
        H[:, :, c, :] = h
    T = np.einsum('spk,bikQ->bQspi', G, H, optimize=True).reshape(nb, nQ, 6, 162)
    FX = {n: [fnx[n][:, d][:, :, None, None] for d in range(3)] for n in fnx}
    FY = {n: [fny[n][:, d][:, None, :, None] for d in range(3)] for n in fny}
    FZ = {n: [fnz[n][:, d][:, None, None, :] for d in range(3)] for n in fnz}
    tp = particular_tau(P.pterms, P.sterms, FX, FY, FZ, False)
    taup = np.zeros((nb, nQ, 6))
    for s in range(6):
        if tp[s] is not None:
            taup[:, :, s] = np.broadcast_to(tp[s], (nb, nq, nq, nq)).reshape(nb, nQ)
    return T, taup


def primal_core_float(P, Bx, By, Bz):
    nb, nq = Bx.shape[0], Bx.shape[-1]
    Bd = basis_derivs(Bx, By, Bz, [VAR[i] for i in FIRST])            # [nb, 3(k), 27, nQ]
    return np.einsum('sak,bkiQ->bQsai', GS, Bd, optimize=True).reshape(nb, nq ** 3, 6, 81)


def _gather_float(P, tb, bx):
    ljx, ljy, ljz = _lj(bx)
    Bx, By, Bz = tb.bas.mid[ljx], tb.bas.mid[ljy], tb.bas.mid[ljz]
    fnx = {n: tb.fn[n].mid[bx[:, 0], ljx] for n in tb.fn}
    fny = {n: tb.fn[n].mid[bx[:, 1], ljy] for n in tb.fn}
    fnz = {n: tb.fn[n].mid[bx[:, 2], ljz] for n in tb.fn}
    W = (tb.w.mid[ljx][:, :, None, None] * tb.w.mid[ljy][:, None, :, None] * tb.w.mid[ljz][:, None, None, :]).reshape(len(bx), -1)
    return Bx, By, Bz, fnx, fny, fnz, W


def _by_cell(bx, M, *arrs):
    """sum the local arrays of consecutive boxes that lie in the same cell (boxes are sorted by cell)"""
    key = (bx[:, 0] * M + bx[:, 1]) * M + bx[:, 2]
    start = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    return bx[start], [np.add.reduceat(a, start, axis=0) for a in arrs]


class _Coo:
    def __init__(self, n, flush=6_000_000):
        self.n, self.flush_at = n, flush
        self.r, self.c, self.v, self.count = [], [], [], 0
        self.A = sp.csr_matrix((n, n))

    def add(self, gidx, Aloc):
        nb, nl = gidx.shape
        self.r.append(np.repeat(gidx, nl, axis=1).ravel())
        self.c.append(np.tile(gidx, (1, nl)).ravel())
        self.v.append(Aloc.ravel())
        self.count += nb * nl * nl
        if self.count > self.flush_at:
            self.flush()

    def flush(self):
        if self.count:
            self.A = self.A + sp.coo_matrix((np.concatenate(self.v), (np.concatenate(self.r), np.concatenate(self.c))),
                                            shape=(self.n, self.n)).tocsr()
            self.r, self.c, self.v, self.count = [], [], [], 0
        return self.A


def assemble_dual(P, boxes, batch=24, verbose=False):
    """quadratic form of  int tau:C^-1:tau  over the boxes (octant), tau = taup + T c:  returns A (csr), b, a0"""
    tb = P.tables(7)
    Lc = np.linalg.cholesky(P.compliance())
    nfull = 6 * P.M ** 3
    coo = _Coo(nfull)
    b = np.zeros(nfull)
    a0 = 0.0
    order = np.lexsort((boxes[:, 2], boxes[:, 1], boxes[:, 0]))
    boxes = boxes[order]
    t0 = time.time()
    for s in range(0, len(boxes), batch):
        bx = boxes[s:s + batch]
        Bx, By, Bz, fnx, fny, fnz, W = _gather_float(P, tb, bx)
        T, taup = dual_core_float(P, Bx, By, Bz, fnx, fny, fnz)
        sw = np.sqrt(W)
        Tp = (np.einsum('ts,bQsn->bQtn', Lc.T, T, optimize=True) * sw[:, :, None, None]).reshape(len(bx), -1, 162)
        tpp = (np.einsum('ts,bQs->bQt', Lc.T, taup) * sw[:, :, None]).reshape(len(bx), -1)
        Aloc = np.matmul(Tp.transpose(0, 2, 1), Tp)
        bloc = np.einsum('bqn,bq->bn', Tp, tpp)
        a0 += float(np.sum(tpp * tpp))
        bxc, (Aloc, bloc) = _by_cell(bx, P.M, Aloc, bloc)
        gidx = local_to_global(bxc[:, 0], bxc[:, 1], bxc[:, 2], 6, P.M)
        coo.add(gidx, Aloc)
        np.add.at(b, gidx.ravel(), bloc.ravel())
        if verbose and (s // batch) % 20 == 0:
            print(f"    dual assembly {s}/{len(boxes)}  {time.time() - t0:.0f}s", flush=True)
    return coo.flush(), b, a0


def assemble_primal(P, boxes, batch=512):
    """int eps(v):C:eps(v) over the boxes (octant), v = (x/2) e_x + spline(c):  returns K (csr), f, k0"""
    tb = P.tables(3)
    C = P.stiffness()
    Lc = np.linalg.cholesky(C)
    e0 = np.array([0.5, 0, 0, 0, 0, 0])
    nfull = 3 * P.M ** 3
    coo = _Coo(nfull)
    f = np.zeros(nfull)
    k0 = 0.0
    boxes = boxes[np.lexsort((boxes[:, 2], boxes[:, 1], boxes[:, 0]))]
    for s in range(0, len(boxes), batch):
        bx = boxes[s:s + batch]
        Bx, By, Bz, fnx, fny, fnz, W = _gather_float(P, tb, bx)
        Bm = primal_core_float(P, Bx, By, Bz)
        sw = np.sqrt(W)
        Bp = (np.einsum('ts,bQsn->bQtn', Lc.T, Bm, optimize=True) * sw[:, :, None, None]).reshape(len(bx), -1, 81)
        e0p = np.broadcast_to((Lc.T @ e0)[None, None, :], (len(bx), W.shape[1], 6)) * sw[:, :, None]
        e0p = e0p.reshape(len(bx), -1)
        Kloc = np.matmul(Bp.transpose(0, 2, 1), Bp)
        floc = np.einsum('bqn,bq->bn', Bp, e0p)
        k0 += float(np.sum(e0p * e0p))
        bxc, (Kloc, floc) = _by_cell(bx, P.M, Kloc, floc)
        gidx = local_to_global(bxc[:, 0], bxc[:, 1], bxc[:, 2], 3, P.M)
        coo.add(gidx, Kloc)
        np.add.at(f, gidx.ravel(), floc.ravel())
    return coo.flush(), f, k0


def solve_reduced(A, rhs, Pt, dense_limit=7000, tol=1e-10, maxiter=20000, verbose=False):
    """minimise  c^T A c / 2 - rhs^T c  over c = Pt y  (A positive semi-definite); returns y"""
    Ar = (Pt.T @ A @ Pt).tocsr()
    br = Pt.T @ rhs
    d = Ar.diagonal()
    scale = d.max() if d.size else 1.0
    dead = d <= 1e-14 * scale
    if Ar.shape[0] <= dense_limit:
        Ad = Ar.toarray()
        Ad[dead, :] = 0
        Ad[:, dead] = 0
        Ad[dead, dead] = scale
        br = np.where(dead, 0.0, br)
        Ad[np.diag_indices_from(Ad)] += 1e-13 * scale
        y = np.linalg.solve(Ad, br)
        return y
    dd = np.where(dead, 1.0, d)
    br = np.where(dead, 0.0, br)
    keep = sp.diags((~dead).astype(float))
    Ak = keep @ Ar @ keep + sp.diags(dead.astype(float) * 1.0 + 1e-13 * scale)
    Mi = spla.LinearOperator(Ak.shape, matvec=lambda x: x / dd)
    its = [0]

    def cb(xk):
        its[0] += 1
        if verbose and its[0] % 500 == 0:
            print(f"    cg iteration {its[0]}", flush=True)

    y, info = spla.cg(Ak, br, rtol=tol, maxiter=maxiter, M=Mi, callback=cb)
    if verbose:
        print(f"    cg stopped after {its[0]} iterations, info={info}", flush=True)
    return y


# ---------------------------------------------------------------------------------------------------------------------
# rigorous evaluation (intervals).  The same code runs on plain floats when rigorous=False (for tests).
def _tab(enc, idx, rigorous):
    return I(enc.lo[idx], enc.hi[idx]) if rigorous else enc.mid[idx]


def _contract(C, BX, BY, BZ, variants):
    """sum-factorised spline derivatives at the tensor grid of points.
    C [nb, nc, 3,3,3] floats (exact); BX.. [nb, 3(d), 3(m), nq] (I or float) -> dict variant -> [nb, nc, nq,nq,nq]"""
    need1 = sorted({v[0] for v in variants})
    need2 = sorted({(v[0], v[1]) for v in variants})
    A1 = {}
    for d in need1:
        acc = None
        for m in range(3):
            term = BX[:, d, m, :][:, None, None, None, :] * C[:, :, m, :, :, None]
            acc = term if acc is None else acc + term
        A1[d] = acc                                              # [nb, nc, 3(i2), 3(i3), nq1]
    A2 = {}
    for (d1, d2) in need2:
        acc = None
        for m in range(3):
            term = A1[d1][:, :, m, :, :, None] * BY[:, d2, m, :][:, None, None, None, :]
            acc = term if acc is None else acc + term
        A2[(d1, d2)] = acc                                       # [nb, nc, 3(i3), nq1, nq2]
    A3 = {}
    for v in variants:
        acc = None
        for m in range(3):
            term = A2[(v[0], v[1])][:, :, m, :, :, None] * BZ[:, v[2], m, :][:, None, None, None, :]
            acc = term if acc is None else acc + term
        A3[v] = acc                                              # [nb, nc, nq1, nq2, nq3]
    return A3


def _local_coeffs(cfull, bx):
    """cfull [nc, M, M, M] -> [nb, nc, 3,3,3]"""
    i = np.arange(3)
    ix = bx[:, 0][:, None, None, None] + i[None, :, None, None]
    iy = bx[:, 1][:, None, None, None] + i[None, None, :, None]
    iz = bx[:, 2][:, None, None, None] + i[None, None, None, :]
    return np.moveaxis(cfull[:, ix, iy, iz], 0, 1)


def _weights_hi(tb, bx, rigorous):
    ljx, ljy, ljz = _lj(bx)
    if rigorous:
        w = up(up(tb.w.hi[ljx][:, :, None, None] * tb.w.hi[ljy][:, None, :, None]) * tb.w.hi[ljz][:, None, None, :])
    else:
        w = tb.w.mid[ljx][:, :, None, None] * tb.w.mid[ljy][:, None, :, None] * tb.w.mid[ljz][:, None, None, :]
    return w


def dual_stress(P, tb, bx, cfull, rigorous=True, scale_p=1.0):
    """the six stress components of  scale_p * inc(Phi_p) + inc(rho^2 Psi)  at the Gauss points of the boxes"""
    ljx, ljy, ljz = _lj(bx)
    BX, BY, BZ = _tab(tb.bas, ljx, rigorous), _tab(tb.bas, ljy, rigorous), _tab(tb.bas, ljz, rigorous)
    fx = {n: _tab(tb.fn[n], (bx[:, 0], ljx), rigorous) for n in tb.fn}
    fy = {n: _tab(tb.fn[n], (bx[:, 1], ljy), rigorous) for n in tb.fn}
    fz = {n: _tab(tb.fn[n], (bx[:, 2], ljz), rigorous) for n in tb.fn}
    FX = {n: [fx[n][:, d][:, :, None, None] for d in range(3)] for n in fx}
    FY = {n: [fy[n][:, d][:, None, :, None] for d in range(3)] for n in fy}
    FZ = {n: [fz[n][:, d][:, None, None, :] for d in range(3)] for n in fz}
    tau = [None] * 6
    if cfull is not None:
        A3 = _contract(_local_coeffs(cfull, bx), BX, BY, BZ, VAR)
        f3 = (fx, fy, fz)
        E = [None] * 3
        for k in range(3):
            if P.eta[k] is not None:
                tabk = f3[k][P.eta[k]]
                sl = [None, None, None, None]
                sl[1 + k] = slice(None)
                E[k] = [tabk[:, d][(slice(None),) + tuple(sl)] for d in range(3)]      # [nb, 1, nq at axis k, ...]
        a2 = fr(P.a2) if rigorous else float(P.a2)
        rho = None
        for k in range(3):
            if E[k] is not None:
                rho = E[k][0] if rho is None else rho + E[k][0]
        rho = rho - a2
        R = sqr(rho)
        dR = [(rho * E[k][1]) * 2.0 if E[k] is not None else None for k in range(3)]
        Hphi = []
        for c, (k, m) in enumerate(COMPS):
            h = R * A3[VAR[SECOND[c]]]
            if dR[k] is not None:
                h = h + dR[k] * A3[VAR[FIRST[m]]]
            if dR[m] is not None:
                h = h + dR[m] * A3[VAR[FIRST[k]]]
            if E[k] is not None and E[m] is not None:
                d2 = (E[k][1] * E[m][1]) * 2.0
                if k == m:
                    d2 = d2 + (rho * E[k][2]) * 2.0
                h = h + d2 * A3[VAR[0]]
            Hphi.append(h)                                        # [nb, 6(p), nq,nq,nq]
        for s in range(6):
            for p in range(6):
                for c in range(6):
                    g = G[s, p, c]
                    if g == 0:
                        continue
                    term = Hphi[c][:, p] * float(g)
                    tau[s] = term if tau[s] is None else tau[s] + term
    if scale_p != 0.0:
        tp = particular_tau(P.pterms, P.sterms, FX, FY, FZ, rigorous)
        for s in range(6):
            if tp[s] is not None:
                term = tp[s] if scale_p == 1.0 else tp[s] * float(scale_p)
                tau[s] = term if tau[s] is None else tau[s] + term
    return tau


def dual_compliance(P, boxes, cfull, L=None, rigorous=True, batch=256, verbose=False, nq=7):
    """upper bound (rigorous=True) of  int_{boxes} tau:C^-1:tau  for tau = inc(Phi_p + rho^2 Psi), one octant.
    The integrand has degree <= 12 in each variable on every box, so the 7-point Gauss rule (exact to degree 13) is exact."""
    tb = P.tables(nq)
    nu = fr(P.nu) if rigorous else float(P.nu)
    tot = []
    t0 = time.time()
    for s0 in range(0, len(boxes), batch):
        bx = boxes[s0:s0 + batch]
        tau = dual_stress(P, tb, bx, cfull, rigorous)
        z = I(np.zeros(1)) if rigorous else 0.0
        t = [x if x is not None else z for x in tau]
        tr = t[0] + t[1] + t[2]
        q = (sqr(t[0]) + sqr(t[1]) + sqr(t[2]) + (sqr(t[3]) + sqr(t[4]) + sqr(t[5])) * 2.0) * (1 + nu) - sqr(tr) * nu
        w = _weights_hi(tb, bx, rigorous)
        if rigorous:
            qh = np.maximum(np.broadcast_to(q.hi, w.shape), 0.0)
            tot.append(sum_up(up(qh * w).reshape(len(bx), -1), axis=1))
        else:
            tot.append(np.sum((np.broadcast_to(q, w.shape) * w).reshape(len(bx), -1), axis=1))
        if verbose and (s0 // batch) % 10 == 0:
            print(f"    dual certificate {s0}/{len(boxes)}  {time.time() - t0:.0f}s", flush=True)
    tot = np.concatenate(tot)
    return float(sum_up(tot)) if rigorous else float(np.sum(tot))


def primal_strain(P, tb, bx, cfull, rigorous=True):
    ljx, ljy, ljz = _lj(bx)
    BX, BY, BZ = _tab(tb.bas, ljx, rigorous), _tab(tb.bas, ljy, rigorous), _tab(tb.bas, ljz, rigorous)
    A3 = _contract(_local_coeffs(cfull, bx), BX, BY, BZ, [VAR[i] for i in FIRST])
    dv = [[A3[VAR[FIRST[k]]][:, a] for k in range(3)] for a in range(3)]       # dv[a][k] = d_k v_a
    return [dv[0][0] + 0.5, dv[1][1], dv[2][2], dv[0][1] + dv[1][0], dv[1][2] + dv[2][1], dv[0][2] + dv[2][0]]


def primal_energy(P, boxes, cfull, rigorous=True, batch=2048, nq=3):
    """upper bound (rigorous=True) of  int_{boxes} eps(v):C:eps(v), v = (x/2) e_x + spline, one octant.
    The integrand has degree <= 4 in each variable on every box, so the 3-point Gauss rule (exact to degree 5) is exact."""
    tb = P.tables(nq)
    lam = fr(P.lam) if rigorous else float(P.lam)
    mu = fr(P.mu) if rigorous else float(P.mu)
    tot = []
    for s0 in range(0, len(boxes), batch):
        bx = boxes[s0:s0 + batch]
        e = primal_strain(P, tb, bx, cfull, rigorous)
        tr = e[0] + e[1] + e[2]
        q = sqr(tr) * lam + (sqr(e[0]) + sqr(e[1]) + sqr(e[2])) * (mu * 2.0) + (sqr(e[3]) + sqr(e[4]) + sqr(e[5])) * mu
        w = _weights_hi(tb, bx, rigorous)
        if rigorous:
            qh = np.maximum(np.broadcast_to(q.hi, w.shape), 0.0)
            tot.append(sum_up(up(qh * w).reshape(len(bx), -1), axis=1))
        else:
            tot.append(np.sum((np.broadcast_to(q, w.shape) * w).reshape(len(bx), -1), axis=1))
    tot = np.concatenate(tot)
    return float(sum_up(tot)) if rigorous else float(np.sum(tot))


# ---------------------------------------------------------------------------------------------------------------------
def certify(P, L_solve=2, L_cert=3, verbose=True, maxiter=3000, Lp_solve=None, Lp_cert=None):
    """rigorous bounds K_low <= K <= K_up for the full body (unit Young's modulus).
    L_solve / L_cert: octree depth used to find / to certify the stress field;  Lp_*: the same for the displacement"""
    t0 = time.time()
    out = {}
    Lp_solve = L_solve if Lp_solve is None else Lp_solve
    Lp_cert = L_cert if Lp_cert is None else Lp_cert
    cache = {}

    def boxes(L):
        if L not in cache:
            cache[L] = P.boxes(L)
        return cache[L]

    # ---- upper bound
    bs, bc = boxes(Lp_solve), boxes(Lp_cert)
    v_all, v_in = box_volume_bounds(bc, P.N)
    if verbose:
        print(f"  N={P.N} upper bound: {len(bs)} boxes to solve (L={Lp_solve}), {len(bc)} to certify (L={Lp_cert}); "
              f"octant volume between {float(v_in):.6f} and {float(v_all):.6f}", flush=True)
    K, f, k0 = assemble_primal(P, bs)
    Pv = v_ties(P.N, P.grip)
    y = solve_reduced(K, -f, Pv, verbose=verbose, maxiter=maxiter)
    cv = (Pv @ y).reshape(3, P.M, P.M, P.M)
    Kup_float = 8 * primal_energy(P, bc, cv, rigorous=False)
    Kup = float(np.nextafter(8 * primal_energy(P, bc, cv, rigorous=True), np.inf))
    out['K_up'], out['K_up_float'], out['v'] = Kup, Kup_float, cv
    if verbose:
        print(f"  upper bound  K <= {Kup:.8f}   (float value {Kup_float:.8f})   {time.time() - t0:.0f}s", flush=True)
    # ---- lower bound
    bs, bc = boxes(L_solve), boxes(L_cert)
    v_all, v_in = box_volume_bounds(bc, P.N)
    out['volume_superset'], out['volume_inside'] = float(v_all), float(v_in)
    if verbose:
        print(f"  N={P.N} lower bound: {len(bs)} boxes to solve (L={L_solve}), {len(bc)} to certify (L={L_cert}); "
              f"octant volume between {float(v_in):.6f} and {float(v_all):.6f}", flush=True)
    A, b, a0 = assemble_dual(P, bs, verbose=verbose)
    Pp = psi_ties(P.N, P.grip)
    y = solve_reduced(A, -b, Pp, verbose=verbose, maxiter=maxiter)
    del A
    cp = (Pp @ y).reshape(6, P.M, P.M, P.M)
    Q_float = 8 * dual_compliance(P, bc, cp, rigorous=False)
    Q_up = float(np.nextafter(8 * dual_compliance(P, bc, cp, rigorous=True, verbose=verbose), np.inf))
    alo = frac_enclosure(P.area[0])[0]
    Klow = float(dn(dn(alo * alo) / Q_up))
    out['K_low'], out['K_low_float'], out['psi'] = Klow, float(P.area[0]) ** 2 / Q_float, cp
    out['gap_percent'] = 100 * (Kup - Klow) / Klow
    if verbose:
        print(f"  lower bound  K >= {Klow:.8f}   (float value {out['K_low_float']:.8f})   {time.time() - t0:.0f}s", flush=True)
        print(f"  certified:  {Klow:.6f} <= K <= {Kup:.6f}    gap {out['gap_percent']:.3f} %", flush=True)
    return out
