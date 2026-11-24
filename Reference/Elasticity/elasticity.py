#!/usr/bin/env python3
"""CertiField elasticity -- reference implementation (Python/NumPy), two dimensions, plane strain.

Certified two-sided bounds for linear elasticity on a rectangle,

    E[u] = int_Omega  1/2 e(u)^T D e(u) - b.u dx ,     e(u) = (d1 u1, d2 u2, d1 u2 + d2 u1) ,
    u = g on the Dirichlet faces, no traction on the others.

docs/ELASTICITY.md has the theorem and its proof.

  upper bound : the energy of the continuous piecewise polynomial (Q_p, p = 1 or 2) displacement with the given nodal
                values, enclosed rigorously.
  lower bound : for ANY statically admissible stress field sigma (polynomial on each cell, continuous normal traction,
                div sigma + b = 0 exactly, zero traction on the traction-free faces exactly), every admissible u
                satisfies
                    E[u]  >=  int_{Gamma_D} (sigma nu).g  -  1/2 int_Omega s^T S s ,        S = D^-1 .
                The stress field is sigma = (d22 phi, d11 phi, -d12 phi) + (-b1 x1, -b2 x2, 0) with phi a C^1 piecewise
                bicubic function (four nodal values phi, d1 phi, d2 phi, d12 phi per node).  At the nodes of the
                traction-free faces the nodal values are those of a polynomial psi_b + sum_k P_k psi_k found in
                exact rational arithmetic, which makes the traction vanish there exactly.

Everything in upper_bound and lower_bound is rigorous: outward-rounded interval arithmetic (ia.py) on exact inputs
(binary floats and rationals).  The two linear solves (solve_candidate, solve_stress) are plain floating point and
are never trusted.

Both evaluators form their quadratic forms from DIFFERENCES of nodal values (from an affine function for the Airy
function, from a constant for the displacement).  The element matrices annihilate affine Airy functions and constant
displacements exactly (asserted in rational arithmetic), so nothing changes mathematically; but the products of the
nodal values themselves are larger than the result by the factor (|phi| / (h^2 |stress|))^2, and outward rounding of
terms that cancel would cost accuracy growing like h^-4 under refinement.
"""
import itertools
import json
import math
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ia import I, fr, dn, up

__version__ = "0.1.0-reference"

# ================================================================================================================
# exact one-dimensional tools: polynomials are ascending lists of Fractions
# ================================================================================================================
def pmul(p, q):
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, c in enumerate(q):
            r[i + j] += a * c
    return r


def pder(p, k=1):
    for _ in range(k):
        p = [i * c for i, c in enumerate(p)][1:] or [Fr(0)]
    return p


def pint01(p):
    return sum(Fr(c) / (i + 1) for i, c in enumerate(p))


# cubic Hermite basis on [0, 1]; index k = 2 e + d: e = end (0 or 1), d = 1 for the slope
HERM = [[Fr(c) for c in row] for row in ([1, 0, -3, 2], [0, 1, -2, 1], [0, 0, 3, -2], [0, 0, -1, 1])]


def lagrange(p):
    """Lagrange basis of degree p on [0, 1] at the nodes j/p"""
    if p == 1:
        return [[Fr(1), Fr(-1)], [Fr(0), Fr(1)]]
    if p == 2:
        return [[Fr(1), Fr(-3), Fr(2)], [Fr(0), Fr(4), Fr(-4)], [Fr(0), Fr(-1), Fr(2)]]
    raise ValueError("degree must be 1 or 2")


def imat(A, da, B, db, k=0):
    """exact matrix of  int_0^1 xi^k A_i^(da) B_j^(db) dxi"""
    mono = [Fr(0)] * k + [Fr(1)]
    return [[pint01(pmul(mono, pmul(pder(a, da), pder(b_, db)))) for b_ in B] for a in A]


IH = [[imat(HERM, a, HERM, c) for c in range(3)] for a in range(3)]          # IH[a][c][k][l] = int H_k^(a) H_l^(c)
JH = [[[pint01(pmul([Fr(0)] * k + [Fr(1)], pder(Hk, a))) for Hk in HERM] for a in range(3)] for k in range(2)]  # JH[k][a][i] = int xi^k H_i^(a)

# strains: e_m contains d_i u_a for (m, a, i) in
STRAIN = [(0, 0, 0), (1, 1, 1), (2, 0, 1), (2, 1, 0)]
# stress of an Airy function: s_m = sign * d1^alpha d2^beta phi for (m, alpha, beta, sign) in
AIRY = [(0, 0, 2, 1), (1, 2, 0, 1), (2, 1, 1, -1)]
# local unknowns of a cell that carry phi, d1 phi, d2 phi at its lower corner (they fix the affine part)
AFF = {(0, 0), (1, 0), (0, 1)}


def inv3(D):
    """exact inverse of a 3x3 matrix of Fractions"""
    a = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(3)] for i, row in enumerate(D)]
    for c in range(3):
        piv = next(r for r in range(c, 3) if a[r][c] != 0)
        a[c], a[piv] = a[piv], a[c]
        a[c] = [x / a[c][c] for x in a[c]]
        for r in range(3):
            if r != c and a[r][c] != 0:
                a[r] = [x - a[r][c] * y for x, y in zip(a[r], a[c])]
    return [row[3:] for row in a]


def nullspace_and_particular(rows, rhs, ncol):
    """exact solution set of  rows . c = rhs  (lists of Fractions): (particular solution or None, nullspace basis)"""
    m = [list(r) + [v] for r, v in zip(rows, rhs)]
    piv = []
    r = 0
    for c in range(ncol):
        k = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if k is None:
            continue
        m[r], m[k] = m[k], m[r]
        m[r] = [x / m[r][c] for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                m[i] = [x - m[i][c] * y for x, y in zip(m[i], m[r])]
        piv.append(c)
        r += 1
    if any(row[-1] != 0 for row in m[r:]):
        return None, []
    part = [Fr(0)] * ncol
    for i, c in enumerate(piv):
        part[c] = m[i][-1]
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for f in free:
        v = [Fr(0)] * ncol
        v[f] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -m[i][f]
        basis.append(v)
    return part, basis


# ================================================================================================================
# the problem
# ================================================================================================================
class Problem:
    """box: [(a1, b1), (a2, b2)] rationals;  N: cells per direction;  D: 3x3 stiffness (rationals, symmetric positive
    definite) or lam, mu (isotropic, plane strain);  b: constant body force (rationals);
    bc[i] = (type at x_i = a_i, type at x_i = b_i), 'D' (displacement prescribed) or 'N' (no traction);
    degree: 1 or 2 (the displacement space Q_p);
    g: function of the node coordinate arrays returning (g1, g2), the Dirichlet nodal values (each float is taken as
    the exact binary number it is).  THE DIRICHLET DATUM OF THE CERTIFIED PROBLEM IS THE TRACE OF THE Q_p INTERPOLANT
    OF THESE NODAL VALUES on the boundary grid."""

    def __init__(self, box, N, D=None, lam=None, mu=None, b=(0, 0), bc=(('D', 'D'), ('D', 'D')), degree=1, g=None,
                 gval=None, name=''):
        self.name = name
        self.box = [(Fr(a), Fr(c)) for a, c in box]
        self.N = tuple(int(x) for x in N)
        assert len(self.box) == 2 == len(self.N) and all(c > a for a, c in self.box) and min(self.N) >= 1
        if D is None:
            lam, mu = Fr(lam), Fr(mu)
            D = [[lam + 2 * mu, lam, 0], [lam, lam + 2 * mu, 0], [0, 0, mu]]
        self.D = [[Fr(x) for x in row] for row in D]
        d = self.D
        assert all(d[i][j] == d[j][i] for i in range(3) for j in range(3)), "D must be symmetric"
        m2 = d[0][0] * d[1][1] - d[0][1] ** 2
        det = (d[0][0] * (d[1][1] * d[2][2] - d[1][2] ** 2) - d[0][1] * (d[0][1] * d[2][2] - d[1][2] * d[0][2])
               + d[0][2] * (d[0][1] * d[1][2] - d[1][1] * d[0][2]))
        assert d[0][0] > 0 and m2 > 0 and det > 0, "D must be positive definite"
        self.S = inv3(self.D)
        self.b = [Fr(x) for x in b]
        self.bc = [tuple(x) for x in bc]
        assert any(t == 'D' for pair in self.bc for t in pair), "at least one Dirichlet face is needed"
        self.p = int(degree)
        self.h = [(c - a) / Ni for (a, c), Ni in zip(self.box, self.N)]
        self.M = tuple(self.p * Ni + 1 for Ni in self.N)                    # displacement nodes per direction
        self.unodes = [[a + j * h / self.p for j in range(Mi)] for (a, c), h, Mi in zip(self.box, self.h, self.M)]
        self.vnodes = [[a + j * h for j in range(Ni + 1)] for (a, c), h, Ni in zip(self.box, self.h, self.N)]
        fixed = np.zeros(self.M, dtype=bool)
        natural = np.zeros([Ni + 1 for Ni in self.N], dtype=bool)          # nodes of the Airy grid on Gamma_N
        for i in range(2):
            for side in (0, 1):
                sl = [slice(None)] * 2
                sl[i] = 0 if side == 0 else -1
                if self.bc[i][side] == 'D':
                    fixed[tuple(sl)] = True
                else:
                    natural[tuple(sl)] = True
        self.fixed, self.natural = fixed, natural
        if gval is not None:
            self.gval = np.array(gval, dtype=float)
            assert self.gval.shape == (2,) + self.M
        else:
            self.gval = np.zeros((2,) + self.M)
            if g is not None:
                X = np.meshgrid(*[np.array([float(x) for x in ax]) for ax in self.unodes], indexing='ij')
                gv = g(*X)
                for a in range(2):
                    self.gval[a][fixed] = (np.asarray(gv[a], dtype=float) * np.ones(self.M))[fixed]
        self.gval[:, ~fixed] = 0.0
        assert np.all(np.isfinite(self.gval))

    def spec(self):
        return dict(name=self.name, box=[[str(a), str(c)] for a, c in self.box], N=list(self.N),
                    D=[[str(x) for x in row] for row in self.D], b=[str(x) for x in self.b],
                    bc=[list(x) for x in self.bc], degree=self.p)

    @staticmethod
    def from_spec(s, gval):
        return Problem(box=[(Fr(a), Fr(c)) for a, c in s['box']], N=s['N'], D=[[Fr(x) for x in row] for row in s['D']],
                       b=[Fr(x) for x in s['b']], bc=s['bc'], degree=s['degree'], name=s.get('name', ''), gval=gval)

    # ---- exact local matrices -----------------------------------------------------------------------------------
    def stiffness_local(self):
        """exact element matrix K[(a, P), (c, Q)] of  int_cell e^T D e  and load vector f[(a, P)] = int_cell b_a phi_P,
        P = (i1, i2) the local node"""
        if not hasattr(self, '_Kloc'):
            L = lagrange(self.p)
            IL = [[imat(L, a, L, c) for c in range(2)] for a in range(2)]
            h = self.h
            loc = list(itertools.product(range(self.p + 1), repeat=2))
            K = {}
            for (m, a, i) in STRAIN:
                for (n, c, j) in STRAIN:
                    if self.D[m][n] == 0:
                        continue
                    for P in loc:
                        for Q in loc:
                            v = h[0] * h[1] / (h[i] * h[j])
                            for k in range(2):
                                v *= IL[int(k == i)][int(k == j)][P[k]][Q[k]]
                            if v != 0:
                                key = ((a, P), (c, Q))
                                K[key] = K.get(key, Fr(0)) + self.D[m][n] * v
            intL = [pint01(q) for q in L]
            f = {(a, P): self.b[a] * h[0] * h[1] * intL[P[0]] * intL[P[1]] for a in range(2) for P in loc}
            # a constant displacement has no strain: K annihilates it (upper_bound relies on this)
            for a in range(2):
                for row in [(c, Q) for c in range(2) for Q in loc]:
                    assert sum(K.get((row, (a, Pn)), Fr(0)) for Pn in loc) == 0
            self._Kloc = ({k: v for k, v in K.items() if v != 0}, f, loc)
        return self._Kloc

    def airy_local(self):
        """exact element matrix A[k, l] of  int_cell s^T S s  for the stress of the 16 bicubic basis functions
        (k = (k1, k2), k_i = 2 e_i + d_i), and the vectors w0, w1, w2 with
        int_cell s^T S s_b = w0 + x_c1 w1 + x_c2 w2  (x_c the lower corner of the cell, s_b = (-b1 x1, -b2 x2, 0))"""
        if not hasattr(self, '_Aloc'):
            h, S, b = self.h, self.S, self.b
            loc = list(itertools.product(range(4), repeat=2))
            scale = {k: h[0] ** (k[0] % 2) * h[1] ** (k[1] % 2) for k in loc}
            A = {}
            for (m, am, bm, sm) in AIRY:
                for (n, an, bn, sn) in AIRY:
                    if S[m][n] == 0:
                        continue
                    c = S[m][n] * sm * sn * h[0] * h[1] / (h[0] ** (am + an) * h[1] ** (bm + bn))
                    for k in loc:
                        for l in loc:
                            v = c * scale[k] * scale[l] * IH[am][an][k[0]][l[0]] * IH[bm][bn][k[1]][l[1]]
                            if v != 0:
                                A[(k, l)] = A.get((k, l), Fr(0)) + v
            w0, w1, w2 = {}, {}, {}
            for k in loc:
                a0 = a1 = a2 = Fr(0)
                for (m, am, bm, sm) in AIRY:
                    c = scale[k] * sm * h[0] * h[1] / (h[0] ** am * h[1] ** bm)
                    j0 = JH[0][am][k[0]] * JH[0][bm][k[1]]
                    a1 += c * (-b[0] * S[m][0]) * j0
                    a2 += c * (-b[1] * S[m][1]) * j0
                    a0 += c * ((-b[0] * S[m][0]) * h[0] * JH[1][am][k[0]] * JH[0][bm][k[1]]
                               + (-b[1] * S[m][1]) * h[1] * JH[0][am][k[0]] * JH[1][bm][k[1]])
                w0[k], w1[k], w2[k] = a0, a1, a2
            # the stress of an affine Airy function is zero: A and w annihilate the nodal values of 1, x1 - xc1 and
            # x2 - xc2 (lower_bound relies on this)
            for vec in ({k: Fr(int(k[0] % 2 == 0 and k[1] % 2 == 0)) for k in loc},
                        {k: (h[0] * (k[0] // 2) if (k[0] % 2, k[1] % 2) == (0, 0) else Fr(int((k[0] % 2, k[1] % 2) == (1, 0)))) for k in loc},
                        {k: (h[1] * (k[1] // 2) if (k[0] % 2, k[1] % 2) == (0, 0) else Fr(int((k[0] % 2, k[1] % 2) == (0, 1)))) for k in loc}):
                assert all(sum(A.get((k, l), Fr(0)) * vec[l] for l in loc) == 0 for k in loc)
                assert all(sum(w[k] * vec[k] for k in loc) == 0 for w in (w0, w1, w2))
            self._Aloc = ({k: v for k, v in A.items() if v != 0}, w0, w1, w2, loc)
        return self._Aloc

    def face_local(self):
        """exact matrices for the boundary term: HL2[k][q] = int H_k'' L_q, HL1[k][q] = int H_k' L_q, LI[q] = int L_q"""
        if not hasattr(self, '_Floc'):
            L = lagrange(self.p)
            HL2, HL1 = imat(HERM, 2, L, 0), imat(HERM, 1, L, 0)
            # the second derivative of an affine function and the first derivative of a constant vanish
            for q in range(self.p + 1):
                assert HL2[0][q] + HL2[2][q] == 0 and HL2[1][q] + HL2[2][q] + HL2[3][q] == 0 and HL1[0][q] + HL1[2][q] == 0
            self._Floc = (HL2, HL1, [pint01(q) for q in L])
        return self._Floc

    # ---- the traction-free faces ----------------------------------------------------------------------------------
    def traces(self):
        """exact nodal values on Gamma_N of the polynomials psi_b, psi_1, ..., psi_r of docs/ELASTICITY.md, Section 3:
        returns (nodes, Tb, T) with nodes the list of Airy nodes (i, j) on Gamma_N, Tb[node][t] and T[k][node][t]
        Fractions, t = 0..3 for (psi, d1 psi, d2 psi, d12 psi)"""
        if hasattr(self, '_traces'):
            return self._traces
        nodes = [tuple(ix) for ix in np.argwhere(self.natural)]
        if not nodes:
            self._traces = ([], [], [])
            return self._traces
        # unknown polynomial psi = sum c[m][n] x1^m x2^n, m, n = 0..3 (index 4 m + n)
        rows, rhs = [], []
        for i in range(2):
            for side in (0, 1):
                if self.bc[i][side] != 'N':
                    continue
                a = self.box[i][side]
                for r in range(2):                      # the second tangential derivative: coefficient of x_t^r
                    row = [Fr(0)] * 16
                    for q in range(4):                  # q: power of the normal variable
                        idx = 4 * q + (r + 2) if i == 0 else 4 * (r + 2) + q
                        row[idx] = a ** q * (r + 2) * (r + 1)
                    rows.append(row)
                    rhs.append(self.b[i] * a if r == 0 else Fr(0))
                for r in range(3):                      # the mixed derivative: coefficient of x_t^r
                    row = [Fr(0)] * 16
                    for q in range(1, 4):
                        idx = 4 * q + (r + 1) if i == 0 else 4 * (r + 1) + q
                        row[idx] = q * a ** (q - 1) * (r + 1)
                    rows.append(row)
                    rhs.append(Fr(0))
        part, basis = nullspace_and_particular(rows, rhs, 16)
        assert part is not None, "no polynomial particular solution on the traction-free faces"

        def pattern(c):
            out = []
            for (i, j) in nodes:
                x, y = self.vnodes[0][i], self.vnodes[1][j]
                v = [Fr(0)] * 4
                for m in range(4):
                    for n in range(4):
                        cm = c[4 * m + n]
                        if cm == 0:
                            continue
                        v[0] += cm * x ** m * y ** n
                        if m:
                            v[1] += cm * m * x ** (m - 1) * y ** n
                        if n:
                            v[2] += cm * x ** m * n * y ** (n - 1)
                        if m and n:
                            v[3] += cm * m * x ** (m - 1) * n * y ** (n - 1)
                out.append(v)
            return out

        def flat(pt):
            return [x for v in pt for x in v]

        # keep the polynomials whose traces are independent of each other and of the traces of the affine functions
        aff = []
        for idx in (0, 4, 1):                           # 1, x1, x2
            c = [Fr(0)] * 16
            c[idx] = Fr(1)
            aff.append(c)
        reduced = []                                    # rows in echelon form: (pivot column, row)

        def residual(v):
            v = list(v)
            for pc, row in reduced:
                if v[pc] != 0:
                    f = v[pc]
                    v = [x - f * y for x, y in zip(v, row)]
            return v

        def push(v):
            v = residual(v)
            pc = next((i for i, x in enumerate(v) if x != 0), None)
            if pc is None:
                return False
            reduced.append((pc, [x / v[pc] for x in v]))
            return True

        for c in aff:
            push(flat(pattern(c)))
        T = []
        for c in basis:
            pt = pattern(c)
            if push(flat(pt)):
                T.append(pt)
        self._traces = (nodes, pattern(part), T)
        return self._traces


# ================================================================================================================
# sums
# ================================================================================================================
def _fsum_up(a):
    return float(up(up(math.fsum(np.asarray(a, dtype=np.float64).ravel().tolist()))))


def _fsum_dn(a):
    return float(dn(dn(math.fsum(np.asarray(a, dtype=np.float64).ravel().tolist()))))


def _total(x):
    """enclosure of the sum of an interval array"""
    return I(_fsum_dn(x.lo * np.ones(x.shape)), _fsum_up(x.hi * np.ones(x.shape)))


def _quadratic(K, X):
    """enclosure, cell by cell, of  sum_{k, l} K[k, l] X[k] X[l]  (K exact and symmetric, X a dict of interval arrays)"""
    acc = None
    for (k, l), v in K.items():
        if k == l:
            term = X[k].sqr() * fr(v)
        elif (l, k) in K and repr(k) < repr(l):
            term = X[k] * X[l] * fr(2 * v)
        elif (l, k) in K:
            continue
        else:
            term = X[k] * X[l] * fr(v)
        acc = term if acc is None else acc + term
    return acc


# ================================================================================================================
# the displacement: candidate (heuristic) and its energy (RIGOROUS)
# ================================================================================================================
def _ucells(P, U):
    """for every cell, the nodal values (a, (i1, i2)) as arrays over the cells"""
    p = P.p
    return {(a, (i, j)): U[a, i:i + p * (P.N[0] - 1) + 1:p, j:j + p * (P.N[1] - 1) + 1:p]
            for a in range(2) for i in range(p + 1) for j in range(p + 1)}


def upper_bound(P, U):
    """RIGOROUS enclosure (lo, hi) of E[u_h], u_h the Q_p function with nodal values U (shape (2,) + P.M).
    E[u_h] = sum over the cells of  1/2 u^T K u - f.u  with the exact element matrix K and load vector f.
    K annihilates constants, so u^T K u = d^T K d with d the differences from the nodal value at the lower corner of
    the cell, component by component; f.u = f.d + (sum of f over the nodes of the component) u(lower corner)."""
    assert U.shape == (2,) + P.M and np.all(np.isfinite(U))
    assert np.array_equal(U[:, P.fixed], P.gval[:, P.fixed]), "Dirichlet nodal values are not imposed exactly"
    K, f, loc = P.stiffness_local()
    X = {k: I(v) for k, v in _ucells(P, U).items()}
    D = {k: X[k] - X[(k[0], (0, 0))] for k in X if k[1] != (0, 0)}
    en = _quadratic({kl: v for kl, v in K.items() if kl[0][1] != (0, 0) and kl[1][1] != (0, 0)}, D) * fr(Fr(1, 2))
    for k, v in f.items():
        if v != 0 and k[1] != (0, 0):
            en = en - D[k] * fr(v)
    for a in range(2):
        fa = sum(v for k, v in f.items() if k[0] == a)
        if fa != 0:
            en = en - X[(a, (0, 0))] * fr(fa)
    tot = _total(en)
    return float(tot.lo), float(tot.hi)


def _stiffness_global(P):
    K, f, loc = P.stiffness_local()
    M1, M2 = P.M
    p = P.p
    ci, cj = np.meshgrid(np.arange(P.N[0]), np.arange(P.N[1]), indexing='ij')
    ci, cj = ci.ravel(), cj.ravel()

    def gidx(a, Pn):
        return a * M1 * M2 + (p * ci + Pn[0]) * M2 + (p * cj + Pn[1])

    rows, cols, vals = [], [], []
    for (k, l), v in K.items():
        rows.append(gidx(*k))
        cols.append(gidx(*l))
        vals.append(np.full(ci.size, float(v)))
    Kg = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(2 * M1 * M2,) * 2).tocsr()
    F = np.zeros(2 * M1 * M2)
    for k, v in f.items():
        np.add.at(F, gidx(*k), float(v))
    return Kg, F


def solve_candidate(P):
    """heuristic: the finite-element solution in Q_p (any nodal array is a valid candidate)"""
    Kg, F = _stiffness_global(P)
    U = P.gval.copy().ravel()
    fixed = np.concatenate([P.fixed.ravel(), P.fixed.ravel()])
    free = ~fixed
    rhs = F - Kg @ U
    U[free] = spla.spsolve(Kg[free][:, free].tocsc(), rhs[free])
    return U.reshape((2,) + P.M)


# ================================================================================================================
# the stress field
# ================================================================================================================
class Stress:
    """X: the nodal values (phi, d1 phi, d2 phi, d12 phi) of the Airy function at the nodes that are not on Gamma_N
    (shape (4, N1 + 1, N2 + 1); the entries at the nodes of Gamma_N are ignored);  par: the parameters P_1 ... P_r."""

    def __init__(self, P, X=None, par=None):
        nodes, Tb, T = P.traces()
        self.X = np.zeros((4, P.N[0] + 1, P.N[1] + 1)) if X is None else np.array(X, dtype=float)
        self.par = np.zeros(len(T)) if par is None else np.array(par, dtype=float)
        assert self.X.shape == (4, P.N[0] + 1, P.N[1] + 1) and self.par.shape == (len(T),)

    def copy(self, P):
        return Stress(P, self.X.copy(), self.par.copy())


def _pinned(P):
    """if there is no traction-free face the affine part of the Airy function is fixed: phi, d1 phi, d2 phi = 0 at one node"""
    return [] if P.natural.any() else [(0, 0, 0), (1, 0, 0), (2, 0, 0)]


def nodal_enclosure(P, st):
    """RIGOROUS: the nodal values of the Airy function as an interval array.  Nodes not on Gamma_N: the stored floats
    (exact).  Nodes on Gamma_N: enclosures of the exact values of psi_b + sum_k P_k psi_k."""
    nodes, Tb, T = P.traces()
    lo, hi = st.X.copy(), st.X.copy()
    for t, i, j in _pinned(P):
        lo[t, i, j] = hi[t, i, j] = 0.0
    for q, (i, j) in enumerate(nodes):
        for t in range(4):
            v = fr(Tb[q][t])
            for k in range(len(T)):
                if T[k][q][t] != 0:
                    v = v + I(st.par[k]) * fr(T[k][q][t])
            lo[t, i, j], hi[t, i, j] = float(v.lo), float(v.hi)
    return I(lo, hi)


def _xcells(P, X):
    """for every cell, the 16 nodal values k = (k1, k2), k_i = 2 e_i + d_i, as (interval) arrays over the cells"""
    N1, N2 = P.N
    return {(k1, k2): X[(k1 % 2) + 2 * (k2 % 2), (k1 // 2):(k1 // 2) + N1, (k2 // 2):(k2 // 2) + N2]
            for k1 in range(4) for k2 in range(4)}


def _dirichlet_faces(P):
    """(i, side, sign) for every Dirichlet face: normal direction i, outward sign"""
    return [(i, side, -1 if side == 0 else 1) for i in range(2) for side in (0, 1) if P.bc[i][side] == 'D']


def _face_data(P, X, i, side):
    """on the Dirichlet face (i, side): for every face cell the four coefficients cA of the cubic that phi is along the
    face, the four coefficients cB of the cubic that the normal derivative of phi is, and the nodal values of the
    normal and tangential components of g"""
    j = 1 - i
    Nt, ht, p = P.N[j], P.h[j], P.p
    pos = 0 if side == 0 else P.N[i]
    tt, tn = (2, 1) if i == 0 else (1, 2)               # type of the tangential / normal first derivative

    def line(t):
        return X[t, pos, :] if i == 0 else X[t, :, pos]

    cA = [line(0)[0:Nt], line(tt)[0:Nt] * fr(ht), line(0)[1:Nt + 1], line(tt)[1:Nt + 1] * fr(ht)]
    cB = [line(tn)[0:Nt], line(3)[0:Nt] * fr(ht), line(tn)[1:Nt + 1], line(3)[1:Nt + 1] * fr(ht)]
    gpos = 0 if side == 0 else P.M[i] - 1
    G = [P.gval[a, gpos, :] if i == 0 else P.gval[a, :, gpos] for a in range(2)]
    gn = [I(G[i][q:q + p * (Nt - 1) + 1:p]) for q in range(p + 1)]
    gt = [I(G[j][q:q + p * (Nt - 1) + 1:p]) for q in range(p + 1)]
    return cA, cB, gn, gt


def lower_bound(P, st, return_parts=False):
    """RIGOROUS lower bound  int_{Gamma_D} (sigma nu).g - 1/2 int s^T S s  for the stress field st.
    All coefficients are exact rationals; the nodal values are exact floats or enclosures (nodal_enclosure)."""
    A, w0, w1, w2, loc = P.airy_local()
    HL2, HL1, LI = P.face_local()
    X = nodal_enclosure(P, st)
    Xc = _xcells(P, X)
    # differences from the affine function that has the values (phi, d1 phi, d2 phi) of the lower corner of the cell;
    # A and w annihilate affine functions, so X^T A X = d^T A d and w.X = w.d
    a0, a1, a2 = Xc[(0, 0)], Xc[(1, 0)], Xc[(0, 1)]
    Dc = {}
    for k in loc:
        if k in AFF:
            continue
        kind = (k[0] % 2, k[1] % 2)
        if kind == (0, 0):
            Dc[k] = Xc[k] - (a0 + a1 * fr((k[0] // 2) * P.h[0]) + a2 * fr((k[1] // 2) * P.h[1]))
        elif kind == (1, 0):
            Dc[k] = Xc[k] - a1
        elif kind == (0, 1):
            Dc[k] = Xc[k] - a2
        else:
            Dc[k] = Xc[k]
    # the cell integrals:  int s^T S s = Q_phi + 2 Q_cross + Q_b
    q = _quadratic({kl: v for kl, v in A.items() if kl[0] not in AFF and kl[1] not in AFF}, Dc)
    if any(x != 0 for x in P.b):
        # the lower corners of the cells are rationals; enclose them
        x1 = I(np.array([[fr(P.vnodes[0][i]).lo] for i in range(P.N[0])]) * np.ones(P.N),
               np.array([[fr(P.vnodes[0][i]).hi] for i in range(P.N[0])]) * np.ones(P.N))
        x2 = I(np.array([[fr(P.vnodes[1][j]).lo for j in range(P.N[1])]]) * np.ones(P.N),
               np.array([[fr(P.vnodes[1][j]).hi for j in range(P.N[1])]]) * np.ones(P.N))
        for k in loc:
            if k not in AFF and (w0[k] != 0 or w1[k] != 0 or w2[k] != 0):
                q = q + Dc[k] * (fr(w0[k]) + x1 * fr(w1[k]) + x2 * fr(w2[k])) * fr(2)
    (a1, c1), (a2, c2) = P.box
    S, b = P.S, P.b
    qb = (S[0][0] * b[0] ** 2 * (c1 ** 3 - a1 ** 3) / 3 * (c2 - a2)
          + 2 * S[0][1] * b[0] * b[1] * (c1 ** 2 - a1 ** 2) / 2 * (c2 ** 2 - a2 ** 2) / 2
          + S[1][1] * b[1] ** 2 * (c1 - a1) * (c2 ** 3 - a2 ** 3) / 3)
    intQ = _total(q) + fr(qb)
    # the boundary term
    bnd = I(0.0)
    for i, side, sign in _dirichlet_faces(P):
        j = 1 - i
        ht = P.h[j]
        a = P.box[i][side]
        cA, cB, gn, gt = _face_data(P, X, i, side)
        # differences again: the second derivative annihilates the affine function through the value and the slope at
        # the lower end of the cell face, the first derivative the constant through the value there
        dA = {2: cA[2] - (cA[0] + cA[1]), 3: cA[3] - cA[1]}
        dB = {1: cB[1], 2: cB[2] - cB[0], 3: cB[3]}
        acc = None
        for qn in range(P.p + 1):
            for k in range(4):
                if k in dA and HL2[k][qn] != 0:
                    term = dA[k] * gn[qn] * fr(HL2[k][qn] / ht)
                    acc = term if acc is None else acc + term
                if k in dB and HL1[k][qn] != 0:
                    term = dB[k] * gt[qn] * fr(-HL1[k][qn])
                    acc = term if acc is None else acc + term
            if P.b[i] != 0 and a != 0:
                term = gn[qn] * fr(-P.b[i] * a * ht * LI[qn])
                acc = term if acc is None else acc + term
        if acc is not None:
            bnd = bnd + _total(acc) * fr(sign)
    lb = bnd - intQ * fr(Fr(1, 2))
    if return_parts:
        return float(lb.lo), dict(boundary=(float(bnd.lo), float(bnd.hi)), intQ=(float(intQ.lo), float(intQ.hi)))
    return float(lb.lo)


def solve_stress(P):
    """heuristic: the maximiser of the lower bound over the free nodal values and the parameters (one sparse solve)"""
    A, w0, w1, w2, loc = P.airy_local()
    HL2, HL1, LI = P.face_local()
    nodes, Tb, T = P.traces()
    N1, N2 = P.N
    n1, n2 = N1 + 1, N2 + 1
    nd = 4 * n1 * n2

    def gid(t, i, j):
        return (t * n1 + i) * n2 + j

    ci, cj = np.meshgrid(np.arange(N1), np.arange(N2), indexing='ij')
    ci, cj = ci.ravel(), cj.ravel()

    def cid(k):
        return gid((k[0] % 2) + 2 * (k[1] % 2), ci + k[0] // 2, cj + k[1] // 2)

    rows, cols, vals = [], [], []
    for (k, l), v in A.items():
        rows.append(cid(k))
        cols.append(cid(l))
        vals.append(np.full(ci.size, float(v)))
    Ag = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nd, nd)).tocsr()
    wg = np.zeros(nd)
    if any(x != 0 for x in P.b):
        xc1 = np.array([float(P.vnodes[0][i]) for i in ci])
        xc2 = np.array([float(P.vnodes[1][j]) for j in cj])
        for k in loc:
            np.add.at(wg, cid(k), float(w0[k]) + xc1 * float(w1[k]) + xc2 * float(w2[k]))
    beta = np.zeros(nd)                                  # the boundary term is beta . x + const
    for i, side, sign in _dirichlet_faces(P):
        j = 1 - i
        Nt, ht, p = P.N[j], float(P.h[j]), P.p
        pos = 0 if side == 0 else P.N[i]
        tt, tn = (2, 1) if i == 0 else (1, 2)
        gpos = 0 if side == 0 else P.M[i] - 1
        G = [P.gval[a, gpos, :] if i == 0 else P.gval[a, :, gpos] for a in range(2)]
        cells = np.arange(Nt)
        for qn in range(p + 1):
            gn = G[i][qn:qn + p * (Nt - 1) + 1:p]
            gt = G[j][qn:qn + p * (Nt - 1) + 1:p]
            for k in range(4):
                e, d = k // 2, k % 2
                tA, tB = (tt if d else 0), (3 if d else tn)
                sc = ht if d else 1.0
                ia_ = gid(tA, pos, cells + e) if i == 0 else gid(tA, cells + e, pos)
                ib_ = gid(tB, pos, cells + e) if i == 0 else gid(tB, cells + e, pos)
                np.add.at(beta, ia_, sign * sc * float(HL2[k][qn]) / ht * gn)
                np.add.at(beta, ib_, -sign * sc * float(HL1[k][qn]) * gt)
    # x = Tm z + x0: free nodal values, then the parameters
    x0 = np.zeros(nd)
    isnat = np.zeros(nd, dtype=bool)
    for q, (i, j) in enumerate(nodes):
        for t in range(4):
            isnat[gid(t, i, j)] = True
            x0[gid(t, i, j)] = float(Tb[q][t])
    for t, i, j in _pinned(P):
        isnat[gid(t, i, j)] = True
    free = np.flatnonzero(~isnat)
    nz = free.size + len(T)
    Tm = sp.lil_matrix((nd, nz))
    Tm[free, np.arange(free.size)] = 1.0
    for k in range(len(T)):
        for q, (i, j) in enumerate(nodes):
            for t in range(4):
                if T[k][q][t] != 0:
                    Tm[gid(t, i, j), free.size + k] = float(T[k][q][t])
    Tm = Tm.tocsr()
    H = (Tm.T @ Ag @ Tm).tocsc()
    rhs = Tm.T @ (beta - wg - Ag @ x0)
    z = spla.spsolve(H, rhs)
    X = np.zeros(nd)
    X[free] = z[:free.size]
    return Stress(P, X.reshape(4, n1, n2), z[free.size:])


# ================================================================================================================
# drivers and certificates on disk
# ================================================================================================================
def certify(P, U=None, st=None):
    """both bounds: candidate (finite-element solution unless given) and stress field (dual solve unless given)"""
    import time
    t0 = time.time()
    U = solve_candidate(P) if U is None else U
    st = solve_stress(P) if st is None else st
    ub = upper_bound(P, U)
    lb = lower_bound(P, st)
    return dict(N=list(P.N), degree=P.p, lower=lb, upper=ub[1], upper_enclosure=list(ub),
                rel_gap=(ub[1] - lb) / abs(ub[1]) if ub[1] != 0 else float('nan'), seconds=time.time() - t0), U, st


def claim_strings(lower, upper, digits=10):
    """the two claimed numbers as exact decimal strings: the computed bounds weakened by one part in 10^9 of the larger
    modulus and rounded outward to `digits` significant digits"""
    lo, hi = Fr(float(lower)), Fr(float(upper))
    scale = max(abs(lo), abs(hi))
    if scale == 0:
        return '0', '0'
    ex = math.floor(math.log10(float(scale))) - digits + 1
    q = Fr(10) ** ex
    mlo = math.floor((lo - scale / 10 ** 9) / q)
    mhi = math.ceil((hi + scale / 10 ** 9) / q)

    def dec(m):
        if ex >= 0:
            return str(m * 10 ** ex)
        sgn, a = ('-' if m < 0 else ''), str(abs(m)).rjust(-ex + 1, '0')
        return f'{sgn}{a[:ex]}.{a[ex:]}'

    return dec(mlo), dec(mhi)


def save_certificate(stem, P, U, st, lower, upper, extra=None):
    lo, hi = claim_strings(lower, upper)
    meta = dict(format='CertiField-elasticity-certificate-1', problem=P.spec(), lower=lo, upper=hi, reference=__version__)
    if extra:
        meta['info'] = extra
    json.dump(meta, open(stem + '.json', 'w'), indent=1)
    np.savez_compressed(stem + '.npz', gval=P.gval, U=U, X=st.X, par=st.par)


def load_certificate(stem):
    meta = json.load(open(stem + '.json'))
    z = np.load(stem + '.npz')
    P = Problem.from_spec(meta['problem'], z['gval'])
    return meta, P, np.array(z['U'], dtype=float), Stress(P, z['X'], z['par'])


def verify_certificate(stem, verbose=True):
    """recompute both bounds from the stored exact data; True if they are at least as strong as the stored claims"""
    meta, P, U, st = load_certificate(stem)
    lb = lower_bound(P, st)
    ub = upper_bound(P, U)[1]
    ok = bool(math.isfinite(lb) and math.isfinite(ub) and Fr(lb) >= Fr(meta['lower']) and Fr(ub) <= Fr(meta['upper'])
              and Fr(meta['lower']) <= Fr(meta['upper']))
    if verbose:
        print(f"{stem}: recomputed {lb!r} <= inf E <= {ub!r};  claimed {meta['lower']} <= inf E <= {meta['upper']};  {'OK' if ok else 'FAILED'}")
    return ok, lb, ub


def export_dir(d, P, U, st, expected=None):
    """raw export for the Wolfram Language model: little-endian float64, C order, plus field.json"""
    import os
    os.makedirs(d, exist_ok=True)
    meta = dict(format='CertiField-elasticity-field-1', problem=P.spec(), NodeDims=list(P.M), AiryDims=[P.N[0] + 1, P.N[1] + 1],
                Parameters=len(st.par))
    if expected is not None:
        meta['expected'] = expected
    json.dump(meta, open(os.path.join(d, 'field.json'), 'w'), indent=1)
    np.asarray(P.gval, dtype='<f8').tofile(os.path.join(d, 'gval.f64'))
    np.asarray(U, dtype='<f8').tofile(os.path.join(d, 'U.f64'))
    np.asarray(st.X, dtype='<f8').tofile(os.path.join(d, 'X.f64'))
    np.asarray(st.par, dtype='<f8').tofile(os.path.join(d, 'par.f64'))
