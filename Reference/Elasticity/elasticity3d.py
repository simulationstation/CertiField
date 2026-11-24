#!/usr/bin/env python3
"""CertiField elasticity in three dimensions -- reference implementation (Python/NumPy).

Certified two-sided bounds for linear elasticity on a box,

    E[u] = int_Omega  1/2 e(u)^T D e(u) - b.u dx ,
    e(u) = (d1 u1, d2 u2, d3 u3, d2 u3 + d3 u2, d1 u3 + d3 u1, d1 u2 + d2 u1) ,
    u = g on the Dirichlet faces, no traction on the others.

docs/ELASTICITY_3D.md has the statements and their proofs.

  upper bound : the energy of the continuous piecewise polynomial (Q_p, p = 1 or 2) displacement with the given nodal
                values, enclosed rigorously.
  lower bound : for ANY statically admissible stress field sigma (polynomial on each cell, single-valued traction,
                div sigma + b = 0 exactly, zero traction on the traction-free faces exactly), every admissible u
                satisfies
                    E[u]  >=  int_{Gamma_D} (sigma nu).g  -  1/2 int_Omega s^T S s ,        S = D^-1 .
                The stress field is the Beltrami field of six stress functions Phi_11 ... Phi_12, each a C^1 piecewise
                tricubic function (eight nodal values per node), plus the particular field -diag(b_i x_i).  At the
                nodes of the traction-free faces the nodal values that enter the traction are those of polynomials
                Psi_b + sum_k P_k Psi_k found in exact rational arithmetic, which makes the traction vanish there
                exactly.

Everything in upper_bound and lower_bound is rigorous: outward-rounded interval arithmetic (ia.py) on exact inputs
(binary floats and rationals).  The linear solvers (solve_candidate, solve_stress) are plain floating point and are
never trusted.
"""
import itertools
import json
import math
import time
from fractions import Fraction as Fr

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ia import I, fr, dn, up

__version__ = "0.1.0-reference-3d"

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
JH = [[[pint01(pmul([Fr(0)] * q + [Fr(1)], pder(Hk, a))) for Hk in HERM] for a in range(3)] for q in range(2)]  # JH[q][a][k] = int xi^q H_k^(a)

# ---- the operators ----------------------------------------------------------------------------------------------
# strains: e_m contains d_i u_a for (m, a, i) in
STRAIN = [(0, 0, 0), (1, 1, 1), (2, 2, 2), (3, 2, 1), (3, 1, 2), (4, 2, 0), (4, 0, 2), (5, 1, 0), (5, 0, 1)]
# position of sigma_ij in the vectors of six components (11, 22, 33, 23, 13, 12)
VOIGT = [[0, 5, 4], [5, 1, 3], [4, 3, 2]]
PAIRS = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
FNAMES = ['Phi11', 'Phi22', 'Phi33', 'Phi23', 'Phi13', 'Phi12']


def _levi(i, j, k):
    return (i - j) * (j - k) * (k - i) // 2


def _beltrami_table():
    """sigma_ij = sum eps_ikl eps_jmn d_k d_m Phi_ln:  for every component m a list of (f, alpha, coefficient),
    meaning  s_m = sum coefficient * d^alpha Phi_f"""
    tab = []
    for (i, j) in PAIRS:
        acc = {}
        for k in range(3):
            for l in range(3):
                e1 = _levi(i, k, l)
                if not e1:
                    continue
                for m in range(3):
                    for n in range(3):
                        e2 = _levi(j, m, n)
                        if not e2:
                            continue
                        al = [0, 0, 0]
                        al[k] += 1
                        al[m] += 1
                        key = (VOIGT[l][n], tuple(al))
                        acc[key] = acc.get(key, 0) + e1 * e2
        tab.append([(f, al, c) for (f, al), c in sorted(acc.items()) if c != 0])
    return tab


BELTRAMI = _beltrami_table()
FUNCTION_SETS = {'beltrami': [0, 1, 2, 3, 4, 5], 'maxwell': [0, 1, 2], 'morera': [3, 4, 5]}


def ff(n, k):
    """falling factorial n (n-1) ... (n-k+1)"""
    r = 1
    for q in range(k):
        r *= (n - q)
    return r


def inv_exact(D):
    """exact inverse of a square matrix of Fractions"""
    n = len(D)
    a = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(D)]
    for c in range(n):
        piv = next(r for r in range(c, n) if a[r][c] != 0)
        a[c], a[piv] = a[piv], a[c]
        a[c] = [x / a[c][c] for x in a[c]]
        for r in range(n):
            if r != c and a[r][c] != 0:
                fct = a[r][c]
                a[r] = [x - fct * y for x, y in zip(a[r], a[c])]
    return [row[n:] for row in a]


def is_spd_exact(D):
    """symmetric elimination without pivoting: all pivots positive"""
    n = len(D)
    a = [[Fr(x) for x in row] for row in D]
    if any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        return False
    for c in range(n):
        if a[c][c] <= 0:
            return False
        for r in range(c + 1, n):
            fct = a[r][c] / a[c][c]
            a[r] = [x - fct * y for x, y in zip(a[r], a[c])]
    return True


# ---- exact sparse linear algebra: rows are dicts {column: Fraction} --------------------------------------------------
class Echelon:
    """the reduced row echelon form of the rows inserted so far (it is unique, so the order of insertion does not
    matter); a right-hand side is carried along"""

    def __init__(self):
        self.rows = {}           # pivot column -> (row, rhs); the row has 1 at the pivot and 0 at the other pivots
        self.consistent = True

    def insert(self, row, rhs=Fr(0)):
        row = {c: v for c, v in row.items() if v != 0}
        for c in [c for c in row if c in self.rows]:
            fct = row.pop(c)
            prow, prhs = self.rows[c]
            for cc, v in prow.items():
                if cc != c:
                    nv = row.get(cc, Fr(0)) - fct * v
                    if nv == 0:
                        row.pop(cc, None)
                    else:
                        row[cc] = nv
            rhs = rhs - fct * prhs
        if not row:
            if rhs != 0:
                self.consistent = False
            return False
        pc = min(row)
        piv = row[pc]
        row = {c: v / piv for c, v in row.items()}
        rhs = rhs / piv
        for c, (prow, prhs) in list(self.rows.items()):
            fct = prow.get(pc)
            if fct is not None:
                for cc, v in row.items():
                    nv = prow.get(cc, Fr(0)) - fct * v
                    if nv == 0:
                        prow.pop(cc, None)
                    else:
                        prow[cc] = nv
                self.rows[c] = (prow, prhs - fct * rhs)
        self.rows[pc] = (row, rhs)
        return True

    def particular(self):
        """the solution with all free unknowns zero"""
        return {c: rhs for c, (row, rhs) in self.rows.items() if rhs != 0}

    def nullspace(self, ncol):
        """one vector per free column, in increasing order of the column"""
        basis = []
        for f in range(ncol):
            if f in self.rows:
                continue
            v = {f: Fr(1)}
            for c, (row, rhs) in self.rows.items():
                x = row.get(f)
                if x is not None:
                    v[c] = -x
            basis.append(v)
        return basis


_PRIMES = (2147483647, 2147483629)


class EchelonMod:
    """row echelon form modulo a prime, for rank decisions on dense integer vectors"""

    def __init__(self, n, p):
        self.n, self.p = n, p
        self.piv = []

    def push(self, v):
        p = self.p
        v = np.array(v, dtype=np.int64) % p
        for pc, row in self.piv:
            if v[pc]:
                v = (v - int(v[pc]) * row) % p
        nz = np.flatnonzero(v)
        if nz.size == 0:
            return False
        pc = int(nz[0])
        self.piv.append((pc, (v * pow(int(v[pc]), p - 2, p)) % p))
        return True


def _residue(q, p):
    q = Fr(q)
    return (q.numerator % p) * pow(q.denominator % p, p - 2, p) % p


# ================================================================================================================
# the problem
# ================================================================================================================
class Problem3D:
    """box: [(a1, b1), (a2, b2), (a3, b3)] rationals;  N: cells per direction;  D: 6x6 stiffness (rationals, symmetric
    positive definite, components in the order 11, 22, 33, 23, 13, 12 with engineering shear strains) or lam, mu
    (isotropic);  b: constant body force (rationals);
    bc[i] = (type at x_i = a_i, type at x_i = b_i), 'D' (displacement prescribed) or 'N' (no traction);
    degree: 1 or 2 (the displacement space Q_p);
    functions: 'beltrami' (six stress functions), 'maxwell' or 'morera' (three), or a list of indices into
    (Phi11, Phi22, Phi33, Phi23, Phi13, Phi12);
    g: function of the node coordinate arrays returning (g1, g2, g3), the Dirichlet nodal values (each float is taken
    as the exact binary number it is).  THE DIRICHLET DATUM OF THE CERTIFIED PROBLEM IS THE TRACE OF THE Q_p
    INTERPOLANT OF THESE NODAL VALUES on the boundary grid."""

    def __init__(self, box, N, D=None, lam=None, mu=None, b=(0, 0, 0), bc=(('D', 'D'), ('D', 'D'), ('D', 'D')),
                 degree=1, g=None, gval=None, functions='beltrami', states='resultants', name=''):
        self.name = name
        assert states in ('resultants', 'all')
        self.state_rule = states
        self.box = [(Fr(a), Fr(c)) for a, c in box]
        self.N = tuple(int(x) for x in N)
        assert len(self.box) == 3 == len(self.N) and all(c > a for a, c in self.box) and min(self.N) >= 1
        if D is None:
            lam, mu = Fr(lam), Fr(mu)
            D = [[(lam + 2 * mu if i == j else lam) if (i < 3 and j < 3) else (mu if i == j else 0) for j in range(6)]
                 for i in range(6)]
        self.D = [[Fr(x) for x in row] for row in D]
        assert len(self.D) == 6 and all(len(r) == 6 for r in self.D)
        assert is_spd_exact(self.D), "D must be symmetric positive definite"
        self.S = inv_exact(self.D)
        self.b = [Fr(x) for x in b]
        self.bc = [tuple(x) for x in bc]
        assert all(t in ('D', 'N') for pair in self.bc for t in pair)
        assert any(t == 'D' for pair in self.bc for t in pair), "at least one Dirichlet face is needed"
        self.p = int(degree)
        self.funcs = list(FUNCTION_SETS[functions]) if isinstance(functions, str) else [int(f) for f in functions]
        self.functions = functions if isinstance(functions, str) else list(self.funcs)
        self.nf = len(self.funcs)
        # the table of the stress operator restricted to the functions in use: (local function index, alpha, coef)
        self.stress = [[(self.funcs.index(f), al, c) for (f, al, c) in BELTRAMI[m] if f in self.funcs] for m in range(6)]
        self.h = [(c - a) / Ni for (a, c), Ni in zip(self.box, self.N)]
        self.M = tuple(self.p * Ni + 1 for Ni in self.N)                    # displacement nodes per direction
        self.n = tuple(Ni + 1 for Ni in self.N)                             # nodes of the grid per direction
        self.unodes = [[a + j * h / self.p for j in range(Mi)] for (a, c), h, Mi in zip(self.box, self.h, self.M)]
        self.vnodes = [[a + j * h for j in range(Ni + 1)] for (a, c), h, Ni in zip(self.box, self.h, self.N)]
        fixed = np.zeros(self.M, dtype=bool)
        for i in range(3):
            for side in (0, 1):
                if self.bc[i][side] == 'D':
                    sl = [slice(None)] * 3
                    sl[i] = 0 if side == 0 else -1
                    fixed[tuple(sl)] = True
        self.fixed = fixed
        if gval is not None:
            self.gval = np.array(gval, dtype=float)
            assert self.gval.shape == (3,) + self.M
        else:
            self.gval = np.zeros((3,) + self.M)
            if g is not None:
                X = np.meshgrid(*[np.array([float(x) for x in ax]) for ax in self.unodes], indexing='ij')
                gv = g(*X)
                for a in range(3):
                    self.gval[a][fixed] = (np.asarray(gv[a], dtype=float) * np.ones(self.M))[fixed]
        self.gval[:, ~fixed] = 0.0
        assert np.all(np.isfinite(self.gval))

    def spec(self):
        return dict(name=self.name, box=[[str(a), str(c)] for a, c in self.box], N=list(self.N),
                    D=[[str(x) for x in row] for row in self.D], b=[str(x) for x in self.b],
                    bc=[list(x) for x in self.bc], degree=self.p, functions=self.functions, states=self.state_rule)

    @staticmethod
    def from_spec(s, gval):
        return Problem3D(box=[(Fr(a), Fr(c)) for a, c in s['box']], N=s['N'], D=[[Fr(x) for x in row] for row in s['D']],
                         b=[Fr(x) for x in s['b']], bc=s['bc'], degree=s['degree'], functions=s.get('functions', 'beltrami'),
                         states=s.get('states', 'resultants'), name=s.get('name', ''), gval=gval)

    def free_faces(self):
        return [(i, side) for i in range(3) for side in (0, 1) if self.bc[i][side] == 'N']

    def dirichlet_faces(self):
        """(i, side, sign): normal direction, side, sign of the outward normal"""
        return [(i, side, -1 if side == 0 else 1) for i in range(3) for side in (0, 1) if self.bc[i][side] == 'D']

    def relevant(self, i):
        """the pairs (function, order m of the normal derivative) that enter the traction on a plane normal to e_i
        (the list (R_i) of docs/ELASTICITY_3D.md)"""
        out = set()
        for j in range(3):
            for (f, al, c) in self.stress[VOIGT[i][j]]:
                assert al[i] <= 1
                out.add((f, al[i]))
        return sorted(out)

    # ---- exact local matrices -----------------------------------------------------------------------------------
    def stiffness_local(self):
        """exact element matrix K[(a, P), (c, Q)] of  int_cell e^T D e  and load vector f[(a, P)] = int_cell b_a phi_P,
        P = (i1, i2, i3) the local node"""
        if not hasattr(self, '_Kloc'):
            L = lagrange(self.p)
            IL = [[imat(L, a, L, c) for c in range(2)] for a in range(2)]
            h = self.h
            vol = h[0] * h[1] * h[2]
            loc = list(itertools.product(range(self.p + 1), repeat=3))
            K = {}
            for (m, a, i) in STRAIN:
                for (n, c, j) in STRAIN:
                    if self.D[m][n] == 0:
                        continue
                    tabs = [IL[int(k == i)][int(k == j)] for k in range(3)]
                    pref = self.D[m][n] * vol / (h[i] * h[j])
                    for P in loc:
                        for Q in loc:
                            v = tabs[0][P[0]][Q[0]]
                            if v == 0:
                                continue
                            v = v * tabs[1][P[1]][Q[1]] * tabs[2][P[2]][Q[2]]
                            if v != 0:
                                key = ((a, P), (c, Q))
                                K[key] = K.get(key, Fr(0)) + pref * v
            intL = [pint01(q) for q in L]
            f = {(a, P): self.b[a] * vol * intL[P[0]] * intL[P[1]] * intL[P[2]] for a in range(3) for P in loc}
            self._Kloc = ({k: v for k, v in K.items() if v != 0}, f, loc)
        return self._Kloc

    def stress_local(self):
        """exact element matrix A[(f, k), (g, l)] of  int_cell s^T S s  for the stress of the nf * 64 tricubic basis
        functions (f the function, k = (k1, k2, k3), k_i = 2 e_i + d_i), and vectors w0, w1, w2, w3 with
        int_cell s^T S s_b = w0 + x_c1 w1 + x_c2 w2 + x_c3 w3  (x_c the lower corner of the cell,
        s_b = (-b1 x1, -b2 x2, -b3 x3, 0, 0, 0))"""
        if not hasattr(self, '_Aloc'):
            h, S, b = self.h, self.S, self.b
            loc = [(f, k) for f in range(self.nf) for k in itertools.product(range(4), repeat=3)]
            # scaled one-dimensional tables: T[i][a][c][k][l] = h_i^(1 + k%2 + l%2 - a - c) int H_k^(a) H_l^(c)
            T = [[[[[h[i] ** (1 + k % 2 + l % 2 - a - c) * IH[a][c][k][l] for l in range(4)] for k in range(4)]
                   for c in range(3)] for a in range(3)] for i in range(3)]
            A = {}
            ks = list(itertools.product(range(4), repeat=3))
            for m in range(6):
                for n in range(6):
                    if S[m][n] == 0:
                        continue
                    for (f, al, c1) in self.stress[m]:
                        for (g, be, c2) in self.stress[n]:
                            cc = S[m][n] * c1 * c2
                            t0, t1, t2 = T[0][al[0]][be[0]], T[1][al[1]][be[1]], T[2][al[2]][be[2]]
                            for k in ks:
                                r0, r1, r2 = t0[k[0]], t1[k[1]], t2[k[2]]
                                for l in ks:
                                    v = r0[l[0]]
                                    if v == 0:
                                        continue
                                    v = v * r1[l[1]] * r2[l[2]]
                                    if v != 0:
                                        key = ((f, k), (g, l))
                                        A[key] = A.get(key, Fr(0)) + cc * v
            # J[i][q][a][k] = h_i^(1 + k%2 - a) int xi^q H_k^(a)
            J = [[[[h[i] ** (1 + k % 2 - a) * JH[q][a][k] for k in range(4)] for a in range(3)] for q in range(2)] for i in range(3)]
            w = [dict() for _ in range(4)]
            for m in range(6):
                for (f, al, c1) in self.stress[m]:
                    for n in range(3):
                        if S[m][n] == 0 or b[n] == 0:
                            continue
                        cc = c1 * S[m][n] * (-b[n])
                        for k in ks:
                            j0 = [J[i][0][al[i]][k[i]] for i in range(3)]
                            v = cc * j0[0] * j0[1] * j0[2]
                            if v != 0:
                                w[n + 1][(f, k)] = w[n + 1].get((f, k), Fr(0)) + v
                            j1 = list(j0)
                            j1[n] = J[n][1][al[n]][k[n]]
                            v = cc * h[n] * j1[0] * j1[1] * j1[2]
                            if v != 0:
                                w[0][(f, k)] = w[0].get((f, k), Fr(0)) + v
            self._Aloc = ({k: v for k, v in A.items() if v != 0}, w, loc)
        return self._Aloc

    def face_local(self):
        """exact tables for the boundary term: HL[a][k][q] = int H_k^(a) L_q (a = 0, 1, 2), LI[q] = int L_q"""
        if not hasattr(self, '_Floc'):
            L = lagrange(self.p)
            self._Floc = ([imat(HERM, a, L, 0) for a in range(3)], [pint01(q) for q in L])
        return self._Floc

    # ---- the traction-free faces ----------------------------------------------------------------------------------
    def _col(self, f, ga):
        return f * 64 + 16 * ga[0] + 4 * ga[1] + ga[2]

    def states(self):
        """the polynomials Psi_b, Psi_1, ..., Psi_r of docs/ELASTICITY_3D.md, Section 3.4, as sparse vectors
        {column: Fraction} over the nf * 64 coefficients (column = 64 f + 16 g1 + 4 g2 + g3 for the coefficient of
        x1^g1 x2^g2 x3^g3 in function f).  Exact, except that the choice of Psi_1 ... Psi_r among the echelon basis
        of the homogeneous solutions is made with ranks computed modulo two primes (the choice affects tightness only;
        every vector returned satisfies the traction conditions exactly)."""
        if hasattr(self, '_states'):
            return self._states
        faces = self.free_faces()
        ncol = self.nf * 64
        if not faces:
            self._states = ({}, [])
            return self._states
        # 1. the traction conditions
        ech = Echelon()
        for (i, side) in faces:
            a = self.box[i][side]
            t1, t2 = [t for t in range(3) if t != i]
            for j in range(3):
                for r1 in range(4):
                    for r2 in range(4):
                        row = {}
                        for (f, al, c) in self.stress[VOIGT[i][j]]:
                            g1, g2 = r1 + al[t1], r2 + al[t2]
                            if g1 > 3 or g2 > 3:
                                continue
                            for gi in range(al[i], 4):
                                ga = [0, 0, 0]
                                ga[i], ga[t1], ga[t2] = gi, g1, g2
                                col = self._col(f, ga)
                                row[col] = row.get(col, Fr(0)) + c * ff(g1, al[t1]) * ff(g2, al[t2]) * ff(gi, al[i]) * a ** (gi - al[i])
                        rhs = self.b[i] * a if (j == i and r1 == 0 and r2 == 0) else Fr(0)
                        ech.insert(row, rhs)
        if not ech.consistent:
            raise ValueError("no polynomial particular field for this body force and these traction-free faces")
        psi_b = ech.particular()
        basis = ech.nullspace(ncol)
        # 2. the null sextuples (stress identically zero)
        null = Echelon()
        for m in range(6):
            for ga in itertools.product(range(4), repeat=3):
                row = {}
                for (f, al, c) in self.stress[m]:
                    tg = [ga[k] + al[k] for k in range(3)]
                    if max(tg) > 3:
                        continue
                    col = self._col(f, tg)
                    row[col] = row.get(col, Fr(0)) + c * ff(tg[0], al[0]) * ff(tg[1], al[1]) * ff(tg[2], al[2])
                null.insert(row)
        nullbasis = null.nullspace(ncol)
        # 3. the face data: for every traction-free face and every pair (f, m) of the list (R_i) the polynomial
        #    d_i^m Psi_f restricted to the face
        cols = {}                                   # column -> list of (row of the face data, value)
        nrow = 0
        for (i, side) in faces:
            a = self.box[i][side]
            t1, t2 = [t for t in range(3) if t != i]
            for (f, m) in self.relevant(i):
                for r1 in range(4):
                    for r2 in range(4):
                        for gi in range(m, 4):
                            ga = [0, 0, 0]
                            ga[i], ga[t1], ga[t2] = gi, r1, r2
                            cols.setdefault(self._col(f, ga), []).append((nrow, ff(gi, m) * a ** (gi - m)))
                        nrow += 1

        def facedata(v, p):
            out = np.zeros(nrow, dtype=np.int64)
            for c, x in v.items():
                xr = _residue(x, p)
                for (r, val) in cols.get(c, ()):
                    out[r] = (out[r] + xr * _residue(val, p)) % p
            return out

        chosen = None
        for p in _PRIMES:
            em = EchelonMod(nrow, p)
            for v in nullbasis:
                em.push(facedata(v, p))
            sel = [q for q, v in enumerate(basis) if em.push(facedata(v, p))]
            assert chosen is None or chosen == sel, "the two modular rank computations disagree"
            chosen = sel
        self._all_states = (psi_b, [basis[q] for q in chosen])
        if self.state_rule == 'all':
            self._states = self._all_states
            return self._states
        # the default rule: only the states that carry a resultant force or moment through a ring of
        # traction-free faces (exact rank decisions on vectors of six rationals)
        axis = self.ring_axis()
        keep = []
        if axis is not None:
            ech6 = Echelon()
            for v in basis:
                res = self.resultants(v, axis)
                if ech6.insert({q: x for q, x in enumerate(res) if x != 0}):
                    keep.append(v)
                    if len(keep) == 6:
                        break
        self._states = (psi_b, keep)
        return self._states

    def ring_axis(self):
        """a if the two faces normal to e_a are Dirichlet faces and the four others are traction-free, else None"""
        for a in range(3):
            if all(self.bc[i] == (('D', 'D') if i == a else ('N', 'N')) for i in range(3)):
                return a
        return None

    def stress_polynomials(self, v):
        """the six stress components of the sextuple with coefficient vector v, as dicts {(g1, g2, g3): Fraction}"""
        out = [dict() for _ in range(6)]
        byf = {}
        for c, x in v.items():
            byf.setdefault(c // 64, []).append((((c % 64) // 16, (c % 16) // 4, c % 4), x))
        for m in range(6):
            for (f, al, c) in self.stress[m]:
                for ga, x in byf.get(f, ()):
                    if any(ga[k] < al[k] for k in range(3)):
                        continue
                    key = tuple(ga[k] - al[k] for k in range(3))
                    val = out[m].get(key, Fr(0)) + c * x * ff(ga[0], al[0]) * ff(ga[1], al[1]) * ff(ga[2], al[2])
                    if val == 0:
                        out[m].pop(key, None)
                    else:
                        out[m][key] = val
        return out

    def resultants(self, v, axis):
        """force (three components) and moment about the origin (three components) of the traction sigma e_axis of
        the sextuple v on the cross-section x_axis = lower end of the box; exact"""
        s = self.stress_polynomials(v)
        c0 = self.box[axis][0]
        t1, t2 = [t for t in range(3) if t != axis]

        def integral(poly, extra):
            """int over the cross-section of poly * x^extra (extra: exponents added to the monomials)"""
            tot = Fr(0)
            for ga, x in poly.items():
                e = [ga[k] + extra[k] for k in range(3)]
                val = x * c0 ** e[axis]
                for t in (t1, t2):
                    lo, hi = self.box[t]
                    val *= (hi ** (e[t] + 1) - lo ** (e[t] + 1)) / (e[t] + 1)
                tot += val
            return tot

        trac = [s[VOIGT[axis][j]] for j in range(3)]
        force = [integral(trac[j], (0, 0, 0)) for j in range(3)]
        unit = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        moment = [sum(_levi(i, j, k) * integral(trac[k], unit[j]) for j in range(3) for k in range(3)) for i in range(3)]
        return force + moment

    def prescribed(self):
        """the nodal values fixed by the traction-free faces:
        mask[f, t, i, j, k] (t = d1 + 2 d2 + 4 d3 the type), and for the entries of the mask in C order the exact
        values of Psi_b (Tb) and of Psi_1 ... Psi_r (T[k]) as lists of Fractions"""
        if hasattr(self, '_presc'):
            return self._presc
        n = self.n
        mask = np.zeros((self.nf, 8) + n, dtype=bool)
        for (i, side) in self.free_faces():
            sl = [slice(None)] * 3
            sl[i] = 0 if side == 0 else -1
            for (f, m) in self.relevant(i):
                for t in range(8):
                    if (t >> i) & 1 == m:
                        mask[(f, t) + tuple(sl)] = True
        psi_b, psis = self.states()
        entries = np.argwhere(mask)
        # pw[dir][node][d][g] = d^d/dx^d x^g at the node
        pw = [[[[(ff(g, d) * x ** (g - d) if g >= d else Fr(0)) for g in range(4)] for d in range(2)] for x in self.vnodes[k]]
              for k in range(3)]

        def values(v):
            byf = {}
            for c, x in v.items():
                byf.setdefault(c // 64, []).append(((c % 64) // 16, (c % 16) // 4, c % 4, x))
            out = []
            for (f, t, i, j, k) in entries:
                d = (t & 1, (t >> 1) & 1, (t >> 2) & 1)
                s = Fr(0)
                for (g1, g2, g3, x) in byf.get(int(f), ()):
                    y = pw[0][i][d[0]][g1]
                    if y == 0:
                        continue
                    y = y * pw[1][j][d[1]][g2] * pw[2][k][d[2]][g3]
                    if y != 0:
                        s += x * y
                out.append(s)
            return out

        self._presc = (mask, values(psi_b) if psi_b else [Fr(0)] * len(entries), [values(v) for v in psis])
        return self._presc

    def prescribed_enclosures(self):
        """float enclosures (lo, hi arrays over the entries of the mask) of Tb and of each T[k]"""
        if not hasattr(self, '_presc_enc'):
            mask, Tb, T = self.prescribed()
            cache = {}

            def enc(vals):
                lo, hi = np.empty(len(vals)), np.empty(len(vals))
                for q, v in enumerate(vals):
                    if v not in cache:
                        e = fr(v)
                        cache[v] = (float(e.lo), float(e.hi))
                    lo[q], hi[q] = cache[v]
                return lo, hi

            self._presc_enc = (enc(Tb), [enc(t) for t in T])
        return self._presc_enc


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


def _pairs(K):
    """the quadratic form sum K[k, l] X[k] X[l] with K exact and symmetric, as a list of rows
    (k, [(l, enclosure of the coefficient), ...]) over l >= k, the coefficient doubled for l > k"""
    keys = sorted({k for (k, l) in K} | {l for (k, l) in K})
    pos = {k: q for q, k in enumerate(keys)}
    rows = {}
    for (k, l), v in K.items():
        assert K.get((l, k)) == v, "the matrix is not symmetric"
        if pos[l] < pos[k]:
            continue
        rows.setdefault(k, []).append((l, fr(v if l == k else 2 * v)))
    return [(k, rows[k]) for k in keys if k in rows]


def _quadratic(rows, X):
    """enclosure, entry by entry, of  sum_{k, l} K[k, l] X[k] X[l]  (rows from _pairs, X a dict of interval arrays)"""
    acc = None
    for k, lst in rows:
        inner = None
        for l, c in lst:
            if l == k:
                term = X[k].sqr() * c
                acc = term if acc is None else acc + term
            else:
                term = X[l] * c
                inner = term if inner is None else inner + term
        if inner is not None:
            term = X[k] * inner
            acc = term if acc is None else acc + term
    return acc


# ================================================================================================================
# the displacement: candidate (heuristic) and its energy (RIGOROUS)
# ================================================================================================================
def _ucells(P, U):
    """for every cell, the nodal values (a, (i1, i2, i3)) as arrays over the cells"""
    p, N = P.p, P.N
    return {(a, (i, j, k)): U[a, i:i + p * (N[0] - 1) + 1:p, j:j + p * (N[1] - 1) + 1:p, k:k + p * (N[2] - 1) + 1:p]
            for a in range(3) for i in range(p + 1) for j in range(p + 1) for k in range(p + 1)}


def upper_bound(P, U):
    """RIGOROUS enclosure (lo, hi) of E[u_h], u_h the Q_p function with nodal values U (shape (3,) + P.M).
    E[u_h] = sum over the cells of  1/2 u^T K u - f.u  with the exact element matrix K and load vector f."""
    assert U.shape == (3,) + P.M and np.all(np.isfinite(U))
    assert np.array_equal(U[:, P.fixed], P.gval[:, P.fixed]), "Dirichlet nodal values are not imposed exactly"
    K, f, loc = P.stiffness_local()
    if not hasattr(P, '_Krows'):
        P._Krows = _pairs(K)
    X = {k: I(v) for k, v in _ucells(P, U).items()}
    en = _quadratic(P._Krows, X) * fr(Fr(1, 2))
    for k, v in f.items():
        if v != 0:
            en = en - X[k] * fr(v)
    tot = _total(en)
    return float(tot.lo), float(tot.hi)


def _stiffness_global(P):
    K, f, loc = P.stiffness_local()
    M1, M2, M3 = P.M
    p = P.p
    ci, cj, ck = [x.ravel() for x in np.meshgrid(*[np.arange(Ni) for Ni in P.N], indexing='ij')]

    def gidx(a, Pn):
        return ((a * M1 + (p * ci + Pn[0])) * M2 + (p * cj + Pn[1])) * M3 + (p * ck + Pn[2])

    nloc = 3 * len(loc)
    lidx = {(a, Pn): q for q, (a, Pn) in enumerate([(a, Pn) for a in range(3) for Pn in loc])}
    Kd = np.zeros((nloc, nloc))
    for (k, l), v in K.items():
        Kd[lidx[k], lidx[l]] = float(v)
    G = np.stack([gidx(a, Pn) for a in range(3) for Pn in loc], axis=1)            # cells x nloc
    nz = np.argwhere(Kd != 0)
    rows = G[:, nz[:, 0]].ravel()
    cols = G[:, nz[:, 1]].ravel()
    vals = np.tile(Kd[nz[:, 0], nz[:, 1]], G.shape[0])
    nd = 3 * M1 * M2 * M3
    Kg = sp.coo_matrix((vals, (rows, cols)), shape=(nd, nd)).tocsr()
    F = np.zeros(nd)
    for k, v in f.items():
        if v != 0:
            np.add.at(F, gidx(*k), float(v))
    return Kg, F


def solve_candidate(P, U0=None, rtol=1e-11):
    """heuristic: the finite-element solution in Q_p (any nodal array is a valid candidate)"""
    Kg, F = _stiffness_global(P)
    U = P.gval.copy().ravel()
    fixed = np.concatenate([P.fixed.ravel()] * 3)
    free = np.flatnonzero(~fixed)
    rhs = (F - Kg @ U)[free]
    Kff = Kg[free][:, free].tocsr()
    if free.size <= 20000:
        U[free] = spla.spsolve(Kff.tocsc(), rhs)
    else:
        d = Kff.diagonal()
        x0 = None if U0 is None else U0.ravel()[free]
        sol, info = spla.cg(Kff, rhs, x0=x0, rtol=rtol, maxiter=20000, M=sp.diags(1.0 / d))
        U[free] = sol
    return U.reshape((3,) + P.M)


# ================================================================================================================
# the stress field
# ================================================================================================================
class Stress:
    """X[f, t, i, j, k]: the nodal value of type t = d1 + 2 d2 + 4 d3 (the derivative d1^d1 d2^d2 d3^d3) of stress
    function f at node (i, j, k); the entries fixed by the traction-free faces (P.prescribed()) are ignored.
    par: the parameters P_1 ... P_r."""

    def __init__(self, P, X=None, par=None):
        mask, Tb, T = P.prescribed()
        shape = (P.nf, 8) + P.n
        self.X = np.zeros(shape) if X is None else np.array(X, dtype=float)
        self.par = np.zeros(len(T)) if par is None else np.array(par, dtype=float)
        assert self.X.shape == shape and self.par.shape == (len(T),)
        assert np.all(np.isfinite(self.X)) and np.all(np.isfinite(self.par))

    def copy(self, P):
        return Stress(P, self.X.copy(), self.par.copy())


def nodal_enclosure(P, st):
    """RIGOROUS: the nodal values of the stress functions as an interval array.  Free entries: the stored floats
    (exact).  Prescribed entries: enclosures of the exact values of Psi_b + sum_k P_k Psi_k."""
    mask, Tb, T = P.prescribed()
    lo, hi = st.X.copy(), st.X.copy()
    if mask.any():
        (blo, bhi), encT = P.prescribed_enclosures()
        v = I(blo, bhi)
        for k in range(len(T)):
            v = v + I(st.par[k]) * I(*encT[k])
        lo[mask], hi[mask] = v.lo * np.ones(blo.shape), v.hi * np.ones(blo.shape)
    return I(lo, hi)


def _xcells(P, X):
    """for every cell, the nodal values (f, (k1, k2, k3)), k_i = 2 e_i + d_i, as (interval) arrays over the cells"""
    N = P.N
    out = {}
    for f in range(P.nf):
        for k in itertools.product(range(4), repeat=3):
            t = (k[0] % 2) + 2 * (k[1] % 2) + 4 * (k[2] % 2)
            out[(f, k)] = X[f, t, (k[0] // 2):(k[0] // 2) + N[0], (k[1] // 2):(k[1] // 2) + N[1], (k[2] // 2):(k[2] // 2) + N[2]]
    return out


def _face_slices(P, i, side):
    """helpers for the face (i, side): tangential directions, and functions giving the face-cell arrays"""
    t1, t2 = [t for t in range(3) if t != i]
    pos = 0 if side == 0 else P.N[i]
    gpos = 0 if side == 0 else P.M[i] - 1
    N, p = P.N, P.p

    def xface(X, f, m, k1, k2):
        """nodal values of d_i^m Phi_f of tangential local index (k1, k2), over the face cells"""
        d = [0, 0, 0]
        d[i], d[t1], d[t2] = m, k1 % 2, k2 % 2
        t = d[0] + 2 * d[1] + 4 * d[2]
        sl = [None, None, None]
        sl[i] = pos
        sl[t1] = slice(k1 // 2, k1 // 2 + N[t1])
        sl[t2] = slice(k2 // 2, k2 // 2 + N[t2])
        return X[(f, t) + tuple(sl)]

    def gface(G, a, q1, q2):
        """nodal values of g_a at the local node (q1, q2) of the face cells"""
        sl = [None, None, None]
        sl[i] = gpos
        sl[t1] = slice(q1, q1 + p * (N[t1] - 1) + 1, p)
        sl[t2] = slice(q2, q2 + p * (N[t2] - 1) + 1, p)
        return G[(a,) + tuple(sl)]

    return t1, t2, xface, gface


def lower_bound(P, st, return_parts=False):
    """RIGOROUS lower bound  int_{Gamma_D} (sigma nu).g - 1/2 int s^T S s  for the stress field st.
    All coefficients are exact rationals; the nodal values are exact floats or enclosures (nodal_enclosure)."""
    A, w, loc = P.stress_local()
    HL, LI = P.face_local()
    if not hasattr(P, '_Arows'):
        P._Arows = _pairs(A)
    X = nodal_enclosure(P, st)
    Xc = _xcells(P, X)
    # the cell integrals:  int s^T S s = Q_phi + 2 Q_cross + Q_b
    q = _quadratic(P._Arows, Xc)
    if q is None:
        q = I(np.zeros(P.N))
    if any(x != 0 for x in P.b):
        shp = [(-1, 1, 1), (1, -1, 1), (1, 1, -1)]
        xc = []
        for n in range(3):                      # the lower corners of the cells are rationals; enclose them
            enc = [fr(P.vnodes[n][j]) for j in range(P.N[n])]
            xc.append(I(np.array([float(e.lo) for e in enc]).reshape(shp[n]) * np.ones(P.N),
                        np.array([float(e.hi) for e in enc]).reshape(shp[n]) * np.ones(P.N)))
        for k in loc:
            cf = [w[r].get(k, Fr(0)) for r in range(4)]
            if any(c != 0 for c in cf):
                lin = fr(cf[0]) + xc[0] * fr(cf[1]) + xc[1] * fr(cf[2]) + xc[2] * fr(cf[3])
                q = q + Xc[k] * lin * fr(2)
    S, b = P.S, P.b
    L = [c - a for a, c in P.box]
    m1 = [(c ** 2 - a ** 2) / 2 for a, c in P.box]
    m2 = [(c ** 3 - a ** 3) / 3 for a, c in P.box]
    qb = Fr(0)
    for i in range(3):
        for j in range(3):
            if i == j:
                mom = m2[i] * L[(i + 1) % 3] * L[(i + 2) % 3]
            else:
                mom = m1[i] * m1[j] * L[3 - i - j]
            qb += S[i][j] * b[i] * b[j] * mom
    intQ = _total(q) + fr(qb)
    # the boundary term
    bnd = I(0.0)
    G = I(P.gval)
    for i, side, sign in P.dirichlet_faces():
        a = P.box[i][side]
        t1, t2, xface, gface = _face_slices(P, i, side)
        h1, h2 = P.h[t1], P.h[t2]
        acc = None
        for j in range(3):
            gq = {(q1, q2): gface(G, j, q1, q2) for q1 in range(P.p + 1) for q2 in range(P.p + 1)}
            for (f, al, c) in P.stress[VOIGT[i][j]]:
                for k1 in range(4):
                    for k2 in range(4):
                        comb = None
                        for (q1, q2), gv in gq.items():
                            cf = c * h1 ** (1 + k1 % 2 - al[t1]) * h2 ** (1 + k2 % 2 - al[t2]) * HL[al[t1]][k1][q1] * HL[al[t2]][k2][q2]
                            if cf != 0:
                                term = gv * fr(cf)
                                comb = term if comb is None else comb + term
                        if comb is not None:
                            term = xface(X, f, al[i], k1, k2) * comb
                            acc = term if acc is None else acc + term
            if j == i and P.b[i] != 0 and a != 0:
                for (q1, q2), gv in gq.items():
                    cf = -P.b[i] * a * h1 * h2 * LI[q1] * LI[q2]
                    if cf != 0:
                        term = gv * fr(cf)
                        acc = term if acc is None else acc + term
        if acc is not None:
            bnd = bnd + _total(acc) * fr(sign)
    lb = bnd - intQ * fr(Fr(1, 2))
    if return_parts:
        return float(lb.lo), dict(boundary=(float(bnd.lo), float(bnd.hi)), intQ=(float(intQ.lo), float(intQ.hi)))
    return float(lb.lo)


# ================================================================================================================
# the heuristic for the stress field: maximise the lower bound (floating point, not trusted)
# ================================================================================================================
def _dual_setup(P):
    """the floating-point data of the dual problem: LB = beta.x - wg.x - 1/2 x^T A x + const over the nodal values x"""
    if hasattr(P, '_dual'):
        return P._dual
    A, w, loc = P.stress_local()
    HL, LI = P.face_local()
    nf, n, N = P.nf, P.n, P.N
    nd = nf * 8 * n[0] * n[1] * n[2]
    nloc = len(loc)
    lidx = {k: q for q, k in enumerate(loc)}
    Af = np.zeros((nloc, nloc))
    for (k, l), v in A.items():
        Af[lidx[k], lidx[l]] = float(v)

    def gid(f, t, i, j, k):
        return (((f * 8 + t) * n[0] + i) * n[1] + j) * n[2] + k

    ci, cj, ck = [x.ravel() for x in np.meshgrid(*[np.arange(Ni) for Ni in N], indexing='ij')]
    cidx = np.empty((ci.size, nloc), dtype=np.int64)
    for (f, k), q in lidx.items():
        t = (k[0] % 2) + 2 * (k[1] % 2) + 4 * (k[2] % 2)
        cidx[:, q] = gid(f, t, ci + k[0] // 2, cj + k[1] // 2, ck + k[2] // 2)
    wg = np.zeros(nd)
    if any(x != 0 for x in P.b):
        xc = [np.array([float(P.vnodes[0][i]) for i in ci]), np.array([float(P.vnodes[1][j]) for j in cj]),
              np.array([float(P.vnodes[2][k]) for k in ck])]
        for k in loc:
            cf = [float(w[r].get(k, 0)) for r in range(4)]
            if any(cf):
                np.add.at(wg, cidx[:, lidx[k]], cf[0] + xc[0] * cf[1] + xc[1] * cf[2] + xc[2] * cf[3])
    beta = np.zeros(nd)
    for i, side, sign in P.dirichlet_faces():
        t1, t2 = [t for t in range(3) if t != i]
        pos = 0 if side == 0 else N[i]
        gpos = 0 if side == 0 else P.M[i] - 1
        h1, h2 = float(P.h[t1]), float(P.h[t2])
        c1, c2 = [x.ravel() for x in np.meshgrid(np.arange(N[t1]), np.arange(N[t2]), indexing='ij')]
        for j in range(3):
            Gj = np.take(P.gval[j], gpos, axis=i)                        # (M_t1, M_t2)
            for (f, al, c) in P.stress[VOIGT[i][j]]:
                for k1 in range(4):
                    for k2 in range(4):
                        comb = 0.0
                        for q1 in range(P.p + 1):
                            for q2 in range(P.p + 1):
                                cf = float(c * HL[al[t1]][k1][q1] * HL[al[t2]][k2][q2]) * h1 ** (1 + k1 % 2 - al[t1]) * h2 ** (1 + k2 % 2 - al[t2])
                                if cf:
                                    comb = comb + cf * Gj[P.p * c1 + q1, P.p * c2 + q2]
                        if np.ndim(comb) == 0:
                            continue
                        d = [0, 0, 0]
                        d[i], d[t1], d[t2] = al[i], k1 % 2, k2 % 2
                        nidx = [None, None, None]
                        nidx[i], nidx[t1], nidx[t2] = pos, c1 + k1 // 2, c2 + k2 // 2
                        np.add.at(beta, gid(f, d[0] + 2 * d[1] + 4 * d[2], nidx[0], nidx[1], nidx[2]), sign * comb)
    mask, Tb, T = P.prescribed()
    isp = mask.ravel()
    x0 = np.zeros(nd)
    Tm = np.zeros((len(T), nd))
    if isp.any():
        x0[isp] = [float(v) for v in Tb]
        for k in range(len(T)):
            Tm[k, isp] = [float(v) for v in T[k]]
    # scale of each nodal value: h^d (the nodal values are derivatives)
    hh = [float(x) for x in P.h]
    sc = np.ones((nf, 8) + n)
    for t in range(8):
        sc[:, t] = hh[0] ** (t & 1) * hh[1] ** ((t >> 1) & 1) * hh[2] ** ((t >> 2) & 1)
    P._dual = dict(Af=Af, cidx=cidx, wg=wg, beta=beta, x0=x0, Tm=Tm, free=np.flatnonzero(~isp), nd=nd, scale=sc.ravel())
    return P._dual


def _apply_A(d, x):
    """A x for the assembled matrix of the cell integrals"""
    y = x[d['cidx']] @ d['Af']
    return np.bincount(d['cidx'].ravel(), weights=y.ravel(), minlength=d['nd'])


def dual_value(P, st):
    """floating-point value of the lower bound (for monitoring only)"""
    d = _dual_setup(P)
    x = st.X.ravel().copy()
    isp = np.ones(d['nd'], dtype=bool)
    isp[d['free']] = False
    x[isp] = d['x0'][isp] + st.par @ d['Tm'][:, isp] if len(st.par) else d['x0'][isp]
    return float((d['beta'] - d['wg']) @ x - 0.5 * x @ _apply_A(d, x))


def solve_stress(P, st0=None, rtol=1e-9, maxiter=20000, eps=1e-10, method='auto', verbose=False):
    """heuristic: a maximiser of the lower bound over the free nodal values and the parameters.  The system is
    singular (docs/ELASTICITY_3D.md, Section 3.5); a small multiple of the identity (in scaled unknowns) is added.
    method: 'dense' (least squares), 'cg' (preconditioned conjugate gradients), 'auto'."""
    d = _dual_setup(P)
    free, nd, Tm, sc = d['free'], d['nd'], d['Tm'], d['scale']
    r = Tm.shape[0]
    nz = free.size + r
    lin = d['beta'] - d['wg'] - _apply_A(d, d['x0'])
    # unknowns z = (y, P): x = x0 + scale * y on the free entries, + sum P_k Tm[k]
    tnorm = np.array([max(np.max(np.abs(Tm[k])), 1e-300) for k in range(r)])

    def to_x(z):
        x = np.zeros(nd)
        x[free] = sc[free] * z[:free.size]
        if r:
            x += (z[free.size:] / tnorm) @ Tm
        return x

    def from_x(g):
        out = np.empty(nz)
        out[:free.size] = sc[free] * g[free]
        if r:
            out[free.size:] = (Tm @ g) / tnorm
        return out

    rhs = from_x(lin)
    if method == 'auto':
        method = 'dense' if nz <= 3000 else 'cg'
    if method == 'dense':
        H = np.empty((nz, nz))
        for q in range(nz):
            e = np.zeros(nz)
            e[q] = 1.0
            H[:, q] = from_x(_apply_A(d, to_x(e)))
        H = 0.5 * (H + H.T)
        wv, V = sla.eigh(H)
        keep = wv > 1e-11 * wv[-1]
        z = V[:, keep] @ ((V[:, keep].T @ rhs) / wv[keep])
        info = dict(method='dense', rank=int(keep.sum()), unknowns=nz)
    else:
        # diagonal of the operator in the unknowns z
        dA = np.bincount(d['cidx'].ravel(), weights=np.tile(np.diag(d['Af']), d['cidx'].shape[0]), minlength=nd)
        diag = np.empty(nz)
        diag[:free.size] = sc[free] ** 2 * dA[free]
        for k in range(r):
            diag[free.size + k] = (Tm[k] / tnorm[k]) @ _apply_A(d, Tm[k] / tnorm[k])
        assert np.all(diag > 0)
        its = [0]

        def mv(z):                              # the regularisation is a multiple of the diagonal
            its[0] += 1
            return from_x(_apply_A(d, to_x(z))) + eps * diag * z

        op = spla.LinearOperator((nz, nz), matvec=mv, dtype=float)
        z0 = None
        if st0 is not None:
            z0 = np.empty(nz)
            z0[:free.size] = st0.X.ravel()[free] / sc[free]
            z0[free.size:] = st0.par * tnorm
        z, flag = spla.cg(op, rhs, x0=z0, rtol=rtol, maxiter=maxiter, M=sp.diags(1.0 / diag))
        info = dict(method='cg', iterations=its[0], flag=int(flag), unknowns=nz)
    X = np.zeros(nd)
    X[free] = sc[free] * z[:free.size]
    st = Stress(P, X.reshape((P.nf, 8) + P.n), z[free.size:] / tnorm if r else None)
    if verbose:
        print('   solve_stress:', info, flush=True)
    st.info = info
    return st


# ================================================================================================================
# drivers and certificates on disk
# ================================================================================================================
def certify(P, U=None, st=None, **kw):
    """both bounds: candidate (finite-element solution unless given) and stress field (dual solve unless given)"""
    t0 = time.time()
    U = solve_candidate(P) if U is None else U
    t1 = time.time()
    st = solve_stress(P, **kw) if st is None else st
    t2 = time.time()
    ub = upper_bound(P, U)
    lb = lower_bound(P, st)
    t3 = time.time()
    return dict(N=list(P.N), degree=P.p, lower=lb, upper=ub[1], upper_enclosure=list(ub),
                rel_gap=(ub[1] - lb) / abs(ub[1]) if ub[1] != 0 else float('nan'), seconds=t3 - t0,
                seconds_candidate=t1 - t0, seconds_stress=t2 - t1, seconds_bounds=t3 - t2), U, st


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
    meta = dict(format='CertiField-elasticity3d-certificate-1', problem=P.spec(), lower=lo, upper=hi, reference=__version__)
    if extra:
        meta['info'] = extra
    json.dump(meta, open(stem + '.json', 'w', newline='\n'), indent=1)
    mask = P.prescribed()[0]
    X = np.where(mask, 0.0, st.X)                       # the prescribed entries are not data
    np.savez_compressed(stem + '.npz', gval=P.gval, U=U, X=X, par=st.par)


def load_certificate(stem):
    meta = json.load(open(stem + '.json'))
    z = np.load(stem + '.npz')
    P = Problem3D.from_spec(meta['problem'], z['gval'])
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
    """raw export for a port: little-endian float64, C order, plus field.json.
    X.f64 has the shape (functions, 8, n1, n2, n3): X[f, t, i, j, k] is the nodal value of type t = d1 + 2 d2 + 4 d3 of
    function f at node (i, j, k); the entries fixed by the traction-free faces are stored as zero and are not data."""
    import os
    os.makedirs(d, exist_ok=True)
    mask, Tb, T = P.prescribed()
    meta = dict(format='CertiField-elasticity3d-field-1', problem=P.spec(), NodeDims=list(P.M), GridDims=list(P.n),
                Functions=[FNAMES[f] for f in P.funcs], Parameters=len(st.par), Prescribed=int(mask.sum()))
    if expected is not None:
        meta['expected'] = expected
    json.dump(meta, open(os.path.join(d, 'field.json'), 'w', newline='\n'), indent=1)
    np.asarray(P.gval, dtype='<f8').tofile(os.path.join(d, 'gval.f64'))
    np.asarray(U, dtype='<f8').tofile(os.path.join(d, 'U.f64'))
    np.asarray(np.where(mask, 0.0, st.X), dtype='<f8').tofile(os.path.join(d, 'X.f64'))
    np.asarray(st.par, dtype='<f8').tofile(os.path.join(d, 'par.f64'))
