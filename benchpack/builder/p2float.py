"""P2 (L-bracket, full domain): spaces, float assembly and solves (not trusted; they only FIND the fields).

Fields (certificates/FORMAT.md section 4):
    v   = (0, -(x+1)/2) + tensor B-splines (A, B), A, B zero on the first and last x-index;
    tau = Airy(phi_p + Psi), phi_p = chi(y)(alpha + beta x + gamma y), Psi tensor B-splines with zero sets.
"""
import os
import sys
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', 'problems'))
import definitions as D                                    # noqa: E402
from splines import Spline1D, open_knots, gauss01           # noqa: E402
from p1float import solve_spd                               # noqa: E402

C11, C12, C66 = (float(D.C[0][0]), float(D.C[0][1]), float(D.C[2][2]))
NU = float(D.NU)
C1111 = D.C1111
PATCH = D.P2['patch']                                        # ((0, 1/10), (0, 1/10))
PATCH_AREA = D.P2['patch_area']                              # 1/100
# chi(y) = 1 - 3y^2 + 2y^3 on [0, 1] (1 for y <= 0)
CHI = [Fr(1), Fr(0), Fr(-3), Fr(2)]


def graded_breaks(lo, hi, h, targets, extra=()):
    """sorted rational breakpoints: uniform spacing h, plus geometric layers towards the targets.
    targets: list of (c, sides, L, m): points c +- h 2^-k (1 + j/m), k = 0..L-1, j = 0..m-1, on the given sides"""
    lo, hi, h = Fr(lo), Fr(hi), Fr(h)
    pts = set()
    n = int((hi - lo) / h)
    assert lo + n * h == hi
    pts.update(lo + k * h for k in range(n + 1))
    for (c, sides, L, m) in targets:
        c = Fr(c)
        for k in range(L):
            base = h / 2 ** (k + 1)
            for j in range(m):
                d = base * (1 + Fr(j, m))
                for s in sides:
                    p = c + s * d
                    if lo < p < hi:
                        pts.add(p)
    pts.update(Fr(e) for e in extra)
    return sorted(pts)


def classes(S):
    """index classes of a Spline1D on [-1, 1] relative to 0: 'L' support in [-1,0], 'R' in [0,1], 'S' otherwise"""
    t = S.t
    p = S.p
    out = []
    for i in range(S.n):
        a, b = t[i], t[i + p + 1]
        out.append('L' if b <= 0 else ('R' if a >= 0 else 'S'))
    return np.array(out)


class Setup2:
    def __init__(self, sbx, sby, ps, ubx, uby, pu, smult=None, umult=None):
        self.sx = Spline1D(open_knots(sbx, ps, (smult or {}).get('x')), ps)
        self.sy = Spline1D(open_knots(sby, ps, (smult or {}).get('y')), ps)
        self.ux = Spline1D(open_knots(ubx, pu, (umult or {}).get('x')), pu)
        self.uy = Spline1D(open_knots(uby, pu, (umult or {}).get('y')), pu)
        self.ps, self.pu = ps, pu
        for Sx in (self.sx, self.sy):
            assert sum(1 for k in Sx.t if k == 0) == ps - 1, "knot 0 must have multiplicity ps - 1 in the stress space"
        for Sx in (self.sx, self.sy, self.ux, self.uy):
            assert 0 in Sx.t and Fr(1, 10) in Sx.t

    def cells(self, X, Y):
        """cells of Omega (cells in [0,1] x [-1,0] are outside)"""
        return [(ex, ey) for ex in range(X.ne) for ey in range(Y.ne) if not (X.a[ex] >= 0 and Y.b[ey] <= 0)]

    def patch_cells(self, X, Y):
        p = PATCH
        return [(ex, ey) for ex in range(X.ne) for ey in range(Y.ne)
                if X.a[ex] >= p[0][0] and X.b[ex] <= p[0][1] and Y.a[ey] >= p[1][0] and Y.b[ey] <= p[1][1]]

    # ------------------------------------------------------------------ bases
    def _tables(self, X, Y, cells, nq, nder):
        cells = np.array(cells, dtype=np.int64).reshape(-1, 2)
        g, w = gauss01(nq)
        ex, ey = cells[:, 0], cells[:, 1]
        bx = X.eval_float(np.repeat(ex[:, None], nq, 1), np.broadcast_to(g, (len(ex), nq)), nder)   # [nb, nq, p+1, d]
        by = Y.eval_float(np.repeat(ey[:, None], nq, 1), np.broadcast_to(g, (len(ey), nq)), nder)
        x = X.af[ex][:, None] + X.hf[ex][:, None] * g[None, :]
        y = Y.af[ey][:, None] + Y.hf[ey][:, None] * g[None, :]
        W = (w[None, :, None] * w[None, None, :]) * (X.hf[ex] * Y.hf[ey])[:, None, None]
        gx = X.first[ex][:, None] + np.arange(X.p + 1)[None, :]
        gy = Y.first[ey][:, None] + np.arange(Y.p + 1)[None, :]
        gidx = (gx[:, :, None] * Y.n + gy[:, None, :]).reshape(len(ex), -1)
        return ex, ey, bx, by, x, y, W.reshape(len(ex), -1), gidx

    @staticmethod
    def _tp(bx, by, dx, dy):
        nb, nq = bx.shape[:2]
        return np.einsum('bim,bjn->bijmn', bx[..., dx], by[..., dy]).reshape(nb, nq * nq, -1)

    def dual_basis(self, cells, nq):
        """T [nb, nq^2, 3, nloc] stress of local Psi functions, Tp [nb, nq^2, 3, 3] stress of (alpha, beta, gamma)"""
        X, Y = self.sx, self.sy
        ex, ey, bx, by, x, y, W, gidx = self._tables(X, Y, cells, nq, 2)
        nb = len(ex)
        Pxx, Pyy, Pxy = self._tp(bx, by, 2, 0), self._tp(bx, by, 0, 2), self._tp(bx, by, 1, 1)
        T = np.stack([Pyy, Pxx, -Pxy], axis=2)
        X2 = np.broadcast_to(x[:, :, None], (nb, nq, nq)).reshape(nb, -1)
        Y2 = np.broadcast_to(y[:, None, :], (nb, nq, nq)).reshape(nb, -1)
        up = Y2 >= 0
        c = [float(v) for v in CHI]
        chi = np.where(up, c[0] + c[1] * Y2 + c[2] * Y2 ** 2 + c[3] * Y2 ** 3, 1.0)
        d1 = np.where(up, c[1] + 2 * c[2] * Y2 + 3 * c[3] * Y2 ** 2, 0.0)
        d2 = np.where(up, 2 * c[2] + 6 * c[3] * Y2, 0.0)
        # phi_p = chi (alpha + beta x + gamma y): tau_xx = chi'' L + 2 chi' gamma, tau_yy = 0, tau_xy = -chi' beta
        Tp = np.zeros((nb, nq * nq, 3, 3))
        Tp[:, :, 0, 0] = d2
        Tp[:, :, 0, 1] = d2 * X2
        Tp[:, :, 0, 2] = d2 * Y2 + 2 * d1
        Tp[:, :, 2, 1] = -d1
        return T, Tp, W, gidx

    def primal_basis(self, cells, nq):
        X, Y = self.ux, self.uy
        ex, ey, bx, by, x, y, W, gidx = self._tables(X, Y, cells, nq, 1)
        Nx, Ny = self._tp(bx, by, 1, 0), self._tp(bx, by, 0, 1)
        Z = np.zeros_like(Nx)
        Bm = np.stack([np.concatenate([Nx, Z], -1), np.concatenate([Z, Ny], -1), np.concatenate([Ny, Nx], -1)], axis=2)
        nA = X.n * Y.n
        return Bm, W, np.concatenate([gidx, gidx + nA], axis=1)

    # ------------------------------------------------------------------ assembly
    def assemble_dual(self, batch=256):
        """Q(tau) = z^T A z for z = (alpha, beta, gamma, C.ravel());  F = f.z;  l(tau) = lt.z"""
        X, Y = self.sx, self.sy
        nC = X.n * Y.n
        Smat = np.array([[1, -NU, 0], [-NU, 1, 0], [0, 0, 2 * (1 + NU)]])
        nq = self.ps + 1
        cells = self.cells(X, Y)
        rows, cols, vals = [], [], []
        App = np.zeros((3, 3))
        Apc = np.zeros((3, nC))
        for s in range(0, len(cells), batch):
            T, Tp, W, gidx = self.dual_basis(cells[s:s + batch], nq)
            ST = np.einsum('cd,bqdn->bqcn', Smat, T)
            Al = np.einsum('bq,bqcm,bqcn->bmn', W, T, ST)
            nb, nl = gidx.shape
            rows.append(np.repeat(gidx, nl, axis=1).ravel())
            cols.append(np.tile(gidx, (1, nl)).ravel())
            vals.append(Al.ravel())
            cr = np.einsum('bq,bqck,bqcn->bkn', W, Tp, ST)
            for b in range(nb):
                np.add.at(Apc, (slice(None), gidx[b]), cr[b])
            App += np.einsum('bq,bqck,cd,bqdl->kl', W, Tp, Smat, Tp)
        Acc = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nC, nC)).tocsr()
        A = sp.bmat([[sp.csr_matrix(App), sp.csr_matrix(Apc)], [sp.csr_matrix(Apc.T), Acc]], format='csr')
        f = np.zeros(3 + nC)
        f[1] = -1.0                                         # F = -beta
        # patch functional
        lt = np.zeros(3 + nC)
        pc = self.patch_cells(X, Y)
        T, Tp, W, gidx = self.dual_basis(pc, nq)
        lt[:3] = np.einsum('bq,bqk->k', W, Tp[:, :, 0, :])
        np.add.at(lt, (gidx + 3).ravel(), np.einsum('bq,bqn->bn', W, T[:, :, 0, :]).ravel())
        return A, f, lt / float(PATCH_AREA)

    def assemble_primal(self, batch=512):
        """a(v,v) = z^T K z + 2 b.z + k0 for z = (A.ravel(), B.ravel()); lift (0, -(x+1)/2); l(v) = l0 + lv.z"""
        X, Y = self.ux, self.uy
        n = 2 * X.n * Y.n
        Cm = np.array([[C11, C12, 0], [C12, C11, 0], [0, 0, C66]])
        e0 = np.array([0.0, 0.0, -0.5])                     # strain of the lift: gamma_xy = d/dx(-(x+1)/2) = -1/2
        nq = self.pu + 1
        cells = self.cells(X, Y)
        rows, cols, vals = [], [], []
        b = np.zeros(n)
        k0 = 0.0
        for s in range(0, len(cells), batch):
            Bm, W, gidx = self.primal_basis(cells[s:s + batch], nq)
            CB = np.einsum('cd,bqdn->bqcn', Cm, Bm)
            Kl = np.einsum('bq,bqcm,bqcn->bmn', W, Bm, CB)
            nb, nl = gidx.shape
            rows.append(np.repeat(gidx, nl, axis=1).ravel())
            cols.append(np.tile(gidx, (1, nl)).ravel())
            vals.append(Kl.ravel())
            np.add.at(b, gidx.ravel(), np.einsum('bq,c,bqcn->bn', W, Cm @ e0, Bm).ravel())
            k0 += float(np.sum(W)) * float(e0 @ Cm @ e0)
        K = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
        lv = np.zeros(n)
        Bm, W, gidx = self.primal_basis(self.patch_cells(X, Y), nq)
        sxx = C11 * Bm[:, :, 0, :] + C12 * Bm[:, :, 1, :]
        np.add.at(lv, gidx.ravel(), np.einsum('bq,bqn->bn', W, sxx).ravel())
        l0 = 0.0                                            # the lift has eps_xx = eps_yy = 0
        return K, b, k0, l0, lv / float(PATCH_AREA)

    # ------------------------------------------------------------------ ties
    def stress_free(self):
        """boolean mask [nx, ny] of the Psi coefficients that are NOT forced to zero (FORMAT.md 4.2)"""
        cx, cy = classes(self.sx), classes(self.sy)
        nx, ny = self.sx.n, self.sy.n
        free = np.ones((nx, ny), dtype=bool)
        free[:, ny - 2:] = False
        LSx = (cx == 'L') | (cx == 'S')
        free[np.ix_(LSx, np.array([True, True] + [False] * (ny - 2)))] = False
        free[np.ix_(cx == 'S', (cy == 'L') | (cy == 'S'))] = False
        free[np.ix_(cx == 'R', (cy == 'S') | (cy == 'L'))] = False
        return free

    def stress_ties(self):
        free = self.stress_free().ravel()
        idx = np.flatnonzero(free)
        Pc = sp.csr_matrix((np.ones(len(idx)), (idx, np.arange(len(idx)))), shape=(free.size, len(idx)))
        return sp.block_diag([sp.identity(3, format='csr'), Pc], format='csr')

    def disp_ties(self):
        nx, ny = self.ux.n, self.uy.n
        Ix, Iy = sp.identity(nx, format='csr'), sp.identity(ny, format='csr')
        PA = sp.kron(Ix[:, 1:nx - 1], Iy, format='csr')
        return sp.block_diag([PA, PA], format='csr')
