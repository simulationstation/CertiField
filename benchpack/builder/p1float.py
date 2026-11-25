"""P1 (plate with a hole, quarter model): spaces, float quadrature, assembly and solves.

Nothing in this file needs to be rigorous: it only FINDS the fields.  The certified numbers come from p1rig.py,
which re-evaluates the stored fields with interval / exact arithmetic.  See benchpack/certificates/FORMAT.md for the
definition of the fields and benchpack/NOTES.md for the argument.
"""
import math
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', 'problems'))
import definitions as D                                    # noqa: E402
from splines import Spline1D, open_knots, gauss01          # noqa: E402

A_HOLE = Fr(1, 2)
R2 = Fr(11, 20)
THETA_HALF = Fr(1, 10)
PH_AREA = D.P1_QUARTER['half_patch_area']                 # 21/8000
C11, C12, C66 = (float(D.C[0][0]), float(D.C[0][1]), float(D.C[2][2]))
NU = float(D.NU)
C1111 = D.C1111

# rho = eta(x) + zeta(y)
ETA = [Fr(0), Fr(0), Fr(2), Fr(0), Fr(-1)]                 # 2x^2 - x^4
ZETA = [Fr(-7, 16), Fr(0), Fr(3, 2), Fr(0), Fr(1)]         # y^4 + 3/2 y^2 - 7/16
assert all(D.RHO.get((i, 0), 0) == ETA[i] for i in range(1, 5)) and all(D.RHO.get((0, j), 0) == ZETA[j] for j in range(1, 5))
assert D.RHO[(0, 0)] == ZETA[0]


def _pv(c, x, d=0):
    """float value of the d-th derivative of the polynomial with Fraction coefficients c"""
    c = [float(v) for v in c]
    for _ in range(d):
        c = [k * c[k] for k in range(1, len(c))] or [0.0]
    r = np.zeros_like(x, dtype=float)
    for v in reversed(c):
        r = r * x + v
    return r


class Setup:
    """spaces and discretisation parameters of one P1 certificate.
    sx, sy: Spline1D of the stress space (degree ps);  ux, uy: of the displacement space (degree pu);
    gexp: exponents of the particular field g;  hole: False gives the plate WITHOUT hole (known-answer test)"""

    def __init__(self, sbx, sby, ps, ubx, uby, pu, gexp=None, gbase=Fr(1, 2), hole=True, grip='frictionless', smult=None,
                 umult=None):
        sbx, sby, ubx, uby = [[Fr(b) for b in v] for v in (sbx, sby, ubx, uby)]
        self.gbase = Fr(gbase)
        assert self.gbase in sby and self.gbase >= Fr(1, 2), "g_base must be a knot of the stress space in y, >= 1/2"
        self.sx = Spline1D(open_knots(sbx, ps, (smult or {}).get('x')), ps)
        self.sy = Spline1D(open_knots(sby, ps, (smult or {}).get('y')), ps)
        self.ux = Spline1D(open_knots(ubx, pu, (umult or {}).get('x')), pu)
        self.uy = Spline1D(open_knots(uby, pu, (umult or {}).get('y')), pu)
        self.ps, self.pu = ps, pu
        self.gexp = tuple(int(e) for e in (gexp if gexp is not None else (ps - 2, ps - 1, ps)))
        self.hole = hole
        self.grip = grip
        self.nqs = ps + 9          # tau has degree <= ps + 8 per variable: integrand <= 2 ps + 16 -> ps + 9 Gauss points
        self.nqu = pu + 1          # strain degree <= pu: integrand <= 2 pu -> pu + 1 points
        self.theta0 = math.pi / 2 - 0.1

    # ---------------------------------------------------------------------------------------------- geometry
    def cells(self, X, Y):
        """list of (ex, ey, kind) with kind 0 = inside Omega_q, 1 = cut by the circle (cells in the disc dropped)"""
        out = []
        a2 = A_HOLE * A_HOLE
        for ex in range(X.ne):
            for ey in range(Y.ne):
                if not self.hole:
                    out.append((ex, ey, 0))
                    continue
                x0, x1, y0, y1 = X.a[ex], X.b[ex], Y.a[ey], Y.b[ey]
                if x0 * x0 + y0 * y0 >= a2:
                    out.append((ex, ey, 0))
                elif x1 * x1 + y1 * y1 > a2:
                    out.append((ex, ey, 1))
        return out

    def cut_rule(self, X, Y, ex, ey, ny, nt):
        """float quadrature of (cell minus disc) directly (no subtraction): local coordinates and weights"""
        hx, hy = X.hf[ex], Y.hf[ey]
        x0, y0 = X.af[ex], Y.af[ey]
        x1, y1 = x0 + hx, y0 + hy
        g, w = gauss01(ny)
        a = float(A_HOLE)
        xs = math.sqrt(a * a - y1 * y1) if y1 < a else 0.0          # Y(x) >= y1 for x <= xs: nothing there
        X0 = math.sqrt(a * a - y0 * y0)                               # Y(x) <= y0 for x >= X0: full height
        pts_x, pts_y, ws = [], [], []
        if X0 < x1:                                                   # rectangle [max(x0,X0), x1] x [y0, y1]
            xa = max(x0, X0)
            SX, SY = np.meshgrid(g, g, indexing='ij')
            pts_x.append(xa + (x1 - xa) * SX.ravel())
            pts_y.append(y0 + hy * SY.ravel())
            ws.append(np.outer(w, w).ravel() * (x1 - xa) * hy)
        xa, xb = max(x0, xs), min(x1, X0)
        if xa < xb:                                                   # y from sqrt(a^2 - x^2) to y1
            th_hi, th_lo = math.acos(xa / a), math.acos(xb / a)
            gt, wt = gauss01(nt)
            th = th_lo + (th_hi - th_lo) * gt
            xx = a * np.cos(th)
            ybot = a * np.sin(th)
            jac = a * np.sin(th) * (th_hi - th_lo) * wt
            pts_x.append(np.repeat(xx, ny))
            pts_y.append((ybot[:, None] + (y1 - ybot[:, None]) * g[None, :]).ravel())
            ws.append((jac[:, None] * (y1 - ybot[:, None]) * w[None, :]).ravel())
        px, py, ww = np.concatenate(pts_x), np.concatenate(pts_y), np.concatenate(ws)
        return (px - x0) / hx, (py - y0) / hy, ww

    def rules(self, X, Y, nq, nt=None):
        """per cell: (ex, ey, sx [nq], sy [nq], w [nq]) grouped: full cells batched, cut cells one by one"""
        nt = nt or 2 * nq + 10
        g, w = gauss01(nq)
        SX, SY = np.meshgrid(g, g, indexing='ij')
        full, cut = [], []
        for ex, ey, kind in self.cells(X, Y):
            if kind == 0:
                full.append((ex, ey))
            else:
                cut.append((ex, ey) + self.cut_rule(X, Y, ex, ey, nq, nt))
        full = np.array(full, dtype=np.int64).reshape(-1, 2)
        wf = np.outer(w, w).ravel()
        return dict(full=full, sx=SX.ravel(), sy=SY.ravel(), wf=wf, cut=cut)

    # ---------------------------------------------------------------------------------------------- stress basis
    def dual_basis(self, ex, ey, sx, sy):
        """at points of cells ex, ey [nb] with local coordinates sx, sy [nb, nq]:
        T [nb, nq, 3, nloc] (stress (xx, yy, xy) of each local Psi function), Gx [nb, nq, ng] (tau_xx of each g term),
        global indices [nb, nloc]"""
        X, Y = self.sx, self.sy
        p = self.ps
        nb, nq = sx.shape
        bx = X.eval_float(np.repeat(ex[:, None], nq, 1), sx, 2)          # [nb, nq, p+1, 3]
        by = Y.eval_float(np.repeat(ey[:, None], nq, 1), sy, 2)
        x = X.af[ex][:, None] + X.hf[ex][:, None] * sx
        y = Y.af[ey][:, None] + Y.hf[ey][:, None] * sy

        def B(dx, dy):
            return np.einsum('bqm,bqn->bqmn', bx[..., dx], by[..., dy]).reshape(nb, nq, -1)
        P0, Px, Py, Pxx, Pxy, Pyy = B(0, 0), B(1, 0), B(0, 1), B(2, 0), B(1, 1), B(0, 2)
        if self.hole:
            e0, e1, e2 = _pv(ETA, x), _pv(ETA, x, 1), _pv(ETA, x, 2)
            z0, z1, z2 = _pv(ZETA, y), _pv(ZETA, y, 1), _pv(ZETA, y, 2)
            r = e0 + z0
            R, Rx, Ry = r * r, 2 * r * e1, 2 * r * z1
            Rxx, Ryy, Rxy = 2 * e1 * e1 + 2 * r * e2, 2 * z1 * z1 + 2 * r * z2, 2 * e1 * z1
        else:
            R = np.ones_like(x)
            Rx = Ry = Rxx = Ryy = Rxy = np.zeros_like(x)
        ex_ = lambda a: a[..., None]
        fxx = ex_(Rxx) * P0 + 2 * ex_(Rx) * Px + ex_(R) * Pxx
        fyy = ex_(Ryy) * P0 + 2 * ex_(Ry) * Py + ex_(R) * Pyy
        fxy = ex_(Rxy) * P0 + ex_(Rx) * Py + ex_(Ry) * Px + ex_(R) * Pxy
        T = np.stack([fyy, fxx, -fxy], axis=2)
        gb = float(self.gbase)
        yy = np.maximum(y - gb, 0.0)
        Gx = np.stack([np.where(y >= gb, yy ** e, 0.0) for e in self.gexp], axis=-1)
        gx = X.first[ex][:, None] + np.arange(p + 1)[None, :]
        gy = Y.first[ey][:, None] + np.arange(p + 1)[None, :]
        gidx = (gx[:, :, None] * Y.n + gy[:, None, :]).reshape(nb, -1)
        return T, Gx, gidx

    def primal_basis(self, ex, ey, sx, sy):
        """strain (xx, yy, gamma) of each local function: Bm [nb, nq, 3, 2 nloc] (first the A part, then B), gidx"""
        X, Y = self.ux, self.uy
        p = self.pu
        nb, nq = sx.shape
        bx = X.eval_float(np.repeat(ex[:, None], nq, 1), sx, 1)
        by = Y.eval_float(np.repeat(ey[:, None], nq, 1), sy, 1)
        Nx = np.einsum('bqm,bqn->bqmn', bx[..., 1], by[..., 0]).reshape(nb, nq, -1)
        Ny = np.einsum('bqm,bqn->bqmn', bx[..., 0], by[..., 1]).reshape(nb, nq, -1)
        Z = np.zeros_like(Nx)
        Bm = np.stack([np.concatenate([Nx, Z], -1), np.concatenate([Z, Ny], -1), np.concatenate([Ny, Nx], -1)], axis=2)
        gx = X.first[ex][:, None] + np.arange(p + 1)[None, :]
        gy = Y.first[ey][:, None] + np.arange(p + 1)[None, :]
        gi = (gx[:, :, None] * Y.n + gy[:, None, :]).reshape(nb, -1)
        return Bm, np.concatenate([gi, gi + X.n * Y.n], axis=1)

    # ---------------------------------------------------------------------------------------------- assembly
    def assemble_dual(self, batch=64):
        """Q(tau) = z^T A z for z = (gamma, C.ravel());  returns A (csr), force vector f (F_q = f.z)"""
        X, Y = self.sx, self.sy
        ng = len(self.gexp)
        nC = X.n * Y.n
        n = ng + nC
        Smat = np.array([[1, -NU, 0], [-NU, 1, 0], [0, 0, 2 * (1 + NU)]])
        rl = self.rules(X, Y, self.nqs)
        rows, cols, vals = [], [], []
        Agg = np.zeros((ng, ng))
        Agc = np.zeros((ng, nC))

        def add(T, Gx, gidx, w):
            # T [nb,nq,3,nl], Gx [nb,nq,ng], w [nb,nq]
            ST = np.einsum('cd,bqdn->bqcn', Smat, T)
            Al = np.einsum('bq,bqcm,bqcn->bmn', w, T, ST)
            nb, nl = gidx.shape
            rows.append(np.repeat(gidx, nl, axis=1).ravel())
            cols.append(np.tile(gidx, (1, nl)).ravel())
            vals.append(Al.ravel())
            # g part: tau_g = (Gx, 0, 0): Q cross terms 2 tau_g S tau = 2 Gx (S00 Txx + S01 Tyy)
            cr = np.einsum('bq,bqk,bqn->bkn', w, Gx, ST[:, :, 0, :])
            for b in range(nb):
                np.add.at(Agc, (slice(None), gidx[b]), cr[b])
            Agg[:] += np.einsum('bq,bqk,bql->kl', w, Gx, Gx)

        full = rl['full']
        for s in range(0, len(full), batch):
            ex, ey = full[s:s + batch, 0], full[s:s + batch, 1]
            nb = len(ex)
            sx = np.broadcast_to(rl['sx'], (nb, len(rl['sx'])))
            sy = np.broadcast_to(rl['sy'], (nb, len(rl['sy'])))
            T, Gx, gidx = self.dual_basis(ex, ey, sx, sy)
            w = rl['wf'][None, :] * (X.hf[ex] * Y.hf[ey])[:, None]
            add(T, Gx, gidx, w)
        for (ex, ey, sx, sy, w) in rl['cut']:
            T, Gx, gidx = self.dual_basis(np.array([ex]), np.array([ey]), sx[None], sy[None])
            add(T, Gx, gidx, w[None])
        Acc = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nC, nC)).tocsr()
        A = sp.bmat([[sp.csr_matrix(Agg), sp.csr_matrix(Agc)], [sp.csr_matrix(Agc.T), Acc]], format='csr')
        f = np.zeros(n)
        f[:ng] = [float((1 - self.gbase) ** (e + 1) / (e + 1)) for e in self.gexp]
        return A, f

    def assemble_primal(self, batch=128):
        """a_q(v,v) = z^T K z + 2 b.z + k0 for z = (A.ravel(), B.ravel()); lift v0 = (x/2, 0)"""
        X, Y = self.ux, self.uy
        n = 2 * X.n * Y.n
        Cm = np.array([[C11, C12, 0], [C12, C11, 0], [0, 0, C66]])
        e0 = np.array([0.5, 0, 0])
        rl = self.rules(X, Y, self.nqu)
        rows, cols, vals = [], [], []
        b = np.zeros(n)
        k0 = 0.0

        def add(Bm, gidx, w):
            nonlocal k0
            CB = np.einsum('cd,bqdn->bqcn', Cm, Bm)
            Kl = np.einsum('bq,bqcm,bqcn->bmn', w, Bm, CB)
            nb, nl = gidx.shape
            rows.append(np.repeat(gidx, nl, axis=1).ravel())
            cols.append(np.tile(gidx, (1, nl)).ravel())
            vals.append(Kl.ravel())
            bl = np.einsum('bq,c,bqcn->bn', w, Cm @ e0, Bm)
            np.add.at(b, gidx.ravel(), bl.ravel())
            k0 += float(np.sum(w)) * float(e0 @ Cm @ e0)

        full = rl['full']
        for s in range(0, len(full), batch):
            ex, ey = full[s:s + batch, 0], full[s:s + batch, 1]
            nb = len(ex)
            sx = np.broadcast_to(rl['sx'], (nb, len(rl['sx'])))
            sy = np.broadcast_to(rl['sy'], (nb, len(rl['sy'])))
            Bm, gidx = self.primal_basis(ex, ey, sx, sy)
            w = rl['wf'][None, :] * (X.hf[ex] * Y.hf[ey])[:, None]
            add(Bm, gidx, w)
        for (ex, ey, sx, sy, w) in rl['cut']:
            Bm, gidx = self.primal_basis(np.array([ex]), np.array([ey]), sx[None], sy[None])
            add(Bm, gidx, w[None])
        K = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
        return K, b, k0

    # ---------------------------------------------------------------------------------------------- patch functional
    def patch_boundary(self, X, Y, nt=24):
        """float rule on the boundary of P_h: points (x, y), weights (ds), outward normals; pieces split at cell lines.
        Returns dict piece -> (x, y, w, nx, ny) for 'outer', 'inner', 'ray'."""
        th0 = self.theta0
        xs = sorted(set(float(v) for v in X.t))
        ys = sorted(set(float(v) for v in Y.t))
        out = {}
        g, w = gauss01(nt)
        for name, Rr in (('outer', float(R2)), ('inner', float(A_HOLE))):
            br = [th0, math.pi / 2]
            for xv in xs:
                if 0 < xv < Rr:
                    t = math.acos(xv / Rr)
                    if th0 < t < math.pi / 2:
                        br.append(t)
            for yv in ys:
                if 0 < yv < Rr:
                    t = math.asin(yv / Rr)
                    if th0 < t < math.pi / 2:
                        br.append(t)
            br = sorted(br)
            P = []
            for t0, t1 in zip(br[:-1], br[1:]):
                th = t0 + (t1 - t0) * g
                sgn = 1.0 if name == 'outer' else -1.0
                P.append((Rr * np.cos(th), Rr * np.sin(th), Rr * (t1 - t0) * w, sgn * np.cos(th), sgn * np.sin(th)))
            out[name] = tuple(np.concatenate([p[k] for p in P]) for k in range(5))
        c0, s0 = math.cos(th0), math.sin(th0)
        r0, r1 = float(A_HOLE), float(R2)
        br = [r0, r1]
        for xv in xs:
            r = xv / c0
            if r0 < r < r1:
                br.append(r)
        for yv in ys:
            r = yv / s0
            if r0 < r < r1:
                br.append(r)
        br = sorted(br)
        P = []
        for a_, b_ in zip(br[:-1], br[1:]):
            r = a_ + (b_ - a_) * g
            P.append((r * c0, r * s0, (b_ - a_) * w, np.full(nt, s0), np.full(nt, -c0)))
        out['ray'] = tuple(np.concatenate([p[k] for p in P]) for k in range(5))
        return out

    def _points_basis(self, X, Y, x, y, nder):
        ex, ey = X.cell_of(x), Y.cell_of(y)
        sx = (x - X.af[ex]) / X.hf[ex]
        sy = (y - Y.af[ey]) / Y.hf[ey]
        bx = X.eval_float(ex, sx, nder)
        by = Y.eval_float(ey, sy, nder)
        gx = X.first[ex][:, None] + np.arange(X.p + 1)[None, :]
        gy = Y.first[ey][:, None] + np.arange(Y.p + 1)[None, :]
        gidx = (gx[:, :, None] * Y.n + gy[:, None, :]).reshape(len(x), -1)
        return bx, by, gidx

    def patch_primal(self):
        """l(v) = l0 + lv.z  for z = (A.ravel(), B.ravel())"""
        X, Y = self.ux, self.uy
        nA = X.n * Y.n
        lv = np.zeros(2 * nA)
        pb = self.patch_boundary(X, Y)
        for name in ('outer', 'inner', 'ray'):
            x, y, w, nx, ny = pb[name]
            bx, by, gidx = self._points_basis(X, Y, x, y, 0)
            N = np.einsum('qm,qn->qmn', bx[..., 0], by[..., 0]).reshape(len(x), -1)
            np.add.at(lv, gidx.ravel(), (C11 * (w * nx)[:, None] * N).ravel())
            np.add.at(lv, (gidx + nA).ravel(), (C12 * (w * ny)[:, None] * N).ravel())
        return C11 * 0.5, lv / float(PH_AREA)

    def patch_dual(self):
        """l(tau) = lt.z for z = (gamma, C.ravel()): (1/|P_h|) oint Phi_y n_y over the outer arc and the ray"""
        X, Y = self.sx, self.sy
        ng = len(self.gexp)
        lt = np.zeros(ng + X.n * Y.n)
        pb = self.patch_boundary(X, Y)
        for name in (('outer', 'ray') if self.hole else ('outer', 'ray', 'inner')):
            x, y, w, nx, ny = pb[name]
            bx, by, gidx = self._points_basis(X, Y, x, y, 1)
            P0 = np.einsum('qm,qn->qmn', bx[..., 0], by[..., 0]).reshape(len(x), -1)
            Py = np.einsum('qm,qn->qmn', bx[..., 0], by[..., 1]).reshape(len(x), -1)
            if self.hole:
                r = _pv(ETA, x) + _pv(ZETA, y)
                Phiy = (2 * r * _pv(ZETA, y, 1))[:, None] * P0 + (r * r)[:, None] * Py
            else:
                Phiy = Py
            np.add.at(lt, (gidx + ng).ravel(), ((w * ny)[:, None] * Phiy).ravel())
            yy = np.maximum(y - float(self.gbase), 0.0)
            for k, e in enumerate(self.gexp):
                lt[k] += np.sum(w * ny * yy ** (e + 1) / (e + 1))
        return lt / float(PH_AREA)

    # ---------------------------------------------------------------------------------------------- ties
    def stress_ties(self):
        """reduction P: z_full = P z_red  (gamma free; C ties of FORMAT.md 3.2)"""
        nx, ny = self.sx.n, self.sy.n
        Ex = np.zeros((nx, nx - 2))
        Ex[0, 0] = Ex[1, 0] = 1
        for i in range(2, nx - 2):
            Ex[i, i - 1] = 1
        Ex[nx - 2, nx - 3] = Ex[nx - 1, nx - 3] = 1
        if self.grip == 'clamped':
            Ex = np.zeros((nx, nx - 1))
            Ex[0, 0] = Ex[1, 0] = 1
            for i in range(2, nx):
                Ex[i, i - 1] = 1
        Ey = np.zeros((ny, ny - 3))
        Ey[0, 0] = Ey[1, 0] = 1
        for j in range(2, ny - 2):
            Ey[j, j - 1] = 1
        Pc = sp.kron(sp.csr_matrix(Ex), sp.csr_matrix(Ey), format='csr')
        ng = len(self.gexp)
        return sp.block_diag([sp.identity(ng, format='csr'), Pc], format='csr')

    def disp_ties(self):
        nx, ny = self.ux.n, self.uy.n
        Ix, Iy = sp.identity(nx, format='csr'), sp.identity(ny, format='csr')
        if self.grip == 'clamped':
            # u_x = 0 on x = 0, u_x = 1/2 and u_y = 0 on x = 1, u_y = 0 on y = 0
            PA = sp.kron(Ix[:, 1:nx - 1], Iy, format='csr')
            PB = sp.kron(Ix[:, :nx - 1], Iy[:, 1:], format='csr')
        else:
            PA = sp.kron(Ix[:, 1:nx - 1], Iy, format='csr')
            PB = sp.kron(Ix, Iy[:, 1:], format='csr')
        return sp.block_diag([PA, PB], format='csr')


def solve_spd(A, rhs, Pr, reg=1e-13, refine=6, protect=0):
    """minimise z^T A z / 2 - rhs^T z over z = Pr y;  returns the full z (several right-hand sides allowed)"""
    Ar = (Pr.T @ A @ Pr).tocsc()
    d = Ar.diagonal()
    scale = d.max()
    keep = d > 1e-14 * scale
    keep[:protect] = d[:protect] > 0                    # never drop the parameters of the particular field
    idx = np.flatnonzero(keep)
    Ak = Ar[idx][:, idx]
    dk = np.sqrt(Ak.diagonal())
    Dm = sp.diags(1 / dk)
    As0 = (Dm @ Ak @ Dm).tocsc()
    lu = spla.splu(As0 + reg * sp.identity(len(idx), format='csc'))
    rhs = np.atleast_2d(np.asarray(rhs).T).T if np.ndim(rhs) == 1 else rhs
    rr = (Pr.T @ rhs)[idx] / dk[:, None]
    ys = lu.solve(rr)
    for _ in range(refine):                              # iterative refinement against the unregularised matrix
        res = rr - As0 @ ys
        ys = ys + lu.solve(res)
    y = np.zeros((Pr.shape[1], rhs.shape[1]))
    y[idx] = ys / dk[:, None]
    return Pr @ y


def uniform_breaks(N, lo=0, hi=1):
    return [Fr(lo) + (Fr(hi) - Fr(lo)) * Fr(k, N) for k in range(N + 1)]
