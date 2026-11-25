"""P3 (block with a spherical cavity, one octant) with higher-degree splines: spaces, ties, float assembly, solves.

Nothing here needs to be rigorous; it only FINDS the fields.  The certified numbers come from p3rig.py.

Octant [0,1]^3 minus the ball of radius 1/2; mirror planes x = 0, y = 0, z = 0; grip x = 1 (frictionless or
clamped); free faces y = 1, z = 1; free sphere.  E = 1, nu = 3/10 (3D).  Full-body stiffness K:
    K <= 8 a_oct(v, v)                       for v = (x/2, 0, 0) + splines meeting the displacement ties,
    K >= 2 F_oct(tau)^2 / Q_oct(tau)          for tau = g(y,z) e_x e_x + inc(rho^2 Psi) meeting the stress ties,
F_oct = int over the grip face of tau_xx, Q_oct = int tau : C^-1 : tau over the octant.
rho = eta(x) + y^2 + z^2 - 1/4, eta(x) = x^2 for x <= 3/4 and x^2 - 4 (x - 3/4)^2 beyond (curved3d's weight, C^1,
eta'(1) = 0); the knot 3/4 has multiplicity ps - 1 in the x-direction of the stress space (Psi is C^1 there).
g(y, z) = sum_k gamma_k (y - 1/2)_+^{e_k} + delta_k (z - 1/2)_+^{e_k}.
"""
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import scipy.linalg as sla

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', '..', 'Reference', 'Curved3D'))
from curved3d import G, COMPS, VAR, build_boxes          # noqa: E402  (Beltrami tensor, index conventions, octree)
from splines import Spline1D, open_knots, gauss01           # noqa: E402

NU = Fr(3, 10)
LAM = NU / ((1 + NU) * (1 - 2 * NU))
MU = 1 / (2 * (1 + NU))
A2 = Fr(1, 4)
B_FLAT = Fr(3, 4)
nu = float(NU)
SCOMP = (1 + nu) * np.diag([1.0, 1, 1, 2, 2, 2])          # compliance in Voigt form (engineering factors)
SCOMP[:3, :3] -= nu
CSTIF = np.diag([2 * float(MU)] * 3 + [float(MU)] * 3)
CSTIF[:3, :3] += float(LAM)
VIDX = {v: i for i, v in enumerate(VAR)}
LCS = np.linalg.cholesky(SCOMP)                           # SCOMP = LCS LCS^T
LCC = np.linalg.cholesky(CSTIF)


def eta_poly(e, N):
    """eta on cell e (physical x), ascending coefficients (Fractions), valid for x in [e/N, (e+1)/N]"""
    if Fr(e, N) < B_FLAT:
        return [Fr(0), Fr(0), Fr(1)]
    g = 1 / (1 - B_FLAT)                                    # 4
    b = B_FLAT
    return [-g * b * b, 2 * g * b, 1 - g]                   # x^2 - g (x - b)^2


def pval(c, x, d=0):
    c = [float(v) for v in c]
    for _ in range(d):
        c = [k * c[k] for k in range(1, len(c))] or [0.0]
    r = np.zeros_like(x, dtype=float)
    for v in reversed(c):
        r = r * x + v
    return r


class Setup3:
    def __init__(self, N, ps, pu, grip='frictionless', gexp=None, sym=True):
        self.sym = sym                  # restrict to fields symmetric under y <-> z (the problem is)
        assert N % 4 == 0, "N must be a multiple of 4 (knots at 1/2 and 3/4)"
        self.N, self.ps, self.pu, self.grip = N, ps, pu, grip
        br = [Fr(k, N) for k in range(N + 1)]
        self.S = [Spline1D(open_knots(br, ps, {B_FLAT: ps - 1}), ps), Spline1D(open_knots(br, ps), ps),
                  Spline1D(open_knots(br, ps), ps)]
        self.U = [Spline1D(open_knots(br, pu), pu) for _ in range(3)]
        self.gexp = tuple(gexp if gexp is not None else (ps, ps + 1, ps + 2))
        self.nqs = ps + 5          # tau degree <= ps + 4 per variable (rho^2 has degree 4), g degree <= max(gexp)
        self.nqs = max(self.nqs, (2 * max(self.gexp)) // 2 + 1)
        self.nqu = pu + 1
        self.ns = [S.n for S in self.S]
        self.nu_ = [U.n for U in self.U]

    # ------------------------------------------------------------------ ties (reduction maps full -> reduced index)
    def stress_map(self):
        """red[i] = reduced index of full Psi coefficient i (or -1 if it is zero); full layout [6, nx, ny, nz]"""
        nx, ny, nz = self.ns
        red = -np.ones((6, nx, ny, nz), dtype=np.int64)
        cnt = 0
        for c, (a, b) in enumerate(COMPS):
            # per direction: index map 1D -> reduced 1D (-1: zero), sharing via "owner"
            maps = []
            for k, n in enumerate((nx, ny, nz)):
                s = (-1) ** ((a == k) + (b == k))
                m = np.arange(n)
                zero = np.zeros(n, dtype=bool)
                # mirror plane at 0: even -> c0 = c1, odd -> c0 = 0
                if s == 1:
                    m[0] = 1
                else:
                    zero[0] = True
                if k == 0:
                    if self.grip == 'frictionless':                   # grip x = 1 acts like a mirror plane
                        if s == 1:
                            m[n - 1] = n - 2
                        else:
                            zero[n - 1] = True
                else:                                                 # free faces y = 1, z = 1
                    zero[n - 1] = zero[n - 2] = True
                m = np.where(zero, -1, m)
                maps.append(m)
            mx, my, mz = maps
            own = {}
            for ix in range(nx):
                if mx[ix] < 0:
                    continue
                for iy in range(ny):
                    if my[iy] < 0:
                        continue
                    for iz in range(nz):
                        if mz[iz] < 0:
                            continue
                        key = (mx[ix], my[iy], mz[iz])
                        if key not in own:
                            own[key] = cnt
                            cnt += 1
                        red[c, ix, iy, iz] = own[key]
        if self.sym:
            # y <-> z: Psi_c(x, y, z) = Psi_swap(c)(x, z, y) with xx->xx, yy->zz, zz->yy, xy->xz, yz->yz, xz->xy
            swap = [0, 2, 1, 5, 4, 3]
            partner = red[swap].transpose(0, 1, 3, 2)
            return _merge(red.ravel(), partner.ravel())
        return red.ravel(), cnt

    def disp_map(self):
        nx, ny, nz = self.nu_
        red = -np.ones((3, nx, ny, nz), dtype=np.int64)
        cnt = 0
        for a in range(3):
            ok = np.ones((nx, ny, nz), dtype=bool)
            if a == 0:
                ok[0, :, :] = False          # u_x = 0 on x = 0
                ok[-1, :, :] = False         # u_x = 1/2 on x = 1 (lift)
            if a == 1:
                ok[:, 0, :] = False          # u_y = 0 on y = 0
            if a == 2:
                ok[:, :, 0] = False          # u_z = 0 on z = 0
            if self.grip == 'clamped' and a > 0:
                ok[-1, :, :] = False
            idx = np.flatnonzero(ok.ravel())
            r = red[a].reshape(-1)
            r[idx] = cnt + np.arange(len(idx))
            cnt += len(idx)
        if self.sym:
            partner = red[[0, 2, 1]].transpose(0, 1, 3, 2)          # v_x(x,y,z) = v_x(x,z,y), v_y(x,y,z) = v_z(x,z,y)
            return _merge(red.ravel(), partner.ravel())
        return red.ravel(), cnt

    def param_map(self):
        """reduced index of (gamma_k, delta_k); with the symmetry delta_k = gamma_k"""
        ng = len(self.gexp)
        if self.sym:
            return np.concatenate([np.arange(ng), np.arange(ng)]), ng
        return np.arange(2 * ng), 2 * ng

    def solve_boxes(self, L, centre=False):
        """octree boxes for the float solve: the superset of curved3d (centre=False), or without the cut boxes whose
        centre lies in the ball (centre=True, close to the exact domain; used when the certificate is exact)"""
        b = build_boxes(self.N, L, A2)
        if not centre:
            return b
        n = self.N * 2 ** b[:, 3]
        lo = b[:, :3] * (2 ** b[:, 3])[:, None] + b[:, 4:7]
        c2 = ((2 * lo + 1) ** 2).sum(axis=1)                     # (2 lo + 1)^2 / (2 n)^2 = |centre|^2
        keep = (b[:, 7] == 0) | (c2 * 4 >= (2 * n) ** 2)       # |centre|^2 >= 1/4
        return b[keep]

    # ------------------------------------------------------------------ tables
    def _tab(self, Sp, e, lev, j, nq, nder, cache):
        key = (id(Sp), e, lev, j, nq)
        if key not in cache:
            g, w = gauss01(nq)
            s = (j + g) / 2 ** lev
            cache[key] = (Sp.eval_float(e, s, nder), w / 2 ** lev * Sp.hf[e], Sp.af[e] + Sp.hf[e] * s)
        return cache[key]

    def dual_box(self, box, cache):
        """T [nQ, 6, 6*(p+1)^3] stress of each local Psi function, Tg [nQ, 2 ng] (tau_xx of the g terms), W [nQ],
        cell (ex, ey, ez)"""
        ex, ey, ez, lev, jx, jy, jz = [int(v) for v in box[:7]]
        p, nq = self.ps, self.nqs
        bx, wx, x = self._tab(self.S[0], ex, lev, jx, nq, 2, cache)
        by, wy, y = self._tab(self.S[1], ey, lev, jy, nq, 2, cache)
        bz, wz, z = self._tab(self.S[2], ez, lev, jz, nq, 2, cache)
        nl = (p + 1) ** 3

        def B(d):
            return np.einsum('im,jn,ko->ijkmno', bx[:, :, d[0]], by[:, :, d[1]], bz[:, :, d[2]]).reshape(nq ** 3, nl)
        Bv = {v: B(v) for v in VAR}
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()
        ep = eta_poly(ex, self.N)
        r = pval(ep, X) + Y * Y + Z * Z - 0.25
        d1 = [pval(ep, X, 1), 2 * Y, 2 * Z]
        d2 = [pval(ep, X, 2), 2.0 + 0 * Y, 2.0 + 0 * Z]
        R = r * r
        T = np.zeros((nq ** 3, 6, 6, nl))                      # [Q, s, component, m]
        H = []
        for (k, m) in COMPS:
            dk = tuple(1 if i == k else 0 for i in range(3))
            dm = tuple(1 if i == m else 0 for i in range(3))
            dkm = tuple(dk[i] + dm[i] for i in range(3))
            h = R[:, None] * Bv[dkm] + (2 * r * d1[k])[:, None] * Bv[dm] + (2 * r * d1[m])[:, None] * Bv[dk]
            hh = 2 * d1[k] * d1[m] + (2 * r * d2[k] if k == m else 0.0)
            h = h + hh[:, None] * Bv[(0, 0, 0)]
            H.append(h)
        H = np.stack(H, axis=1)                                 # [Q, c, m]
        T = np.einsum('spc,Qcm->Qspm', G, H).reshape(nq ** 3, 6, 6 * nl)
        ng = len(self.gexp)
        Tg = np.zeros((nq ** 3, 2 * ng))
        for i, e in enumerate(self.gexp):
            Tg[:, i] = np.where(Y > 0.5, np.maximum(Y - 0.5, 0) ** e, 0.0)
            Tg[:, ng + i] = np.where(Z > 0.5, np.maximum(Z - 0.5, 0) ** e, 0.0)
        W = (wx[:, None, None] * wy[None, :, None] * wz[None, None, :]).ravel()
        return T, Tg, W, (ex, ey, ez)

    def stress_gidx(self, cell):
        ex, ey, ez = cell
        nx, ny, nz = self.ns
        p = self.ps
        ix = self.S[0].first[ex] + np.arange(p + 1)
        iy = self.S[1].first[ey] + np.arange(p + 1)
        iz = self.S[2].first[ez] + np.arange(p + 1)
        loc = ((ix[:, None, None] * ny + iy[None, :, None]) * nz + iz[None, None, :]).ravel()
        return (np.arange(6)[:, None] * (nx * ny * nz) + loc[None, :]).ravel()

    def assemble_dual(self, L, verbose=False, centre=False):
        """dense reduced matrix A (Q = z.A.z), force vector f (F_oct = f.z), z = (gamma, delta, reduced Psi)"""
        red, nr = self.stress_map()
        pidx, npar = self.param_map()
        ng2 = npar
        self._pidx = pidx
        n = npar + nr
        A = np.zeros((n, n))
        f = np.zeros(n)
        ng = len(self.gexp)
        for i, e in enumerate(self.gexp):
            f[pidx[i]] += 0.5 ** (e + 1) / (e + 1)
            f[pidx[ng + i]] += 0.5 ** (e + 1) / (e + 1)
        boxes = self.solve_boxes(L, centre)
        order = np.lexsort((boxes[:, 2], boxes[:, 1], boxes[:, 0]))
        boxes = boxes[order]
        cache = {}
        t0 = time.time()
        cur, acc = None, None
        for bi, box in enumerate(boxes):
            T, Tg, W, cell = self.dual_box(box, cache)
            Tall = np.concatenate([Tg[:, None, :] * np.array([1.0, 0, 0, 0, 0, 0])[None, :, None], T], axis=2)
            Tp = np.einsum('ts,Qsn->Qtn', LCS.T, Tall * np.sqrt(W)[:, None, None]).reshape(-1, Tall.shape[2])
            Al = Tp.T @ Tp
            if cell != cur:
                if cur is not None:
                    self._scatter(A, acc, cur, red, ng2)
                cur, acc = cell, Al
            else:
                acc = acc + Al
            if verbose and bi % 200 == 0:
                print(f"    dual assembly {bi}/{len(boxes)} boxes  {time.time() - t0:.0f}s", flush=True)
        self._scatter(A, acc, cur, red, ng2)
        return A, f, red, nr

    def _scatter(self, A, Al, cell, red, ng2):
        g = red[self.stress_gidx(cell)]
        idx = np.concatenate([self._pidx, np.where(g >= 0, g + ng2, -1)])
        keep = idx >= 0
        ii = idx[keep]
        sub = Al[np.ix_(keep, keep)]
        np.add.at(A, (ii[:, None], ii[None, :]), sub)

    # ------------------------------------------------------------------ displacement
    def primal_box(self, box, cache):
        ex, ey, ez, lev, jx, jy, jz = [int(v) for v in box[:7]]
        p, nq = self.pu, self.nqu
        bx, wx, _ = self._tab(self.U[0], ex, lev, jx, nq, 1, cache)
        by, wy, _ = self._tab(self.U[1], ey, lev, jy, nq, 1, cache)
        bz, wz, _ = self._tab(self.U[2], ez, lev, jz, nq, 1, cache)
        nl = (p + 1) ** 3
        D = []
        for k in range(3):
            d = [0, 0, 0]
            d[k] = 1
            D.append(np.einsum('im,jn,ko->ijkmno', bx[:, :, d[0]], by[:, :, d[1]], bz[:, :, d[2]]).reshape(nq ** 3, nl))
        # strain (xx, yy, zz, xy, yz, xz) engineering, of v_a = B e_a
        E = np.zeros((nq ** 3, 6, 3, nl))
        E[:, 0, 0] = D[0]
        E[:, 1, 1] = D[1]
        E[:, 2, 2] = D[2]
        E[:, 3, 0], E[:, 3, 1] = D[1], D[0]
        E[:, 4, 1], E[:, 4, 2] = D[2], D[1]
        E[:, 5, 0], E[:, 5, 2] = D[2], D[0]
        W = (wx[:, None, None] * wy[None, :, None] * wz[None, None, :]).ravel()
        return E.reshape(nq ** 3, 6, 3 * nl), W, (ex, ey, ez)

    def disp_gidx(self, cell):
        ex, ey, ez = cell
        nx, ny, nz = self.nu_
        p = self.pu
        ix = self.U[0].first[ex] + np.arange(p + 1)
        iy = self.U[1].first[ey] + np.arange(p + 1)
        iz = self.U[2].first[ez] + np.arange(p + 1)
        loc = ((ix[:, None, None] * ny + iy[None, :, None]) * nz + iz[None, None, :]).ravel()
        return (np.arange(3)[:, None] * (nx * ny * nz) + loc[None, :]).ravel()

    def assemble_primal(self, L, centre=False):
        """a_oct(v,v) = z.K.z + 2 b.z + k0 for the reduced coefficients z"""
        red, nr = self.disp_map()
        K = np.zeros((nr, nr))
        b = np.zeros(nr)
        k0 = 0.0
        e0 = np.array([0.5, 0, 0, 0, 0, 0])
        boxes = self.solve_boxes(L, centre)
        boxes = boxes[np.lexsort((boxes[:, 2], boxes[:, 1], boxes[:, 0]))]
        cache = {}
        cur, acc, accb = None, None, None
        for box in boxes:
            E, W, cell = self.primal_box(box, cache)
            Ep = np.einsum('ts,Qsn->Qtn', LCC.T, E * np.sqrt(W)[:, None, None]).reshape(-1, E.shape[2])
            Kl = Ep.T @ Ep
            bl = np.einsum('Q,s,Qsn->n', W, CSTIF @ e0, E)
            k0 += float(W.sum()) * float(e0 @ CSTIF @ e0)
            if cell != cur:
                if cur is not None:
                    self._scatter_u(K, b, acc, accb, cur, red)
                cur, acc, accb = cell, Kl, bl
            else:
                acc, accb = acc + Kl, accb + bl
        self._scatter_u(K, b, acc, accb, cur, red)
        return K, b, k0, red, nr

    def _scatter_u(self, K, b, Kl, bl, cell, red):
        g = red[self.disp_gidx(cell)]
        keep = g >= 0
        ii = g[keep]
        np.add.at(K, (ii[:, None], ii[None, :]), Kl[np.ix_(keep, keep)])
        np.add.at(b, ii, bl[keep])


def _merge(red, partner):
    """identify reduced indices red[i] and partner[i] (union-find); returns the renumbered map and the count"""
    n = red.max() + 1
    parent = np.arange(n)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a, b in zip(red, partner):
        if a >= 0 and b >= 0:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
        elif (a >= 0) != (b >= 0):
            raise ValueError("symmetry partner of a free coefficient is fixed to zero")
    roots = np.array([find(a) for a in range(n)])
    uniq, newid = np.unique(roots, return_inverse=True)
    out = np.where(red >= 0, newid[np.maximum(red, 0)], -1)
    return out, len(uniq)


def solve_dense(A, rhs, reg=1e-13, protect=0):
    """minimise z.A.z/2 - rhs.z (several right-hand sides); rows with zero diagonal are dropped"""
    d = np.diag(A).copy()
    scale = d.max()
    keep = d > 1e-14 * scale
    keep[:protect] = d[:protect] > 0
    idx = np.flatnonzero(keep)
    dk = np.sqrt(d[idx])
    As = A[np.ix_(idx, idx)] / dk[:, None] / dk[None, :]
    As[np.diag_indices_from(As)] += reg
    rr = np.atleast_2d(rhs.T).T[idx] / dk[:, None]
    lu = sla.lu_factor(As, overwrite_a=False, check_finite=False)
    ys = sla.lu_solve(lu, rr)
    As[np.diag_indices_from(As)] -= reg
    for _ in range(3):
        ys = ys + sla.lu_solve(lu, rr - As @ ys)
    z = np.zeros((A.shape[0], rr.shape[1]))
    z[idx] = ys / dk[:, None]
    return z
