"""Near-optimal admissible fields from a small floating-point Ritz solver (test use: sharper falsification test, and
coarse bounds to set against the builder's intervals).  Written for the checker; shares no code with the builder.

Only the coefficients come from floating point.  They are stored as float64 and the checker then treats them as exact
rationals: every bound computed from them by check.py is rigorous whatever the quality of this solver.  Quadrature here
is approximate on the cells cut by the hole (sub-cells, points inside the disc dropped) and on the patch (polar Gauss).
"""
from fractions import Fraction as Fr

import numpy as np

import exact as X
import fields as FD
import synth

D = FD.D
Cm = np.array(D.C, dtype=float)
Sm = np.array(D.S, dtype=float)


def gauss(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return (x + 1) / 2, w / 2


class Basis1D:
    """float local power coefficients of the B-spline pieces (from the exact pieces)"""

    def __init__(self, knots, p):
        self.p = p
        self.knots = np.array([float(v) for v in knots])
        self.spans = [(k, float(t0), float(t1), np.array(M, dtype=float)) for (k, t0, t1, M) in X.bspline_pieces(knots, p)]
        self.n = len(knots) - p - 1

    def find(self, x):
        """span index (into self.spans) containing x"""
        lefts = np.array([s[1] for s in self.spans])
        return np.clip(np.searchsorted(lefts, x, side='right') - 1, 0, len(self.spans) - 1)

    def eval(self, si, s):
        """values, 1st and 2nd x-derivatives of the p+1 active functions at local points s (array) of span si"""
        k, t0, t1, M = self.spans[si]
        h = t1 - t0
        p = self.p
        P0 = np.vander(s, p + 1, increasing=True)                    # s^a
        a = np.arange(p + 1)
        P1 = np.zeros_like(P0)
        P1[:, 1:] = P0[:, :-1] * a[1:]
        P2 = np.zeros_like(P0)
        P2[:, 2:] = P0[:, :-2] * (a[2:] * (a[2:] - 1))
        return P0 @ M.T, P1 @ M.T / h, P2 @ M.T / (h * h), k - p


def cell_points(x0, x1, y0, y1, ng, hole, depth=4):
    """quadrature points (x, y, w) of the cell (minus the disc of radius 1/2 if hole); approximate on cut sub-cells"""
    gx, gw = gauss(ng)
    out = []

    def rec(a0, a1, b0, b1, d):
        if hole and a1 * a1 + b1 * b1 <= 0.25:
            return
        full = (not hole) or a0 * a0 + b0 * b0 >= 0.25
        if full or d == 0:
            X_, Y_ = np.meshgrid(a0 + (a1 - a0) * gx, b0 + (b1 - b0) * gx, indexing='ij')
            W = np.outer(gw, gw) * (a1 - a0) * (b1 - b0)
            m = np.ones_like(W, dtype=bool) if full else (X_ ** 2 + Y_ ** 2 >= 0.25)
            out.append((X_[m], Y_[m], W[m]))
            return
        am, bm = (a0 + a1) / 2, (b0 + b1) / 2
        for c in ((a0, am, b0, bm), (am, a1, b0, bm), (a0, am, bm, b1), (am, a1, bm, b1)):
            rec(*c, d - 1)
    rec(x0, x1, y0, y1, depth)
    if not out:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    return tuple(np.concatenate([o[i] for o in out]) for i in range(3))


def patch_points(name, n=48):
    if name == 'P1':
        gx, gw = gauss(n)
        r = 0.5 + 0.05 * gx
        th = np.pi / 2 - 0.1 + 0.1 * gx
        R, T = np.meshgrid(r, th, indexing='ij')
        W = np.outer(gw * 0.05, gw * 0.1) * R
        return (R * np.cos(T)).ravel(), (R * np.sin(T)).ravel(), W.ravel(), float(D.P1_QUARTER['half_patch_area'])
    # P2: [0, 1/10]^2, split at 1/20 (a knot or not, the error only affects the quality of the fields)
    gx, gw = gauss(n)
    xs = np.concatenate([0.05 * gx, 0.05 + 0.05 * gx])
    ws = np.concatenate([0.05 * gw, 0.05 * gw])
    Xp, Yp = np.meshgrid(xs, xs, indexing='ij')
    return Xp.ravel(), Yp.ravel(), np.outer(ws, ws).ravel(), float(D.P2['patch_area'])


def points_by_cell(bx, by, X_, Y_, W):
    ix, iy = bx.find(X_), by.find(Y_)
    groups = {}
    for idx in range(len(X_)):
        groups.setdefault((int(ix[idx]), int(iy[idx])), []).append(idx)
    for key, ids in groups.items():
        ids = np.array(ids)
        i, j = key
        _, x0, x1, _ = bx.spans[i]
        _, y0, y1, _ = by.spans[j]
        yield i, j, (X_[ids] - x0) / (x1 - x0), (Y_[ids] - y0) / (y1 - y0), X_[ids], Y_[ids], W[ids]


# ------------------------------------------------------------------------------------------------- primal
def disp_G(bx, by, i, j, s, t, nyy):
    """strain rows (npts, 3, ndof) of the local functions, dof index = 2*(I*ny + J) + comp"""
    Nx, dNx, _, ox = bx.eval(i, s)
    Ny, dNy, _, oy = by.eval(j, t)
    p1 = Nx.shape[1]
    q1 = Ny.shape[1]
    npts = len(s)
    G = np.zeros((npts, 3, p1 * q1 * 2))
    idx = np.zeros(p1 * q1 * 2, dtype=int)
    for a in range(p1):
        for b in range(q1):
            m = 2 * (a * q1 + b)
            gI = (ox + a) * nyy + (oy + b)
            idx[m], idx[m + 1] = 2 * gI, 2 * gI + 1
            G[:, 0, m] = dNx[:, a] * Ny[:, b]          # exx from A
            G[:, 2, m] = Nx[:, a] * dNy[:, b]          # gxy from A
            G[:, 1, m + 1] = Nx[:, a] * dNy[:, b]      # eyy from B
            G[:, 2, m + 1] = dNx[:, a] * Ny[:, b]      # gxy from B
    return G, idx


def solve_primal(name, variant, kx, ky, p, kappa):
    bx, by = Basis1D(kx, p), Basis1D(ky, p)
    nx, ny = bx.n, by.n
    ndof = 2 * nx * ny
    K = np.zeros((ndof, ndof))
    bvec = np.zeros(ndof)
    E0 = np.array([0.5, 0, 0]) if name == 'P1' else np.array([0, 0, -0.5])
    hole = name == 'P1'
    for i, (_, x0, x1, _) in enumerate(bx.spans):
        for j, (_, y0, y1, _) in enumerate(by.spans):
            if name == 'P2' and x0 >= 0 and y1 <= 0:
                continue
            Xq, Yq, Wq = cell_points(x0, x1, y0, y1, p + 2, hole)
            if len(Wq) == 0:
                continue
            G, idx = disp_G(bx, by, i, j, (Xq - x0) / (x1 - x0), (Yq - y0) / (y1 - y0), ny)
            CG = np.einsum('kl,qlm->qkm', Cm, G)
            K[np.ix_(idx, idx)] += np.einsum('q,qkm,qkn->mn', Wq, G, CG)
            bvec[idx] += np.einsum('q,qkm,k->m', Wq, CG, E0)
    L = np.zeros(ndof)
    Xp, Yp, Wp, meas = patch_points(name)
    for i, j, s, t, _, _, w in points_by_cell(bx, by, Xp, Yp, Wp):
        G, idx = disp_G(bx, by, i, j, s, t, ny)
        L[idx] += np.einsum('q,qkm,k->m', w, G, Cm[0]) / meas
    # fixed dofs (zero)
    fixed = np.zeros(ndof, dtype=bool)
    for I in range(nx):
        for J in range(ny):
            gA, gB = 2 * (I * ny + J), 2 * (I * ny + J) + 1
            if name == 'P1':
                if I in (0, nx - 1):
                    fixed[gA] = True
                if J == 0 or (variant == 'clamped' and I == nx - 1):
                    fixed[gB] = True
            else:
                if I in (0, nx - 1):
                    fixed[gA] = fixed[gB] = True
    # dofs with no support in the domain (P2: functions inside the removed square) are fixed too
    fixed |= np.abs(np.diag(K)) == 0
    free = ~fixed
    sol = {}
    for key, k in (('v', 0.0), ('vp', kappa), ('vm', -kappa)):
        c = np.zeros(ndof)
        c[free] = np.linalg.solve(K[np.ix_(free, free)], (k * L - bvec)[free])
        sol[key] = (c[0::2].reshape(nx, ny), c[1::2].reshape(nx, ny))
    return sol


# ------------------------------------------------------------------------------------------------- dual
def rho2_derivs(x, y):
    """rho^2 and its derivatives up to order 2 (float), rho from definitions.py"""
    rho = sum(float(c) * x ** i * y ** j for (i, j), c in D.RHO.items())
    rx = sum(float(c) * i * x ** (i - 1) * y ** j for (i, j), c in D.RHO.items() if i >= 1)
    ry = sum(float(c) * j * x ** i * y ** (j - 1) for (i, j), c in D.RHO.items() if j >= 1)
    rxx = sum(float(c) * i * (i - 1) * x ** (i - 2) * y ** j for (i, j), c in D.RHO.items() if i >= 2)
    ryy = sum(float(c) * j * (j - 1) * x ** i * y ** (j - 2) for (i, j), c in D.RHO.items() if j >= 2)
    rxy = sum(float(c) * i * j * x ** (i - 1) * y ** (j - 1) for (i, j), c in D.RHO.items() if i >= 1 and j >= 1)
    R = rho * rho
    return R, 2 * rho * rx, 2 * rho * ry, 2 * rx * rx + 2 * rho * rxx, 2 * rx * ry + 2 * rho * rxy, 2 * ry * ry + 2 * rho * ryy


def stress_T(bx, by, i, j, s, t, x, y, ny, weight):
    """stress rows (npts, 3, nloc) of the local Psi functions (Airy of rho^2 N or N), dof = I*ny + J"""
    Nx, dNx, ddNx, ox = bx.eval(i, s)
    Ny, dNy, ddNy, oy = by.eval(j, t)
    if weight:
        R, Rx, Ry, Rxx, Rxy, Ryy = rho2_derivs(x, y)
    else:
        R, Rx, Ry, Rxx, Rxy, Ryy = (np.ones_like(x),) + (np.zeros_like(x),) * 5
    p1, q1 = Nx.shape[1], Ny.shape[1]
    T = np.zeros((len(s), 3, p1 * q1))
    idx = np.zeros(p1 * q1, dtype=int)
    for a in range(p1):
        for b in range(q1):
            m = a * q1 + b
            idx[m] = (ox + a) * ny + (oy + b)
            T[:, 0, m] = Nx[:, a] * (Ryy * Ny[:, b] + 2 * Ry * dNy[:, b] + R * ddNy[:, b])
            T[:, 1, m] = Ny[:, b] * (Rxx * Nx[:, a] + 2 * Rx * dNx[:, a] + R * ddNx[:, a])
            T[:, 2, m] = -(Rxy * Nx[:, a] * Ny[:, b] + Rx * Nx[:, a] * dNy[:, b] + Ry * dNx[:, a] * Ny[:, b]
                           + R * dNx[:, a] * dNy[:, b])
    return T, idx


def particular_T(name, x, y, gexp, gbase):
    """stress rows (npts, 3, npar) of the particular parameters: P1 the g terms, P2 (alpha, beta, gamma)"""
    if name == 'P1':
        T = np.zeros((len(x), 3, len(gexp)))
        on = y >= gbase
        for k, e in enumerate(gexp):
            T[on, 0, k] = (y[on] - gbase) ** e
        return T
    # phi_p = chi(y)(alpha + beta x + gamma y)
    T = np.zeros((len(x), 3, 3))
    up = y > 0
    chi = np.where(up, 1 - 3 * y ** 2 + 2 * y ** 3, 1.0)
    d1 = np.where(up, -6 * y + 6 * y ** 2, 0.0)
    d2 = np.where(up, -6 + 12 * y, 0.0)
    lin = [np.ones_like(x), x, y]
    dlin_x = [0 * x, 1 + 0 * x, 0 * x]
    dlin_y = [0 * x, 0 * x, 1 + 0 * x]
    for k in range(3):
        T[:, 0, k] = d2 * lin[k] + 2 * d1 * dlin_y[k]                         # phi_yy
        T[:, 1, k] = 0                                                         # phi_xx (lin is linear)
        T[:, 2, k] = -(d1 * dlin_x[k] + chi * 0)                               # -phi_xy
    return T


def solve_dual(name, variant, kx, ky, p, kappa, gexp=(), gbase=None):
    bx, by = Basis1D(kx, p), Basis1D(ky, p)
    nx, ny = bx.n, by.n
    npar = len(gexp) if name == 'P1' else 3
    nd = nx * ny + npar
    M = np.zeros((nd, nd))
    hole = name == 'P1'
    weight = name == 'P1'
    for i, (_, x0, x1, _) in enumerate(bx.spans):
        for j, (_, y0, y1, _) in enumerate(by.spans):
            if name == 'P2' and x0 >= 0 and y1 <= 0:
                continue
            Xq, Yq, Wq = cell_points(x0, x1, y0, y1, p + (10 if weight else 4), hole)
            if len(Wq) == 0:
                continue
            T, idx = stress_T(bx, by, i, j, (Xq - x0) / (x1 - x0), (Yq - y0) / (y1 - y0), Xq, Yq, ny, weight)
            Tp = particular_T(name, Xq, Yq, gexp, gbase)
            Tall = np.concatenate([T, Tp], axis=2)
            ia = np.concatenate([idx, nx * ny + np.arange(npar)])
            ST = np.einsum('kl,qlm->qkm', Sm, Tall)
            M[np.ix_(ia, ia)] += np.einsum('q,qkm,qkn->mn', Wq, Tall, ST)
    Lt = np.zeros(nd)
    Xp, Yp, Wp, meas = patch_points(name)
    for i, j, s, t, x, y, w in points_by_cell(bx, by, Xp, Yp, Wp):
        T, idx = stress_T(bx, by, i, j, s, t, x, y, ny, weight)
        Lt[idx] += (T[:, 0, :] * w[:, None]).sum(0) / meas
        Tp = particular_T(name, x, y, gexp, gbase)
        Lt[nx * ny:] += (Tp[:, 0, :] * w[:, None]).sum(0) / meas
    f = np.zeros(nd)
    if name == 'P1':
        for k, e in enumerate(gexp):
            f[nx * ny + k] = (1 - gbase) ** (e + 1) / (e + 1)
        wf = 0.5
    else:
        f[nx * ny + 1] = -1.0
        wf = 1.0
    # ties of FORMAT.md: free parameters
    rep = {}
    for I in range(nx):
        for J in range(ny):
            if name == 'P1':
                if J >= ny - 2:
                    continue
                I2 = max(I, 1)
                if variant != 'clamped':
                    I2 = min(I2, nx - 2)
                rep[(I, J)] = (I2, max(J, 1))
            else:
                rep[(I, J)] = (I, J)
    if name == 'P2':
        Lx, Rx, Sx = synth.p2_classes(kx, p)
        Ly, Ry, Sy = synth.p2_classes(ky, p)
        for (I, J) in list(rep):
            if (J in (ny - 2, ny - 1) or ((I in Lx or I in Sx) and J in (0, 1)) or (I in Sx and (J in Ly or J in Sy))
                    or (I in Rx and (J in Sy or J in Ly))):
                del rep[(I, J)]
    reps = sorted(set(rep.values()))
    col = {r: c for c, r in enumerate(reps)}
    Tm = np.zeros((nd, len(reps) + npar))
    for (I, J), r in rep.items():
        Tm[I * ny + J, col[r]] = 1.0
    for k in range(npar):
        Tm[nx * ny + k, len(reps) + k] = 1.0
    Mr = Tm.T @ M @ Tm
    Mr += 1e-13 * np.trace(Mr) / len(Mr) * np.eye(len(Mr))
    fr, Lr = Tm.T @ f, Tm.T @ Lt
    sol = {}
    for key, rhs in (('tau', fr), ('tp', wf * fr - kappa * Lr), ('tm', wf * fr + kappa * Lr)):
        d = Tm @ np.linalg.solve(Mr, rhs)
        Cc = d[:nx * ny].reshape(nx, ny)
        par = d[nx * ny:]
        sol[key] = (Cc, par)
    return sol


def certificate(name, variant, breaks_x, breaks_y, pu, ps, kappa, gexp=(0, 1, 2, 3), gbase=Fr(1, 2), mult0=True):
    """(meta, arrays) of a certificate made of near-optimal fields (arrays: exact Fractions of the float64 values)"""
    kxu, kyu = synth.open_knots(pu, breaks_x), synth.open_knots(pu, breaks_y)
    m = {Fr(0): ps - 1} if (name == 'P2' and mult0) else None
    kxs, kys = synth.open_knots(ps, breaks_x, m), synth.open_knots(ps, breaks_y, m)
    pr = solve_primal(name, variant, kxu, kyu, pu, float(kappa))
    du = solve_dual(name, variant, kxs, kys, ps, float(kappa), gexp, float(gbase))
    sp = dict(g_exponents=list(gexp), g_base=str(gbase)) if name == 'P1' else {}
    meta = dict(format='benchpack-certificate-1', problem=name, variant=variant,
                model='quarter' if name == 'P1' else 'full',
                displacement_space=synth.space(pu, kxu, kyu), stress_space=synth.space(ps, kxs, kys, **sp),
                fields=dict(Q1=dict(v='v', tau='tau'),
                            Q2=dict(kappa=str(kappa), v_plus='vp', v_minus='vm', tau_plus='tp', tau_minus='tm')),
                label='test')
    arrays = {}
    for k, (A, B) in pr.items():
        arrays[k + '_A'], arrays[k + '_B'] = X.farr(A), X.farr(B)
    for k, (Cc, par) in du.items():
        arrays[k + '_C'], arrays[k + '_par'] = X.farr(Cc), X.farr(par)
    return meta, arrays
