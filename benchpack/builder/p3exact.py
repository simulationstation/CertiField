"""P3: exact treatment of the cells cut by the sphere (integration over the cell minus the ball).

On a cut cell the energy density of a stored field is built as one polynomial with rational coefficients
(python-flint, 3 variables) and integrated over the part of the cell outside the ball with p3cut.outside_integral
(exact rationals, Arb balls for the curved part).  Full cells are left to p3rig (interval Gauss rules).
Results per cell are stored in a checkpoint file (rational upper bounds) so that a stopped run resumes.
"""
import json
import os
import time
from fractions import Fraction as Fr

import numpy as np
from flint import fmpq, ctx as fctx

import p3float as pf
from p3float import G, COMPS, A2, eta_poly
from p3cut import X3, Y3, Z3, ZERO3, uni3, local3, outside_integral
from exact import q, arb_bounds


class Pieces3:
    def __init__(self, Sp, var):
        self.Sp, self.var, self.c = Sp, var, {}

    def __call__(self, e):
        if e not in self.c:
            Sp = self.Sp
            self.c[e] = [local3(Sp.pieces[e][m], Sp.a[e], Sp.h[e], self.var) for m in range(Sp.p + 1)]
        return self.c[e]


def tensor_cell(px, py, pz, coef, fx, fy, fz):
    """sum coef[fx+i, fy+j, fz+k] px[i] py[j] pz[k] as one polynomial"""
    p1 = len(px)
    tot = ZERO3
    for i in range(p1):
        acc_i = ZERO3
        for j in range(p1):
            acc_j = ZERO3
            for k in range(p1):
                c = float(coef[fx + i, fy + j, fz + k])
                if c != 0.0:
                    acc_j += q(c) * pz[k]
            if not acc_j.is_zero():
                acc_i += py[j] * acc_j
        if not acc_i.is_zero():
            tot += px[i] * acc_i
    return tot


def cut_cells(N):
    out = []
    a2 = Fr(1, 4)
    for ex in range(N):
        for ey in range(N):
            for ez in range(N):
                lo = Fr(ex * ex + ey * ey + ez * ez, N * N)
                hi = Fr((ex + 1) ** 2 + (ey + 1) ** 2 + (ez + 1) ** 2, N * N)
                if lo < a2 < hi:
                    out.append((ex, ey, ez))
    return out


def full_cell_boxes(N):
    """the cells entirely outside the open ball, as level-0 boxes in curved3d's box format"""
    b = pf.build_boxes(N, 0, A2)
    return b[b[:, 7] == 0]


def _run(cells, fn_cell, checkpoint, verbose, label):
    done = {}
    if checkpoint and os.path.exists(checkpoint):
        with open(checkpoint) as fh:
            done = json.load(fh)
    t0 = time.time()
    prec = fctx.prec
    fctx.prec = 160
    try:
        for n, cell in enumerate(cells):
            key = '%d,%d,%d' % cell
            if key in done:
                continue
            b = fn_cell(cell)
            done[key] = str(arb_bounds(b)[1])
            if checkpoint:
                with open(checkpoint + '.tmp', 'w') as fh:
                    json.dump(done, fh)
                os.replace(checkpoint + '.tmp', checkpoint)
            if verbose and n % 10 == 0:
                print(f"    {label} cut cells {n}/{len(cells)}  {time.time() - t0:.0f}s", flush=True)
    finally:
        fctx.prec = prec
    return sum((Fr(done['%d,%d,%d' % c]) for c in cells), Fr(0))


def dual_cut_upper(S, Psi, gam, checkpoint=None, verbose=False):
    """rational upper bound of the integral of tau:S:tau over the cut cells minus the ball"""
    px, py, pz = Pieces3(S.S[0], X3), Pieces3(S.S[1], Y3), Pieces3(S.S[2], Z3)
    nu = q(pf.NU)
    ng = len(S.gexp)

    def one(cell):
        ex, ey, ez = cell
        fx, fy, fz = [int(S.S[k].first[e]) for k, e in enumerate(cell)]
        rho = uni3(eta_poly(ex, S.N), X3) + Y3 * Y3 + Z3 * Z3 - q(Fr(1, 4))
        R = rho * rho
        Phi = [R * tensor_cell(px(ex), py(ey), pz(ez), Psi[c], fx, fy, fz) for c in range(6)]
        names = 'xyz'
        H = []
        for (k, m) in COMPS:
            H.append([Phi[c].derivative(names[k]).derivative(names[m]) for c in range(6)])   # H[pair][comp]
        tau = []
        for s in range(6):
            t = ZERO3
            for p in range(6):
                for c in range(6):
                    g = G[s, p, c]
                    if g:
                        t += H[c][p] * int(g)
            tau.append(t)
        if Fr(ey, S.N) >= Fr(1, 2):
            tau[0] += sum((q(float(gam[i])) * (Y3 - q(Fr(1, 2))) ** e for i, e in enumerate(S.gexp)), ZERO3)
        if Fr(ez, S.N) >= Fr(1, 2):
            tau[0] += sum((q(float(gam[ng + i])) * (Z3 - q(Fr(1, 2))) ** e for i, e in enumerate(S.gexp)), ZERO3)
        tr = tau[0] + tau[1] + tau[2]
        dens = (tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]
                + 2 * (tau[3] * tau[3] + tau[4] * tau[4] + tau[5] * tau[5])) * (1 + nu) - nu * tr * tr
        box = [Fr(ex, S.N), Fr(ex + 1, S.N), Fr(ey, S.N), Fr(ey + 1, S.N), Fr(ez, S.N), Fr(ez + 1, S.N)]
        return outside_integral(dens, box)

    return _run(cut_cells(S.N), one, checkpoint, verbose, 'dual')


def primal_cut_upper(S, V, checkpoint=None, verbose=False):
    """rational upper bound of the integral of eps:C:eps over the cut cells minus the ball"""
    px, py, pz = Pieces3(S.U[0], X3), Pieces3(S.U[1], Y3), Pieces3(S.U[2], Z3)
    lam, mu = q(pf.LAM), q(pf.MU)

    def one(cell):
        ex, ey, ez = cell
        fx, fy, fz = [int(S.U[k].first[e]) for k, e in enumerate(cell)]
        v = [tensor_cell(px(ex), py(ey), pz(ez), V[a], fx, fy, fz) for a in range(3)]
        v[0] = v[0] + X3 * q(Fr(1, 2))
        d = [[v[a].derivative(n) for n in 'xyz'] for a in range(3)]
        e = [d[0][0], d[1][1], d[2][2], d[0][1] + d[1][0], d[1][2] + d[2][1], d[0][2] + d[2][0]]
        tr = e[0] + e[1] + e[2]
        dens = lam * tr * tr + 2 * mu * (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]) + mu * (e[3] * e[3] + e[4] * e[4]
                                                                                          + e[5] * e[5])
        box = [Fr(ex, S.N), Fr(ex + 1, S.N), Fr(ey, S.N), Fr(ey + 1, S.N), Fr(ez, S.N), Fr(ez + 1, S.N)]
        return outside_integral(dens, box)

    return _run(cut_cells(S.N), one, checkpoint, verbose, 'primal')
