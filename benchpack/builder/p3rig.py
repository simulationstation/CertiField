"""P3: rigorous upper bounds of Q_oct(tau) and a_oct(v) for stored fields (see p3float.py for the model).

Both integrands are non-negative, so they are integrated over a superset of the octant domain: the octree boxes of
curved3d.build_boxes (a box is dropped only if it lies in the closed ball, decided in exact integer arithmetic).  On
every box the integrand is a polynomial of degree <= 2 ps + 8 (stress) or 2 pu (strain) in each variable; tensor
Gauss rules with ps + 5 (resp. pu + 1) points are exact.  Nodes and weights are rational enclosures (curved3d's
gauss_enclosures); the basis functions and the weight rho are evaluated at the enclosed nodes in exact rational
arithmetic and enclosed in floats; everything else is outward-rounded interval arithmetic (curved3d's ia.I), and the
sums are rounded upwards.  The float evaluation (rigorous=False) uses the same code path on midpoints.
"""
import json
import os
import time
from fractions import Fraction as Fr

import numpy as np

import p3float as pf
from p3float import G, COMPS, VAR, build_boxes, A2, eta_poly
from splines import poly_der, gauss01_rig
from ia import I, frac_enclosure, sqr, sum_up, up
from curved3d import RI, ri_horner

VIDX = {v: i for i, v in enumerate(VAR)}


def _enc(r):
    return frac_enclosure(r.lo)[0], frac_enclosure(r.hi)[1]


class DirTables:
    """per direction: enclosures at the Gauss nodes of every dyadic sub-interval of every cell"""

    def __init__(self, Sp, nq, nder, direction, N, gexp=()):
        self.Sp, self.nq, self.nder, self.dir, self.N, self.gexp = Sp, nq, nder, direction, N, tuple(gexp)
        self.nodes, self.weights = gauss01_rig(nq)
        self.cache = {}

    def get(self, e, lev, j):
        key = (e, lev, j)
        if key in self.cache:
            return self.cache[key]
        Sp, nq = self.Sp, self.nq
        T = [RI((j + n.lo) / 2 ** lev, (j + n.hi) / 2 ** lev) for n in self.nodes]
        tab = Sp.rig_table(e, T, self.nder)                          # I [nq, p+1, nder+1]
        h, a = Sp.h[e], Sp.a[e]
        w = np.array([frac_enclosure(wq.hi * h / 2 ** lev)[1] for wq in self.weights])
        X = [RI(a + h * t.lo, a + h * t.hi) for t in T]               # physical coordinates (exact intervals)
        # weight function pieces: eta (x direction) or the square (y, z), with first and second derivatives
        if self.dir == 0:
            c = eta_poly(e, self.N)
        else:
            c = [Fr(0), Fr(0), Fr(1)]
        rho = []
        for d in range(3):
            lo_, hi_ = zip(*[_enc(ri_horner(c, x)) for x in X])
            rho.append(I(np.array(lo_), np.array(hi_)))
            c = poly_der(c)
        # particular field: (y - 1/2)^e for y >= 1/2 (cells above 1/2 only; 1/2 is a knot)
        gp = []
        if self.dir > 0 and Sp.a[e] >= Fr(1, 2):
            for ex_ in self.gexp:
                vals = []
                for x in X:
                    t = x - Fr(1, 2)
                    p = RI(1)
                    for _ in range(ex_):
                        p = p * t
                    vals.append(_enc(p))
                lo_, hi_ = zip(*vals)
                gp.append(I(np.array(lo_), np.array(hi_)))
        self.cache[key] = (tab, w, rho, gp)
        return self.cache[key]


def _stack(items):
    """stack a list of I arrays along a new first axis"""
    return I(np.stack([it.lo for it in items]), np.stack([it.hi for it in items]))


def _contract(C, BX, BY, BZ, variants):
    """C [nb, nc, p+1, p+1, p+1] floats; BX.. I [nb, nq, p+1, nder+1] -> dict variant -> I [nb, nc, nq, nq, nq]"""
    p1 = C.shape[2]
    need_z = sorted({v[2] for v in variants})
    need_yz = sorted({(v[1], v[2]) for v in variants})
    Az = {}
    for dz in need_z:
        acc = None
        for m in range(p1):
            b = BZ[:, :, m, dz]                                     # [nb, nqz]
            term = I(b.lo[:, None, None, None, :], b.hi[:, None, None, None, :]) * C[:, :, :, :, m, None]
            acc = term if acc is None else acc + term
        Az[dz] = acc                                                # [nb, nc, mx, my, qz]
    Ayz = {}
    for (dy, dz) in need_yz:
        acc = None
        for m in range(p1):
            b = BY[:, :, m, dy]                                     # [nb, nqy]
            term = Az[dz][:, :, :, m, None, :] * I(b.lo[:, None, None, :, None], b.hi[:, None, None, :, None])
            acc = term if acc is None else acc + term
        Ayz[(dy, dz)] = acc                                         # [nb, nc, mx, qy, qz]
    out = {}
    for v in variants:
        acc = None
        for m in range(p1):
            b = BX[:, :, m, v[0]]                                   # [nb, nqx]
            term = Ayz[(v[1], v[2])][:, :, m, None, :, :] * I(b.lo[:, None, :, None, None], b.hi[:, None, :, None, None])
            acc = term if acc is None else acc + term
        out[v] = acc                                                # [nb, nc, qx, qy, qz]
    return out


def _local(cfull, firsts, p1):
    """cfull [nc, nx, ny, nz], firsts [nb, 3] -> [nb, nc, p1, p1, p1]"""
    i = np.arange(p1)
    ix = firsts[:, 0][:, None, None, None] + i[None, :, None, None]
    iy = firsts[:, 1][:, None, None, None] + i[None, None, :, None]
    iz = firsts[:, 2][:, None, None, None] + i[None, None, None, :]
    return np.moveaxis(cfull[:, ix, iy, iz], 0, 1)


def _grid(a, axis):
    """I [nb, nq] -> broadcastable [nb, 1, nq?, nq?, nq?] along the given axis (0, 1, 2)"""
    sl = [slice(None), None, None, None, None]
    sl[2 + axis] = slice(None)
    return I(a.lo[tuple(sl)], a.hi[tuple(sl)])


def dual_upper(S, Psi, gam, L, batch=48, verbose=False, checkpoint=None):
    """rigorous upper bound of Q_oct for tau = g + inc(rho^2 Psi); Psi [6, nx, ny, nz], gam = (gamma_k, delta_k).
    checkpoint: a file to store the partial sums per batch (resume after a stop)."""
    nq, p1 = S.nqs, S.ps + 1
    tabs = [DirTables(S.S[k], nq, 2, k, S.N, S.gexp) for k in range(3)]
    boxes = build_boxes(S.N, L, A2)
    nu = pf.NU
    nu_i = I(*frac_enclosure(nu))
    one_nu = I(*frac_enclosure(1 + nu))
    ng = len(S.gexp)
    gamma, delta = [float(v) for v in gam[:ng]], [float(v) for v in gam[ng:]]
    done = {}
    if checkpoint and os.path.exists(checkpoint):
        with open(checkpoint) as fh:
            done = {int(k): v for k, v in json.load(fh).items()}
    total = []
    t0 = time.time()
    for bi, s0 in enumerate(range(0, len(boxes), batch)):
        if bi in done:
            total.append(done[bi])
            continue
        bx = boxes[s0:s0 + batch]
        nb = len(bx)
        T = [[tabs[k].get(int(b[k]), int(b[3]), int(b[4 + k])) for b in bx] for k in range(3)]
        B = [_stack([t[0] for t in T[k]]) for k in range(3)]                 # I [nb, nq, p1, 3]
        W = [np.stack([t[1] for t in T[k]]) for k in range(3)]               # [nb, nq]
        R = [[_stack([t[2][d] for t in T[k]]) for d in range(3)] for k in range(3)]   # rho pieces [nb, nq]
        firsts = np.stack([S.S[k].first[bx[:, k]] for k in range(3)], axis=1)
        C = _local(Psi, firsts, p1)
        A3 = _contract(C, B[0], B[1], B[2], VAR)
        E = [[_grid(R[k][d], k) for d in range(3)] for k in range(3)]       # E[k][d]: d-th derivative of piece k
        rho = E[0][0] + E[1][0] + E[2][0] - Fr(1, 4)
        Rr = sqr(rho)
        dR = [(rho * E[k][1]) * 2.0 for k in range(3)]
        Hphi = []
        for c, (k, m) in enumerate(COMPS):
            dk = tuple(1 if i == k else 0 for i in range(3))
            dm = tuple(1 if i == m else 0 for i in range(3))
            dkm = tuple(dk[i] + dm[i] for i in range(3))
            h = Rr * A3[dkm] + dR[k] * A3[dm] + dR[m] * A3[dk]
            d2 = (E[k][1] * E[m][1]) * 2.0
            if k == m:
                d2 = d2 + (rho * E[k][2]) * 2.0
            Hphi.append(h + d2 * A3[(0, 0, 0)])                             # [nb, 6(p), q, q, q]
        tau = [None] * 6
        for s in range(6):
            for p in range(6):
                for c in range(6):
                    gsc = G[s, p, c]
                    if gsc:
                        term = Hphi[c][:, p] * float(gsc)
                        tau[s] = term if tau[s] is None else tau[s] + term
        # particular field g(y, z)
        gsum = None
        for k, coefs in ((1, gamma), (2, delta)):
            for i, cf in enumerate(coefs):
                if cf == 0.0:
                    continue
                vals = []
                for t in T[k]:
                    vals.append(t[3][i] if t[3] else I(np.zeros(nq)))
                term = _grid(_stack(vals), k)[:, 0] * cf
                gsum = term if gsum is None else gsum + term
        if gsum is not None:
            tau[0] = tau[0] + gsum
        tr = tau[0] + tau[1] + tau[2]
        q = (sqr(tau[0]) + sqr(tau[1]) + sqr(tau[2]) + (sqr(tau[3]) + sqr(tau[4]) + sqr(tau[5])) * 2.0) * one_nu \
            - sqr(tr) * nu_i
        w = up(up(W[0][:, :, None, None] * W[1][:, None, :, None]) * W[2][:, None, None, :])
        qh = np.maximum(q.hi, 0.0)
        part = float(sum_up(up(qh * w).reshape(-1)))
        total.append(part)
        done[bi] = part
        if checkpoint and (bi % 20 == 0):
            _save(checkpoint, done)
        if verbose and bi % 50 == 0:
            print(f"    dual certificate {s0}/{len(boxes)} boxes  {time.time() - t0:.0f}s", flush=True)
    if checkpoint:
        _save(checkpoint, done)
    return float(sum_up(np.array(total)))


def primal_upper(S, V, L, batch=128, verbose=False, checkpoint=None):
    """rigorous upper bound of a_oct(v) for v = (x/2, 0, 0) + splines V [3, nx, ny, nz]"""
    nq, p1 = S.nqu, S.pu + 1
    tabs = [DirTables(S.U[k], nq, 1, k, S.N) for k in range(3)]
    boxes = build_boxes(S.N, L, A2)
    lam = I(*frac_enclosure(pf.LAM))
    mu = I(*frac_enclosure(pf.MU))
    done = {}
    if checkpoint and os.path.exists(checkpoint):
        with open(checkpoint) as fh:
            done = {int(k): v for k, v in json.load(fh).items()}
    total = []
    t0 = time.time()
    for bi, s0 in enumerate(range(0, len(boxes), batch)):
        if bi in done:
            total.append(done[bi])
            continue
        bx = boxes[s0:s0 + batch]
        T = [[tabs[k].get(int(b[k]), int(b[3]), int(b[4 + k])) for b in bx] for k in range(3)]
        B = [_stack([t[0] for t in T[k]]) for k in range(3)]
        W = [np.stack([t[1] for t in T[k]]) for k in range(3)]
        firsts = np.stack([S.U[k].first[bx[:, k]] for k in range(3)], axis=1)
        C = _local(V, firsts, p1)
        A3 = _contract(C, B[0], B[1], B[2], [(1, 0, 0), (0, 1, 0), (0, 0, 1)])
        dv = [[A3[tuple(1 if i == k else 0 for i in range(3))][:, a] for k in range(3)] for a in range(3)]
        e = [dv[0][0] + 0.5, dv[1][1], dv[2][2], dv[0][1] + dv[1][0], dv[1][2] + dv[2][1], dv[0][2] + dv[2][0]]
        tr = e[0] + e[1] + e[2]
        q = sqr(tr) * lam + (sqr(e[0]) + sqr(e[1]) + sqr(e[2])) * (mu * 2.0) + (sqr(e[3]) + sqr(e[4]) + sqr(e[5])) * mu
        w = up(up(W[0][:, :, None, None] * W[1][:, None, :, None]) * W[2][:, None, None, :])
        part = float(sum_up(up(np.maximum(q.hi, 0.0) * w).reshape(-1)))
        total.append(part)
        done[bi] = part
        if checkpoint and bi % 20 == 0:
            _save(checkpoint, done)
        if verbose and bi % 50 == 0:
            print(f"    primal certificate {s0}/{len(boxes)} boxes  {time.time() - t0:.0f}s", flush=True)
    if checkpoint:
        _save(checkpoint, done)
    return float(sum_up(np.array(total)))


def _save(fn, done):
    with open(fn + '.tmp', 'w') as fh:
        json.dump({str(k): v for k, v in done.items()}, fh)
    os.replace(fn + '.tmp', fn)
