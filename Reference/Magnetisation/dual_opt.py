#!/usr/bin/env python3
"""Heuristic search for a good certificate (phi, u) of Proposition 5 (docs/MAGNETISATION.md), degree Q < 0, kappa > 0.

NOTHING HERE IS TRUSTED.  The output is a magnet.Certificate with rational data; its validity is decided by
Certificate.check() in exact arithmetic.

  1. phi is a C^1 piecewise cubic on a mesh that is geometric in s = 1 - t; the unknowns are its nodal values and
     slopes.  For given phi the best u at each t is the smaller root u_-(t) of Pi_0 = 0, and the objective
     1/2 int u_- dt is convex in phi.  It is minimised by L-BFGS with numerical gradients, with (1 - eps) V in place
     of V (a margin for step 2) and a penalty where no root exists.
  2. On every cell a polynomial u slightly above u_-, with the margin raised until the exact check of the cell
     passes.
"""
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import minimize

import magnet as M

_GX, _GW = np.polynomial.legendre.leggauss(10)
_GX, _GW = 0.5 * (_GX + 1), 0.5 * _GW                      # Gauss points and weights on [0, 1]
_H = np.array([2 * _GX ** 3 - 3 * _GX ** 2 + 1, _GX ** 3 - 2 * _GX ** 2 + _GX, -2 * _GX ** 3 + 3 * _GX ** 2, _GX ** 3 - _GX ** 2])
_dH = np.array([6 * _GX ** 2 - 6 * _GX, 3 * _GX ** 2 - 4 * _GX + 1, -6 * _GX ** 2 + 6 * _GX, 3 * _GX ** 2 - 2 * _GX])


def roots(A, kappa, K, h, t, phi, dphi, eps=0.0):
    """smaller and larger root u of Pi_0 = 0 at the points t, and the discriminant (floats)"""
    a1 = (phi + kappa) + (1 + t) * dphi                 # ((1 + t)(phi + kappa))'
    b1 = -(phi - kappa) + (1 - t) * dphi                # ((1 - t)(phi - kappa))'
    P = ((1 - t) * a1) ** 2 / 4
    R = ((1 + t) * b1) ** 2 / 4
    V = (1 - eps) * (K * (1 - t * t) + h * (1 - t))
    B = 4 * A * V + P - R
    D = B * B - 16 * A * V * P
    sq = np.sqrt(np.maximum(D, 0.0))
    return (B - sq) / (2 * V), (B + sq) / (2 * V), D, B, V


def make_mesh(h_eff, ratio=0.75, s_last=2e-3):
    s = [2.0]
    while s[-1] * ratio > s_last:
        s.append(s[-1] * ratio)
    return np.array([1.0 - x for x in s] + [1.0])


def objective(z, tk, A, kappa, K, h, eps):
    n = tk.size - 1
    hk = np.diff(tk)
    sk = 1.0 - tk
    val = np.concatenate([z[:n - 1], [kappa, kappa]])                 # phi = kappa on the last cell
    slope = np.concatenate([z[n - 1:2 * n - 2] / sk[:n - 1], [0.0, 0.0]])
    t = tk[:-1, None] + hk[:, None] * _GX[None, :]
    phi = val[:-1, None] * _H[0] + (hk * slope[:-1])[:, None] * _H[1] + val[1:, None] * _H[2] + (hk * slope[1:])[:, None] * _H[3]
    dphi = (val[:-1, None] * _dH[0] + (hk * slope[:-1])[:, None] * _dH[1] + val[1:, None] * _dH[2] + (hk * slope[1:])[:, None] * _dH[3]) / hk[:, None]
    um, up, D, B, V = roots(A, kappa, K, h, t, phi, dphi, eps)
    bad = (D < 0) | (B <= 0)
    pen = np.where(bad, np.abs(B) / (2 * V) + 50.0 * np.maximum(-D, 0.0) / (4 * A * V) ** 2 * 4 * A + 10.0 * np.maximum(-B, 0.0) / V, 0.0)
    u = np.where(bad, pen, um)
    return 0.5 * np.sum(hk[:, None] * _GW[None, :] * u), int(bad.sum())


def optimise(model, ratio=0.75, eps=2e-3, maxiter=4000, stages=(2e-3, 2e-4, 2e-5), verbose=False):
    """continuation in the mesh: phi = kappa on the last cell [1 - s_last, 1]; s_last is lowered stage by stage and
    every stage starts from the previous field, which it contains"""
    A, kappa, K, h = float(model.A), float(model.kappa), float(model.K), float(model.h)
    assert kappa > 0
    z = None
    for s_last in stages:
        tk = make_mesh(h + 2 * K, ratio, s_last)
        n = tk.size - 1
        z0 = np.concatenate([np.full(n - 1, kappa), np.zeros(n - 1)])        # phi = kappa
        if z is not None:
            m = (z.size // 2)
            z0[:m] = z[:m]
            z0[n - 1:n - 1 + m] = z[m:]
        f = lambda zz: objective(zz, tk, A, kappa, K, h, eps)[0]
        z = z0
        best = objective(z, tk, A, kappa, K, h, eps)
        for restart in range(6):
            res = minimize(f, z, method='L-BFGS-B', options=dict(maxiter=maxiter, maxfun=400 * maxiter, ftol=1e-15, gtol=1e-10, eps=1e-7))
            cur = objective(res.x, tk, A, kappa, K, h, eps)
            if cur[1] == 0 and cur[0] < best[0] - 1e-9:
                z, best = res.x, cur
            else:
                if cur[1] == 0 and cur[0] < best[0]:
                    z, best = res.x, cur
                break
        if verbose:
            print(f"   stage s_last = {s_last:g}: {n} cells, 1/2 int u = {best[0]:.6f}, infeasible points {best[1]}", flush=True)
    return tk, z, best[0], best[1]


def build_certificate(model, ratio=0.75, eps=2e-3, maxiter=4000, stages=(2e-3, 2e-4, 2e-5), verbose=False):
    """a certificate for degree Q < 0.  Returns (certificate, info); the certificate still has to pass .check()."""
    A, kappa, K, h = float(model.A), float(model.kappa), float(model.K), float(model.h)
    tk, z, val, nbad = optimise(model, ratio, eps, maxiter, stages, verbose)
    n = tk.size - 1
    sk = 1.0 - tk
    mesh = [Fr(float(t)) for t in tk]
    mesh[0], mesh[-1] = Fr(-1), Fr(1)
    best = None
    # damping towards phi = kappa: the set of fields for which a root exists is convex, and phi = kappa lies strictly
    # inside it when h > kappa^2/(2A); a small step towards it repairs cells where the optimiser's field has no root
    for theta in (0.0, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1):
        zv = kappa + (1 - theta) * (z[:n - 1] - kappa)
        zs = (1 - theta) * z[n - 1:2 * n - 2] / sk[:n - 1]
        vals = [Fr(float(v)) for v in zv] + [Fr(model.kappa), Fr(model.kappa)]
        slopes = [Fr(float(v)) for v in zs] + [Fr(0), Fr(0)]
        cert0 = M.Certificate.from_hermite(model, mesh, vals, slopes, [[Fr(0)]] * n)
        us, bad = _fit_u(model, mesh, cert0, A, kappa, K, h)
        cert = M.Certificate(model, mesh, cert0.phi, us)
        if best is None or bad < best[1]:
            best = (cert, bad, theta)
        if bad == 0:
            break
    cert, bad, theta = best
    return cert, dict(optimiser_value=2 * A - val, cells=n, failed_cells=bad, infeasible_points=nbad, damping=theta)


def _fit_u(model, mesh, cert0, A, kappa, K, h):
    n = len(mesh) - 1
    us, bad = [], 0
    for i in range(n):
        t0, t1 = float(mesh[i]), float(mesh[i + 1])
        last = (i == n - 1)
        deg = 4
        nodes = t0 + (t1 - t0) * 0.5 * (1 - np.cos(np.pi * (np.arange(deg + 1) + 0.5) / (deg + 1)))     # Chebyshev points
        pphi = np.array([float(c) for c in cert0.phi[i]])
        dpphi = pphi[1:] * np.arange(1, len(pphi))
        ok = False
        if last:
            # phi = kappa on the last cell: P = kappa^2 (1 - t)^2, R = 0, and u = kappa^2 (1 - t)/v0 with
            # v0 <= h + K (1 + t) on the cell satisfies Pi_0 <= 0
            v0 = model.h + model.K * (1 + (mesh[i] if model.K >= 0 else Fr(1)))
            us.append([model.kappa ** 2 / v0, -model.kappa ** 2 / v0])
            continue
        um, up, D, B, V = roots(A, kappa, K, h, nodes, np.polyval(pphi[::-1], nodes), np.polyval(dpphi[::-1], nodes))
        if np.all(D > 0) and np.all(B > 0):
            # first a small relative margin above the smaller root; if that fails, a point inside the band
            # between the two roots (any u between them satisfies Pi_0 <= 0)
            for margin in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, -0.02, -0.1, -0.3, -0.5):
                lin = [-mesh[i] / (mesh[i + 1] - mesh[i]), 1 / (mesh[i + 1] - mesh[i])]
                target = um * (1 + margin) + 1e-3 * margin * um.max() if margin > 0 else um + (-margin) * (up - um)
                c = np.polyfit((nodes - t0) / (t1 - t0), target, deg)[::-1]
                upoly = [Fr(0)]
                for ck in reversed([Fr(float(v)) for v in c]):
                    upoly = M.padd(M.pmul(upoly, lin), [ck])
                trial = M.Certificate(model, [mesh[i], mesh[i + 1]], [cert0.phi[i]], [upoly])
                Pp, Rr, VG, Pi, fourA_u = trial.polynomials(0)
                if M.nonneg(upoly, mesh[i], mesh[i + 1]) and M.nonneg(fourA_u, mesh[i], mesh[i + 1]) and M.nonneg(M.pscale(Pi, -1), mesh[i], mesh[i + 1]):
                    ok = True
                    break
        if not ok:
            bad += 1
            upoly = [4 * model.A]                          # placeholder: the certificate will be rejected
        us.append(upoly)
    return us, bad


if __name__ == '__main__':
    import sys
    import time
    sys.stdout.reconfigure(encoding="utf-8")
    for h in [Fr(x) for x in (sys.argv[1:] or ['2'])]:
        t0 = time.time()
        model = M.Model(A=Fr(1, 2), kappa=1, h=h)
        cert, info = build_certificate(model, verbose=True)
        t1 = time.time()
        ok = cert.check(verbose=True)
        c, d = cert.bound(-1)
        print(f"h = {h}: optimiser {info['optimiser_value']:.6f};  damping {info['damping']:g};  certificate {'VALID' if ok else 'rejected'}: E/(4 pi) >= {float(c) / 4:.6f}"
              f"   cells {info['cells']}, failed {info['failed_cells']}   ({t1 - t0:.0f} s search, {time.time() - t1:.0f} s exact check)", flush=True)
