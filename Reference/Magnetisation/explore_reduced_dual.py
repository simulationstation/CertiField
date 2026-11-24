"""Feasibility study (floating point, nothing rigorous), second version.
Model  E = int A |grad m|^2 + kappa m.curl m + h (1 - m3) + K (1 - m3^2),  A = 1/2, kappa = 1, degree -1.
 (1) energy of the axisymmetric skyrmion (Newton on the discrete radial energy);
 (2) the translation-invariant, rotation-equivariant dual: minimise I = 1/2 int_{-1}^{1} u dt subject to
        P/u + R/(4A - u) <= V(t),  P = ((1-t) a'/2)^2,  R = (kappa (1+t) - a/(1+t) + (1-t) a'/2)^2,  a(-1) = 0, a(1) = 4 kappa,
     by dynamic programming over (x = ln(1 - t), a);  the bound is  E >= 4 pi (2A - I)."""
import sys
import time
import numpy as np
from scipy.linalg import solveh_banded

sys.stdout.reconfigure(encoding="utf-8")
A = 0.5
kappa = 1.0


def skyrmion(h, K=0.0, n=3000, rmax=None, verbose=False):
    rmax = rmax or 30.0 / np.sqrt(2 * (h + 2 * K))
    r = np.concatenate([[0.0], np.geomspace(1e-4 / (h + 1), rmax, n)])
    rm = 0.5 * (r[1:] + r[:-1])
    dr = np.diff(r)
    w = 2 * np.pi * rm * dr

    def parts(th):
        dth = np.diff(th) / dr
        tm = 0.5 * (th[1:] + th[:-1])
        s, c = np.sin(tm), np.cos(tm)
        dens = A * (dth ** 2 + s ** 2 / rm ** 2) + kappa * (dth + s * c / rm) + h * (1 - c) + K * s * s
        Ld = w * (2 * A * dth + kappa)
        Lm = w * (A * 2 * s * c / rm ** 2 + kappa * (c * c - s * s) / rm + h * s + 2 * K * s * c)
        Ldd = w * 2 * A
        Lmm = w * (2 * A * (c * c - s * s) / rm ** 2 - 4 * kappa * s * c / rm + h * c + 2 * K * (c * c - s * s))
        return np.sum(w * dens), Ld, Lm, Ldd, Lmm

    def energy(th):
        return parts(th)[0]

    R0 = kappa / ((h + 2 * K) * max(np.log(h + 2 * K + 2.0), 1.0))
    th = 2 * np.arctan2(R0, r) * np.exp(-np.sqrt(2 * (h + 2 * K)) * r * 0.7)
    th[0], th[-1] = np.pi, 0.0
    lam = 1e-3
    E = energy(th)
    for it in range(400):
        E, Ld, Lm, Ldd, Lmm = parts(th)
        g = np.zeros_like(th)
        g[1:] += Ld / dr + 0.5 * Lm
        g[:-1] += -Ld / dr + 0.5 * Lm
        # Hessian: diag and off-diagonal
        d = np.zeros_like(th)
        d[1:] += Ldd / dr ** 2 + 0.25 * Lmm
        d[:-1] += Ldd / dr ** 2 + 0.25 * Lmm
        o = -Ldd / dr ** 2 + 0.25 * Lmm
        gi = g[1:-1]
        if np.max(np.abs(gi)) < 1e-11:
            break
        while True:
            ab = np.zeros((2, th.size - 2))
            ab[1] = d[1:-1] + lam * w.mean()
            ab[0, 1:] = o[1:-1]
            try:
                step = solveh_banded(ab, -gi)
            except np.linalg.LinAlgError:
                lam *= 10
                continue
            tn = th.copy()
            tn[1:-1] += step
            En = energy(tn)
            if En <= E + 1e-4 * gi.dot(step):
                th = tn
                lam = max(lam / 5, 1e-12)
                break
            lam *= 10
            if lam > 1e12:
                break
        if lam > 1e12:
            break
    k = np.argmax(th < np.pi / 2)
    if verbose:
        print("   newton its", it, "max grad", np.max(np.abs(gi)))
    return E / (4 * np.pi), r[k], r, th


def dual_dp(h, K=0.0, nx=260, na=1500, smin=1e-8, verbose=False):
    """dynamic programme from s = smin (a = 4 kappa) to s = 2 (a = 0)"""
    x = np.linspace(np.log(smin), np.log(2.0), nx + 1)
    dx = x[1] - x[0]
    # grid of a: fine near 4 kappa
    al = np.concatenate([[0.0], np.geomspace(1e-9, 4 * kappa, na)])           # alpha = 4 kappa - a
    agrid = 4 * kappa - al                                                        # decreasing from 4 kappa to 0
    W = np.full(agrid.size, np.inf)
    W[0] = 0.0                                                                    # start: a = 4 kappa at smin
    back = []
    for j in range(nx):
        xm = x[j] + dx / 2
        s = np.exp(xm)
        V = h * s + K * s * (2 - s)
        a0 = agrid[:, None]                                                       # value at x_j
        a1 = agrid[None, :]                                                       # value at x_{j+1}
        am = 0.5 * (a0 + a1)
        g = -(a1 - a0) / dx / 2
        e = kappa * (2 - s) - am / (2 - s)
        P = g * g
        R = (e + g) ** 2
        B = 4 * A * V + P - R
        D = B * B - 16 * A * V * P
        ok = (D >= 0) & (B > 0)
        u = np.where(ok, (B - np.sqrt(np.where(ok, D, 0.0))) / (2 * V), np.inf)
        cost = W[:, None] + 0.5 * s * dx * u
        best = np.argmin(cost, axis=0)
        W = cost[best, np.arange(agrid.size)]
        back.append(best)
    # end: a = 0  (alpha = 4 kappa): last index
    I = W[-1]
    # recover the path
    idx = agrid.size - 1
    path = [idx]
    for j in range(nx - 1, -1, -1):
        idx = back[j][idx]
        path.append(idx)
    path = path[::-1]
    return I, x, agrid[path]


if __name__ == "__main__":
    print("kappa = 1, A = 1/2;  energies in units of 4 pi")
    print("   h      K   E_skyrmion  radius    1-1/h     dual(DP)   gap(simple)  gap(DP)     (1-E)*h*ln(h)   I*h*ln(h)")
    cases = [(0.9, 0), (1.0, 0), (1.25, 0), (1.5, 0), (2.0, 0), (3.0, 0), (5.0, 0), (10.0, 0), (20.0, 0), (50.0, 0), (100.0, 0)]
    for h, K in cases:
        t0 = time.time()
        E, R, r, th = skyrmion(h, K)
        t1 = time.time()
        I, x, a = dual_dp(h, K)
        simple = 1 - 1 / h
        opt = 2 * A - I
        lg = h * np.log(h) if h > 1 else float("nan")
        print(f"{h:6.2f} {K:5.2f}  {E:10.6f}  {R:7.4f}  {simple:9.5f}  {opt:9.5f}   {100 * (E - simple) / abs(E):8.2f} %  {100 * (E - opt) / abs(E):8.2f} %"
              f"    {(1 - E) * lg:8.4f}     {I * lg:8.4f}    ({t1 - t0:.0f} s, {time.time() - t1:.0f} s)", flush=True)
