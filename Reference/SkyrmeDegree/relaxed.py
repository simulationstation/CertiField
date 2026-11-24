#!/usr/bin/env python3
"""The relaxed level-set problem for Skyrme fields of degree B, in ordinary floating point.

    R_B(mu) = inf  int_0^pi  A psi(B c / A, M / A) + V M  dtau ,     A >= g(v),  v' = M >= 0,  v(0) = 0,
    c = 4 pi sin^2 tau,   V = 2 mu^2 (1 + cos tau),   g(v) = (36 pi)^(1/3) v^(2/3),
    psi(D, nu) = sup_{kappa >= 0} [ h**_kappa(D) - kappa nu ],   h_kappa(D) = 2 sqrt((1 + 2D)(D^2 + 2D + kappa)).

docs/SKYRME_DEGREE_TWO.md proves E(U) >= R_B for every finite-energy map of degree +-B.

NOTHING IN THIS FILE IS RIGOROUS.  It integrates the extremals of the problem (Pontryagin system) and finds the one
that satisfies the end condition by bisection.  Its results are used in two ways only:
  * the volume path v_ref(tau) of the extremal places the grids of the certified dynamic programme
    (make_grids.py); grids affect how tight the certified bound is, never whether it is valid;
  * the cost of the extremal is an independent estimate, from above, of the number the certified programme bounds
    from below (check_paths.py).

The Hamiltonian is H(tau, v, kappa) = min_{A >= g(v)} A h**_kappa(B c / A), kappa = V + dW/dv.  Writing
I(D) = h** - D dh**/dD (the intercept of the tangent at D):
  * where I(B c / g(v)) > 0 the area constraint is active, A = g(v) = 4 pi r^2:
        r' = dh**/dkappa (D),    kappa' = V' - (2/r) I(D),    D = B sin^2 tau / r^2 ;
  * otherwise more area is better and A = B c / D*, D* the point where the tangent passes through the origin:
        v' = B c h_kappa'(D*) / D*   (derivative in kappa),    kappa' = V' .
For B = 1 the extremal starts on the first branch with r ~ tau / a (the hedgehog branch; kappa(0) = a^4 - a^2);
for B >= 2 it starts on the second, with a free multiplier kappa(0) - V(0) = lam0 >= 0.
"""
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

K1 = (3 - math.sqrt(5)) / 2          # h_kappa is convex on [0, oo) iff K1 <= kappa <= K2
K2 = (3 + math.sqrt(5)) / 2
FOURPI = 4 * math.pi
UNIT = 12 * math.pi ** 2


def q(k, D):
    return (1 + 2 * D) * (D * D + 2 * D + k)


def h(k, D):
    return 2 * math.sqrt(q(k, D))


def h_D(k, D):
    return (2 * (D * D + 2 * D + k) + (1 + 2 * D) * (2 * D + 2)) / math.sqrt(q(k, D))


def h_k(k, D):
    """dh/dkappa; the optimal nu"""
    return math.sqrt((1 + 2 * D) / (D * D + 2 * D + k))


def intercept(k, D):
    return h(k, D) - D * h_D(k, D)


def D_infl(k):
    """the inflection point of h_kappa (0 if it is convex)"""
    n0 = -4 * (k * k - 3 * k + 1)
    if n0 >= 0:
        return 0.0
    f = lambda D: 12 * D ** 4 + 40 * D ** 3 + 24 * (1 + k) * D * D + 24 * k * D + n0
    hi = 1.0
    while f(hi) < 0:
        hi *= 2
    return brentq(f, 0.0, hi, xtol=1e-15, rtol=1e-14)


def D_t(k):
    """tangent point of the line through (0, 2 sqrt(kappa)); 0 if h_kappa is convex"""
    if K1 <= k <= K2:
        return 0.0
    if k <= 0:
        return 1.0
    di = D_infl(k)
    f = lambda D: intercept(k, D) - 2 * math.sqrt(k)
    hi = max(2 * di, 1e-6)
    while f(hi) > 0:
        hi *= 2
    return brentq(f, di, hi, xtol=1e-15, rtol=1e-13)


def D_star(k):
    """the point where the tangent of h** passes through the origin (h/D is minimal there)"""
    lo = D_t(k) if not (K1 <= k <= K2) else 1e-12
    lo = max(lo, 1e-12)
    hi = max(2 * lo, 1.0)
    while intercept(k, hi) > 0:
        hi *= 2
    if intercept(k, lo) <= 0:
        return lo
    return brentq(lambda D: intercept(k, D), lo, hi, xtol=1e-15, rtol=1e-13)


def envelope(k, D):
    """(h**, dh**/dkappa, intercept, fraction of the area that carries flux or None)"""
    if k > 0:
        dt = D_t(k)
        if D < dt:
            lam = D / dt
            y = math.sqrt(k)
            return (1 - lam) * 2 * y + lam * h(k, dt), (1 - lam) / y + lam * h_k(k, dt), 2 * y, lam
    else:
        if D < 1.0:
            return 6 * D, float('inf'), 0.0, D
    return h(k, D), h_k(k, D), intercept(k, D), None


def Vpot(mu, tau):
    return 2 * mu * mu * (1 + math.cos(tau))


def shoot(B, par, mu, tau0=1e-4, eps=2e-4, rtol=1e-10, atol=1e-13):
    """integrate the extremal from tau0.  B = 1: par = a (central slope);  B >= 2: par = lam0 (multiplier at 0).
    state: r(tau), kappa(tau), cost(tau).  Stops early if the multiplier kappa - V reaches zero or r blows up."""
    if B == 1:
        a = par
        y = [tau0 / a, a ** 4 - a * a, 0.0]
    else:
        k0 = par + Vpot(mu, tau0)
        ds = D_star(max(k0, 1e-300))
        v0 = B * FOURPI / 3 * tau0 ** 3 * h_k(max(k0, 1e-300), ds) / ds
        y = [(3 * v0 / FOURPI) ** (1 / 3), k0, B * FOURPI / 3 * tau0 ** 3 * h(max(k0, 1e-300), ds) / ds - par * v0]

    def rhs(tau, y):
        r, kap, _ = y
        r = max(r, 1e-12)
        V = Vpot(mu, tau)
        kap = max(kap, V, 1e-300)
        c = B * FOURPI * math.sin(tau) ** 2
        D = c / (FOURPI * r * r)
        hs, hk, I, lam = envelope(kap, D)
        if I < 0:                                # the area constraint is not active
            ds = D_star(kap)
            rho = h(kap, ds) / ds
            vp = c * h_k(kap, ds) / ds
            return [vp / (FOURPI * r * r), -2 * mu * mu * math.sin(tau), c * rho - (kap - V) * vp]
        return [hk, -2 * mu * mu * math.sin(tau) - 2.0 / r * I, FOURPI * r * r * (hs - (kap - V) * hk)]

    t_end = math.pi - eps
    ts, ys = [tau0], [list(y)]
    t = tau0
    early = False
    while t < t_end - 1e-15:
        dt = min(0.02, 0.2 * (math.pi - t), t_end - t)
        sol = solve_ivp(rhs, (t, t + dt), y, method='LSODA', rtol=rtol, atol=atol)
        t = float(sol.t[-1])
        y = [float(x) for x in sol.y[:, -1]]
        ts.append(t)
        ys.append(list(y))
        if (not np.all(np.isfinite(y))) or y[1] - Vpot(mu, t) <= 0.0 or y[0] > 1e7 or sol.status != 0:
            early = True
            break
    ts = np.array(ts)
    ys = np.array(ys).T
    return dict(B=B, par=par, mu=mu, t=ts, r=ys[0], kap=ys[1], cost=ys[2], early=early)


def discriminant(res):
    """sign for the bisection.  In the two-phase tail of the massless problem r sqrt(kappa) - (pi - tau) is
    constant along extremals and must vanish; with a mass: early stop (-1) against a multiplier left over (+1)."""
    if res['early']:
        return -1.0
    if res['mu'] == 0.0:
        return res['r'][-1] * math.sqrt(max(res['kap'][-1], 0.0)) - (math.pi - res['t'][-1])
    return 1.0


def tail_cost(res):
    """the cost beyond the last point: flux tube 6 B int c, and for mu = 0 the radial term 4 pi r l^2"""
    ell = math.pi - res['t'][-1]
    out = 6 * res['B'] * FOURPI * ell ** 3 / 3
    if res['mu'] == 0.0:
        out += FOURPI * res['r'][-1] * ell * ell
    return out


def solve(B, mu, lo=None, hi=None, n_scan=25, verbose=False):
    """the extremal that satisfies the end condition; returns the shot that reaches pi - eps, with the value"""
    if lo is None:
        lo, hi = (1.0, 4.0) if B == 1 else (0.0, 40.0)
    f = lambda p: discriminant(shoot(B, p, mu))
    grid = np.linspace(lo, hi, n_scan)
    vals = [f(p) for p in grid]
    if verbose:
        print("  scan:", " ".join(f"{p:.3f}:{'+' if v > 0 else '-'}" for p, v in zip(grid, vals)), flush=True)
    br = None
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] < 0:
            br = (grid[i], grid[i + 1], vals[i])
            break
    if br is None:
        return None
    a, b, va = br
    for _ in range(60):
        m = 0.5 * (a + b)
        vm = f(m)
        if vm * va > 0:
            a = m
        else:
            b = m
    res = shoot(B, a if va > 0 else b, mu)       # the side that reaches pi - eps
    if res['early']:
        res = shoot(B, b if va > 0 else a, mu)
    for d in (2e-9, -2e-9, 2e-8, -2e-8, 2e-7, -2e-7, 2e-6, -2e-6):
        if not res['early']:
            break
        res = shoot(B, 0.5 * (a + b) * (1 + d) + d, mu)
    total = res['cost'][-1] + tail_cost(res)
    res.update(value=total, value_units=total / UNIT, total=total)
    return res


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    # usage: python relaxed.py [B mu]...     default: the known degree-one values and degree two
    cases = [(1, 0.0), (1, 1.0), (2, 0.0), (2, 1.0)]
    if len(sys.argv) > 2:
        cases = [(int(sys.argv[i]), float(sys.argv[i + 1])) for i in range(1, len(sys.argv) - 1, 2)]
    for B, mu in cases:
        s = solve(B, mu, verbose=True)
        if s is None:
            print(f"B = {B}, mu = {mu}: no bracket")
            continue
        print(f"B = {B}, mu = {mu}:  parameter {s['par']:.10f}   value/12pi^2 = {s['value_units']:.7f}   "
              f"(reached tau = pi - {math.pi - s['t'][-1]:.1e}, r = {s['r'][-1]:.3f}, kappa = {s['kap'][-1]:.3e})", flush=True)
