#!/usr/bin/env python3
"""Known answer for the stray-field bounds of the thin-film model (docs/MAGNETISATION.md, Section 10).

For a magnetisation m(x1, x2) of a film of thickness t (m = e3 far away), the stray-field energy relative to the
uniform state is, by the Fourier formula,
    E_ms[m] - E_ms[e3] = t Int w c + t/2 Int sigma(|k| t) |w^|^2 + t/2 Int (1 - sigma(|k| t)) |k^.m'^|^2,
w = m3 - 1, c = 1 on the infinite film, sigma(s) = (1 - exp(-s))/s.  With Parseval,
    E_ms[m] - E_ms[e3] = -t/2 Int (1 - m3^2) + corr,
    corr = -t/2 Int (1 - sigma) |w^|^2 + t/2 Int (1 - sigma) |k^.m'^|^2.
For axisymmetric textures the transforms are Hankel transforms, computed here by quadrature (floating point,
independent of the bounds).  The certified bounds of the model say
    -t^2/4 |w| |grad m3|  <=  corr  <=  t^2/4 |div m'| |m'|,
and corr / t -> 0 as t -> 0 (the local shape anisotropy).  Bloch textures have div m' = 0, Neel textures do not."""
import sys

import numpy as np
from scipy.special import j0, j1

sys.stdout.reconfigure(encoding="utf-8")


def texture(R, L, kind, n=20000):
    r = np.linspace(0.0, L, n + 1)
    chi = np.clip(1 - (r / L) ** 2, 0, None) ** 2
    th = 2 * np.arctan2(R * chi, r)
    dth = np.gradient(th, r)
    s = np.sin(th)
    w = np.cos(th) - 1
    rr = np.where(r > 0, r, 1.0)
    divm = np.where(r > 0, np.gradient(s, r) + s / rr, 2 * dth) if kind == 'neel' else np.zeros_like(r)
    return r, w, s, dth, divm


def trap(f, r):
    return np.trapezoid(f * 2 * np.pi * r, r)


def transforms(r, w, s, kind, kmax, nk=6000):
    k = np.linspace(1e-6, kmax, nk)
    wk = np.array([np.trapezoid(w * j0(kk * r) * r, r) for kk in k]) * 2 * np.pi
    mk = np.array([np.trapezoid(s * j1(kk * r) * r, r) for kk in k]) * 2 * np.pi if kind == 'neel' else 0 * k
    return k, wk, mk


fails = 0
for kind in ('bloch', 'neel'):
    for R, L in ((0.5, 6.0), (0.2, 4.0)):
        r, w, s, dth, divm = texture(R, L, kind)
        k, wk, mk = transforms(r, w, s, kind, kmax=400 / R)
        parseval = abs(np.trapezoid(wk ** 2 * k, k) / (2 * np.pi) - trap(w ** 2, r)) / trap(w ** 2, r)
        nw, ngm3 = np.sqrt(trap(w ** 2, r)), np.sqrt(trap((s * dth) ** 2, r))
        ndiv, nmp = np.sqrt(trap(divm ** 2, r)), np.sqrt(trap(s ** 2, r))
        for t in (0.4, 0.1, 0.025):
            one_minus_sig = 1 - (1 - np.exp(-k * t)) / (k * t)
            corr = (t / 2) * np.trapezoid(one_minus_sig * (mk ** 2 - wk ** 2) * k, k) / (2 * np.pi)
            lo, hi = -t ** 2 / 4 * nw * ngm3, t ** 2 / 4 * ndiv * nmp
            ok = lo <= corr <= hi
            fails += not ok
            print(f"  {kind:5s} R = {R}, t = {t:5.3f}:  {lo:+.7f} <= corr {corr:+.7f} <= {hi:+.7f}   corr/t {corr / t:+.6f}   "
                  f"(Parseval error {parseval:.1e})   {'ok' if ok else '**** FAILED ****'}", flush=True)
print("ALL TESTS PASSED" if fails == 0 else f"FAILURES: {fails}")
sys.exit(1 if fails else 0)
