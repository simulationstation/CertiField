#!/usr/bin/env python3
"""Tests of the certified competitor (competitor.py): known answers, an independent quadrature, falsification.

1. Exchange only, degree -1: inf E = 8 pi A, not attained.  A cut-off Belavin-Polyakov map has certified degree -1 and
   a certified energy just above 8 pi A.
2. Critical family (A = 1/2, kappa = 1/2, K = -1/2, h = 1): inf E = 4 pi (1 - 2 r^2) = 2 pi on the plane, attained by a
   Belavin-Polyakov map of radius 2 r = 1.  The cut-off map has certified degree -1 and certified energy just above.
3. The certified upper bound is above an independent Gauss quadrature of the energy of the same normalised
   interpolant (a coarse field, all four energy terms), and the certified degree agrees with the quadrature.
4. Falsification: the stored lower-bound certificates (results/zeeman_h*.json) are below the certified energies of
   competitors of degree -1 for the same model.
"""
import glob
import json
import sys
from fractions import Fraction as Fr

import numpy as np

import magnet as M
import competitor as C

sys.stdout.reconfigure(encoding="utf-8")
FOURPI = 4 * np.pi
fails = 0


def report(tag, ok, extra=""):
    global fails
    fails += not ok
    print(f"  {tag:60s} {'ok' if ok else '**** FAILED ****'} {extra}", flush=True)


def quadrature(model, x1, x2, Mn, nq=12):
    """floating-point Gauss quadrature of E[m_h] and of the degree, from the definitions (no Bernstein forms)"""
    g, w = np.polynomial.legendre.leggauss(nq)
    g, w = 0.5 * (g + 1), 0.5 * w
    E = Q = 0.0
    B00, B10, B01, B11 = Mn[:-1, :-1], Mn[1:, :-1], Mn[:-1, 1:], Mn[1:, 1:]
    h1 = np.diff(x1)[:, None, None]
    h2 = np.diff(x2)[None, :, None]
    for xi, wx in zip(g, w):
        for eta, wy in zip(g, w):
            Mv = B00 * (1 - xi) * (1 - eta) + B10 * xi * (1 - eta) + B01 * (1 - xi) * eta + B11 * xi * eta
            Mx = ((B10 - B00) * (1 - eta) + (B11 - B01) * eta) / h1
            My = ((B01 - B00) * (1 - xi) + (B11 - B10) * xi) / h2
            nrm = np.linalg.norm(Mv, axis=-1, keepdims=True)
            m = Mv / nrm
            p = Mx / nrm - m * np.sum(m * Mx / nrm, axis=-1, keepdims=True)
            q = My / nrm - m * np.sum(m * My / nrm, axis=-1, keepdims=True)
            e1, e2 = np.array([1.0, 0, 0]), np.array([0, 1.0, 0])
            dens = (float(model.A) * (np.sum(p * p, -1) + np.sum(q * q, -1))
                    + float(model.kappa) * (np.sum(np.cross(m, e1) * p, -1) + np.sum(np.cross(m, e2) * q, -1))
                    + float(model.K) * (1 - m[..., 2] ** 2) + float(model.h) * (1 - m[..., 2]))
            area = (h1 * h2)[..., 0]
            E += np.sum(wx * wy * dens * area)
            Q += np.sum(wx * wy * np.sum(m * np.cross(p, q), -1) * area)
    return E, Q / FOURPI


print("1. exchange only (inf E = 8 pi A = 4 pi in degree -1, not attained)")
mod = M.Model(A=Fr(1, 2))
for R, n, beta in [(0.1, 40, 3.0), (0.1, 80, 3.0)]:
    x = C.graded_axis(1.0, n, beta)
    Mn = C.radial_field(x, x, C.belavin_polyakov(R, 1.0), 1.0)
    ub, parts = C.upper_bound(mod, x, x, Mn)
    Q, enc = C.degree(x, x, Mn)
    report(f"R = {R}, grid {2 * n}x{2 * n}: degree -1, 4 pi <= E <= 1.05 * 4 pi", Q == -1 and FOURPI <= ub <= 1.05 * FOURPI,
           f"E <= 4 pi * {ub / FOURPI:.5f}")

print("2. critical family r = 1/2 (inf E = 2 pi on the plane, Belavin-Polyakov map of radius 1)")
mod = M.Model(A=Fr(1, 2), kappa=Fr(1, 2), K=Fr(-1, 2), h=1)
cert = M.constant_certificate(mod, Fr(1, 2), [Fr(1, 2)])
lb = float(cert.lower(-1))
prev = None
for L, n, beta in [(12.0, 120, 4.0), (48.0, 200, 6.0)]:
    x = C.graded_axis(L, n, beta)
    Mn = C.radial_field(x, x, C.belavin_polyakov(1.0, L), L)
    ub, parts = C.upper_bound(mod, x, x, Mn)
    Q, enc = C.degree(x, x, Mn)
    ok = cert.check() and Q == -1 and lb <= ub and (prev is None or (ub < prev and ub <= 1.02 * lb))
    report(f"L = {L}, grid {2 * n}x{2 * n}: degree -1, 2 pi <= upper" + ("" if prev is None else ", decreasing, gap < 2 %"), ok,
           f"upper {ub:.6f} = 2 pi * {ub / (FOURPI / 2):.5f}")
    prev = ub

print("3. the enclosure is above an independent quadrature (coarse random field, all four terms)")
rng = np.random.default_rng(3)
determined = 0
for trial, (n, amp) in enumerate([(8, 0.15), (16, 0.02), (16, 0.02), (32, 0.005)]):
    mod = M.Model(A=Fr(1, 2), kappa=Fr(int(rng.integers(0, 5)), 4), K=Fr(int(rng.integers(-2, 3)), 4), h=Fr(int(rng.integers(0, 6)), 2))
    x = C.graded_axis(2.0, n, 1.5)
    Mn = C.radial_field(x, x, C.belavin_polyakov(0.6, 2.0), 2.0)
    Mn[1:-1, 1:-1] += amp * rng.standard_normal(Mn[1:-1, 1:-1].shape)
    ub, parts = C.upper_bound(mod, x, x, Mn)
    Q, enc = C.degree(x, x, Mn)
    Eq, Qq = quadrature(mod, x, x, Mn)
    determined += Q is not None
    report(f"grid {2 * n}, noise {amp}: E_quad <= upper; degree {Q} agrees", Eq <= ub and (Q is None or Q == round(Qq)) and abs(Qq - round(Qq)) < 1e-6,
           f"quad {Eq:.6f}  upper {ub:.6f}  (rel. excess {(ub - Eq) / abs(Eq):.1e})")
report("the degree is determined for the finer fields", determined >= 3)

print("4. falsification: stored lower bounds below certified competitor energies")
for js in sorted(glob.glob('results/zeeman_h*.json')):
    d = json.load(open(js))
    cert = M.Certificate.from_json(d)
    h = float(cert.model.h)
    L = 8.0 / np.sqrt(h)
    x = C.graded_axis(L, 60, 3.0)
    R = 1.0 / h
    Mn = C.radial_field(x, x, C.belavin_polyakov(R, L, power=4), L)
    ub, parts = C.upper_bound(cert.model, x, x, Mn)
    Q, enc = C.degree(x, x, Mn)
    ok = cert.check() and Q == -1 and float(d['lower']) <= ub
    report(f"{js}: lower {float(d['lower']):.6f} <= upper of a crude competitor", ok, f"upper {ub:.6f}")

print("\nALL TESTS PASSED" if fails == 0 else f"\nFAILURES: {fails}")
sys.exit(1 if fails else 0)
