#!/usr/bin/env python3
"""Independent cross-check of the stored rational-map certificates (docs/SKYRME_RATIONAL_MAPS.md).

Floating point, NumPy only, no code shared with the Wolfram Language model.  For every certificate in
certificates/ it recomputes, by a different method,

  1. the angular integral  II = (1/4 pi) int_{S^2} J^4 dA,  J = (1 + |z|^2) |R'(z)| / (1 + |R(z)|^2):
     Gauss-Legendre in the polar angle on the two hemispheres (the upper one in the coordinate 1/z) times the
     trapezoidal rule in the azimuth, and checks that the value lies in the certified enclosure;
  2. the degree, as (1/4 pi) int J^2 dA (it must be the integer B);
  3. the radial energy of the stored profile (piecewise cubic with the stored rational data, linear on the last
     cell) with the certified upper bound of II, by Gauss-Legendre on every cell, and checks that it lies below the
     certified upper bound and close to it.

These are floating-point evaluations: they cannot prove the certified numbers, but an error in the certified
evaluation of the size of the last digits quoted in the note would show here.

usage: python check_rational_maps.py
"""
import glob
import json
import math
import os
import sys
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def cpoly(coefs):
    """ascending complex coefficients from [[re, im], ...] with rational strings"""
    return np.array([complex(float(Fr(a)), float(Fr(b))) for a, b in coefs])


def pval(c, z):
    out = np.zeros_like(z, dtype=complex)
    for a in c[::-1]:
        out = out * z + a
    return out


def pder(c):
    return np.array([k * c[k] for k in range(1, len(c))] or [0j])


def angular(p, q, B, ntheta=400, nphi=1600):
    """(II, degree) by quadrature: both hemispheres in a coordinate in which the unit disk is the hemisphere"""
    x, w = np.polynomial.legendre.leggauss(ntheta)
    theta = 0.25 * math.pi * (x + 1.0)                 # 0 .. pi/2
    wt = 0.25 * math.pi * w
    phi = 2 * math.pi * (np.arange(nphi) + 0.5) / nphi
    T, P = np.meshgrid(theta, phi, indexing="ij")
    z = np.tan(T / 2) * np.exp(1j * P)                 # |z| <= 1
    tot4 = tot2 = 0.0
    pp = np.concatenate([p, np.zeros(B + 1 - len(p))])
    qq = np.concatenate([q, np.zeros(B + 1 - len(q))])
    for a, b in ((pp, qq), (pp[::-1], qq[::-1])):      # the map, and the map in the coordinate 1/z
        pa, qa = pval(a, z), pval(b, z)
        wr = pval(pder(a), z) * qa - pa * pval(pder(b), z)
        J = (1 + np.abs(z) ** 2) * np.abs(wr) / (np.abs(pa) ** 2 + np.abs(qa) ** 2)
        dA = np.sin(T) * wt[:, None] * (2 * math.pi / nphi)
        tot4 += np.sum(J ** 4 * dA)
        tot2 += np.sum(J ** 2 * dA)
    return tot4 / (4 * math.pi), tot2 / (4 * math.pi)


def radial(r, p, dp, B, II, mu2, ng=12):
    """4 pi int e dr / (12 pi^2) for the stored trial profile"""
    x, w = np.polynomial.legendre.leggauss(ng)
    n = len(r) - 1
    total = 0.0
    for k in range(n):
        a, b = r[k], r[k + 1]
        h = b - a
        s = 0.5 * (x + 1.0)
        rr = a + h * s
        if k < n - 1:
            y0, y1, m0, m1 = p[k], p[k + 1], dp[k], dp[k + 1]
            P = (2 * s ** 3 - 3 * s ** 2 + 1) * y0 + (s ** 3 - 2 * s ** 2 + s) * h * m0 + (-2 * s ** 3 + 3 * s ** 2) * y1 + (s ** 3 - s ** 2) * h * m1
            D = ((6 * s ** 2 - 6 * s) * (y0 - y1)) / h + (3 * s ** 2 - 4 * s + 1) * m0 + (3 * s ** 2 - 2 * s) * m1
        else:                                          # last cell: p rises linearly to pi
            P = p[k] + (math.pi - p[k]) * s
            D = np.full_like(s, (math.pi - p[k]) / h)
        e = rr ** 2 * D ** 2 + 2 * B * np.sin(P) ** 2 * (1 + D ** 2) + II * np.sin(P) ** 4 / rr ** 2 + 2 * mu2 * rr ** 2 * (1 + np.cos(P))
        total += 0.5 * h * np.sum(w * e)
    return total / (3 * math.pi)


def main():
    files = sorted(glob.glob(os.path.join(HERE, "certificates", "*.json")))
    if not files:
        print("no certificates found")
        return 1
    bad = 0
    for f in files:
        c = json.load(open(f, encoding="utf-8"))
        B = int(c["degree"])
        mu = Fr(c["mu"])
        p, q = cpoly(c["map"]["p"]), cpoly(c["map"]["q"])
        lo, hi = float(Fr(c["angular_integral"]["lower"])), float(Fr(c["angular_integral"]["upper"]))
        II, deg = angular(p, q, B)
        ok_i = lo * (1 - 1e-12) <= II <= hi * (1 + 1e-12)
        ok_d = abs(deg - B) < 1e-9
        t = c["trial"]
        r = np.array([float(Fr(v)) for v in t["r"]])
        pv = np.array([float(Fr(v)) for v in t["p"]])
        dpv = np.array([float(Fr(v)) for v in t["dp"]])
        ub = float(Fr(c["upper_bound"]))
        e = radial(r, pv, dpv, B, hi, float(mu) ** 2)
        ok_e = e <= ub and ub - e <= 2e-5 * ub
        bad += not (ok_i and ok_d and ok_e)
        print(f"{os.path.basename(f):22s} B={B} mu={mu}:  II = {II:.10f} in [{lo:.10f}, {hi:.10f}] {'ok' if ok_i else 'OUTSIDE'};  "
              f"degree {deg:.10f} {'ok' if ok_d else 'WRONG'};  energy/12pi^2 = {e:.9f} <= {ub:.9f} (difference {ub - e:.1e}) "
              f"{'ok' if ok_e else 'MISMATCH'};  per unit of degree {ub / B:.6f}", flush=True)
    print("ALL CHECKS PASSED" if bad == 0 else f"{bad} PROBLEM(S)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main())
