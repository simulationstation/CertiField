#!/usr/bin/env python3
"""The certified upper bound of the magnetisation track: energy and degree of a normalised bilinear interpolant.
docs/MAGNETISATION.md, Section 5.

    m_h = M_h / |M_h| ,   M_h bilinear on every cell of a tensor-product grid with nodal vectors Mn (floats, taken as
    the exact binary numbers they are).

upper_bound and degree are RIGOROUS: outward-rounded interval arithmetic (ia.py) on exact inputs.  They refuse a
field for which |M_h|^2 > 0 cannot be shown on every cell.  The functions that build nodal values are heuristics.
"""
import math
from fractions import Fraction as Fr

import numpy as np

from ia import I, fr, dn, up

C = math.comb
# exact integrals of products of Bernstein polynomials on [0, 1]
G2 = [[Fr(C(2, i) * C(2, j), C(4, i + j) * 5) for j in range(3)] for i in range(3)]                    # b2_i b2_j
T2 = [[[Fr(C(2, i) * C(2, j) * C(2, l), C(6, i + j + l) * 7) for l in range(3)] for j in range(3)] for i in range(3)]
G12 = [[Fr(C(1, i) * C(2, k), C(3, i + k) * 4) for k in range(3)] for i in range(2)]                   # b1_i b2_k


def _dot(a, b, comps=(0, 1, 2)):
    s = None
    for c in comps:
        t = a[c] * b[c]
        s = t if s is None else s + t
    return s


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _half_sum(a, b):
    return [(x + y) * fr(Fr(1, 2)) for x, y in zip(a, b)]


def _sum_up(x):
    return float(up(up(math.fsum(np.asarray(x, dtype=np.float64).ravel().tolist()))))


def _sum_dn(x):
    return float(dn(dn(math.fsum(np.asarray(x, dtype=np.float64).ravel().tolist()))))


def _cells(x1, x2, Mn):
    x1, x2 = np.asarray(x1, dtype=np.float64), np.asarray(x2, dtype=np.float64)
    Mn = np.asarray(Mn, dtype=np.float64)
    assert Mn.shape == (x1.size, x2.size, 3) and np.all(np.diff(x1) > 0) and np.all(np.diff(x2) > 0) and np.all(np.isfinite(Mn))
    B = {(0, 0): [I(Mn[:-1, :-1, c]) for c in range(3)], (1, 0): [I(Mn[1:, :-1, c]) for c in range(3)],
         (0, 1): [I(Mn[:-1, 1:, c]) for c in range(3)], (1, 1): [I(Mn[1:, 1:, c]) for c in range(3)]}
    h1 = (I(x1[1:]) - I(x1[:-1]))[:, None]
    h2 = (I(x2[1:]) - I(x2[:-1]))[None, :]
    return B, I(h1.lo * np.ones((1, x2.size - 1)), h1.hi * np.ones((1, x2.size - 1))), I(h2.lo * np.ones((x1.size - 1, 1)), h2.hi * np.ones((x1.size - 1, 1)))


def _D_coefficients(B, comps=(0, 1, 2)):
    """the nine Bernstein coefficients c[k][l] (k for xi, l for eta) of sum over comps of M_c^2"""
    d = lambda p, q: _dot(B[p], B[q], comps)
    c = [[None] * 3 for _ in range(3)]
    c[0][0], c[2][0], c[0][2], c[2][2] = d((0, 0), (0, 0)), d((1, 0), (1, 0)), d((0, 1), (0, 1)), d((1, 1), (1, 1))
    c[1][0], c[0][1], c[2][1], c[1][2] = d((0, 0), (1, 0)), d((0, 0), (0, 1)), d((1, 0), (1, 1)), d((0, 1), (1, 1))
    c[1][1] = (d((0, 0), (1, 1)) + d((1, 0), (0, 1))) * fr(Fr(1, 2))
    return c


def _range(coefs):
    lo = np.minimum.reduce([x.lo * np.ones(coefs[0].shape) for x in coefs])
    hi = np.maximum.reduce([x.hi * np.ones(coefs[0].shape) for x in coefs])
    return lo, hi


def upper_bound(model, x1, x2, Mn, dmi_axes=((1, 0, 0), (0, 1, 0))):
    """RIGOROUS upper bound for E[m_h], m_h = M_h/|M_h|.  Returns (value, parts) or raises ValueError if |M_h| > 0
    cannot be shown.  model: magnet.Model with h >= 0."""
    assert model.h >= 0
    B, h1, h2 = _cells(x1, x2, Mn)
    c = _D_coefficients(B)
    flat = [c[k][l] for k in range(3) for l in range(3)]
    a_lo, b_hi = _range(flat)
    if not np.all(a_lo > 0):
        raise ValueError("the interpolant may vanish: |M_h|^2 > 0 could not be shown on every cell")
    a, b = I(a_lo), I(b_hi)                              # D takes values in [a, b] on the cell
    shape = a_lo.shape
    third = fr(Fr(1, 3))
    Cl = [(c[0][l] + c[1][l] + c[2][l]) * third for l in range(3)]          # integral over xi of D: coefficients in eta
    Ck = [(c[k][0] + c[k][1] + c[k][2]) * third for k in range(3)]          # integral over eta of D: coefficients in xi
    # ---- exchange: A [ (h2/h1) int X(eta)/D^2 + (h1/h2) int Y(xi)/D^2 ],   1/D^2 <= alpha2 - beta2 D
    v = [_cross(B[(0, 0)], B[(1, 0)]), _half_sum(_cross(B[(0, 0)], B[(1, 1)]), _cross(B[(0, 1)], B[(1, 0)])), _cross(B[(0, 1)], B[(1, 1)])]
    w = [_cross(B[(0, 0)], B[(0, 1)]), _half_sum(_cross(B[(0, 0)], B[(1, 1)]), _cross(B[(1, 0)], B[(0, 1)])), _cross(B[(1, 0)], B[(1, 1)])]
    a2b2 = a * a * b * b
    alpha2 = (a * a + a * b + b * b) * a2b2.inv()
    beta2 = (a + b) * a2b2.inv()

    def quartic(vv, Cdir):
        """int |vv|^2 and int int |vv|^2 D for the quadratic vector polynomial vv (Bernstein coefficients)"""
        s0 = s1 = None
        for i in range(3):
            for j in range(i, 3):
                dij = _dot(vv[i], vv[j]) * fr(Fr(1 if i == j else 2))
                t0 = dij * fr(G2[i][j])
                s0 = t0 if s0 is None else s0 + t0
                for l in range(3):
                    t1 = dij * Cdir[l] * fr(T2[i][j][l])
                    s1 = t1 if s1 is None else s1 + t1
        return s0, s1

    X0, X1 = quartic(v, Cl)
    Y0, Y1 = quartic(w, Ck)
    ex = ((alpha2 * X0 - beta2 * X1) * h2 * h1.inv() + (alpha2 * Y0 - beta2 * Y1) * h1 * h2.inv()) * fr(model.A)
    total = ex
    parts = dict(exchange=_sum_up(ex.hi))
    # ---- Dzyaloshinskii-Moriya: -kappa [ h2 int a1.v(eta)/D + h1 int a2.w(xi)/D ],  1/D = T1(D) + [0, g]
    if model.kappa != 0:
        d0 = 0.5 * (a_lo + b_hi)
        d0i = I(d0)
        g1 = (a - d0i).sqr() * (a * d0i * d0i).inv()
        g2 = (b - d0i).sqr() * (b * d0i * d0i).inv()
        gap = I(np.zeros(shape), np.maximum(g1.hi, g2.hi))

        def linear(vv, axis, Cdir):
            n = [_dot(vv[i], [fr(Fr(x)) for x in axis]) for i in range(3)]
            n0 = (n[0] + n[1] + n[2]) * third
            n1 = None
            for i in range(3):
                for l in range(3):
                    t = n[i] * Cdir[l] * fr(G2[i][l])
                    n1 = t if n1 is None else n1 + t
            nmax = np.maximum.reduce([np.maximum(np.abs(x.lo), np.abs(x.hi)) * np.ones(shape) for x in n])
            main = n0 * fr(2) * d0i.inv() - n1 * (d0i * d0i).inv()
            return main + I(-nmax, nmax) * gap

        dm = (linear(v, dmi_axes[0], Cl) * h2 + linear(w, dmi_axes[1], Ck) * h1) * fr(-model.kappa)
        total = total + dm
        parts['dmi'] = _sum_up(dm.hi)
    area = h1 * h2
    # ---- anisotropy: K (M1^2 + M2^2)/D
    if model.K != 0:
        q = _D_coefficients(B, (0, 1))
        q0 = None
        q1 = None
        for k in range(3):
            for l in range(3):
                t = q[k][l]
                q0 = t if q0 is None else q0 + t
                for k2 in range(3):
                    for l2 in range(3):
                        t1 = q[k][l] * c[k2][l2] * fr(G2[k][k2] * G2[l][l2])
                        q1 = t1 if q1 is None else q1 + t1
        q0 = q0 * fr(Fr(1, 9))
        if model.K > 0:
            an = (q0 * (a + b) - q1) * (a * b).inv() * area * fr(model.K)
        else:
            d0i = I(0.5 * (a_lo + b_hi))
            an = (q0 * fr(2) * d0i.inv() - q1 * (d0i * d0i).inv()) * area * fr(model.K)
        total = total + an
        parts['anisotropy'] = _sum_up(an.hi)
    # ---- Zeeman: h (1 - M3/sqrt(D)),   (3 - D)/2 <= 1/sqrt(D) <= (1 + (a + b - D)/(ab))/2
    if model.h != 0:
        m3 = {k: B[k][2] for k in B}
        m0 = (m3[(0, 0)] + m3[(1, 0)] + m3[(0, 1)] + m3[(1, 1)]) * fr(Fr(1, 4))
        m1 = None
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    for l in range(3):
                        t = m3[(i, j)] * c[k][l] * fr(G12[i][k] * G12[j][l])
                        m1 = t if m1 is None else m1 + t
        m3min = np.minimum.reduce([m3[k].lo for k in m3])
        neg = I(np.maximum(-m3min, 0.0))
        one = I(np.ones(shape))
        gz = I(np.maximum(((one - a).sqr() * (a * fr(2)).inv()).hi, ((one - b).sqr() * (b * fr(2)).inv()).hi))
        ze = (one - m0 * fr(Fr(3, 2)) + m1 * fr(Fr(1, 2)) + neg * gz) * area * fr(model.h)
        total = total + ze
        parts['zeeman'] = _sum_up(ze.hi)
    parts['min_D'] = float(a_lo.min())
    parts['max_D'] = float(b_hi.max())
    return _sum_up(total.hi), parts


def degree(x1, x2, Mn):
    """RIGOROUS: the degree of m_h = M_h/|M_h| (constant on the boundary), or None if the enclosure of
    int m.(d1 m x d2 m) does not single out one integer.  Returns (Q, (lo, hi)) with the enclosure of the integral."""
    B, h1, h2 = _cells(x1, x2, Mn)
    c = _D_coefficients(B)
    a_lo, b_hi = _range([c[k][l] for k in range(3) for l in range(3)])
    if not np.all(a_lo > 0):
        return None, (float('-inf'), float('inf'))
    a, b = I(a_lo), I(b_hi)

    def det(p, q, r):
        return _dot(B[p], _cross(B[q], B[r]))

    corners = [det((0, 0), (1, 0), (0, 1)), det((0, 0), (1, 0), (1, 1)), det((0, 0), (1, 1), (0, 1)), det((1, 0), (1, 1), (0, 1))]
    dlo, dhi = _range(corners)
    wlo = (fr(2) * (b * (b + fr(1))).inv()).lo           # D^(-3/2) >= 2/(b (1 + b))
    whi = ((a + fr(1)) * (a * a * fr(2)).inv()).hi       # D^(-3/2) <= (1 + a)/(2 a^2)
    cell = I(dlo, dhi) * I(wlo, whi)
    lo, hi = _sum_dn(cell.lo), _sum_up(cell.hi)
    fourpi_lo, fourpi_hi = 4 * Fr(314159265358979323846264338327950288, 10 ** 35), 4 * Fr(314159265358979323846264338327950289, 10 ** 35)
    cands = [n for n in range(math.floor(lo / 12.0) - 2, math.ceil(hi / 12.0) + 3)
             if not (max(fourpi_lo * n, fourpi_hi * n) < Fr(lo) or min(fourpi_lo * n, fourpi_hi * n) > Fr(hi))]
    return (cands[0] if len(cands) == 1 else None), (lo, hi)


# ================================================================================================================
# heuristics: grids and nodal values
# ================================================================================================================
def graded_axis(L, n, beta=3.0):
    """2n + 1 nodes on [-L, L], symmetric, finer near 0 (ratio of the largest to the smallest cell about cosh(beta))"""
    tau = np.arange(n + 1) / n
    half = L * np.sinh(beta * tau) / np.sinh(beta) if beta > 0 else L * tau
    half[0], half[-1] = 0.0, L
    return np.concatenate([-half[:0:-1], half])


def radial_field(x1, x2, theta, L, helicity=np.pi / 2):
    """nodal vectors of the axisymmetric field (sin(theta) e^{i(psi + helicity)}, cos(theta)), theta = theta(r) a
    function with theta(r) = 0 for r >= L; the nodes with r >= L and the boundary nodes are exactly (0, 0, 1)"""
    X, Y = np.meshgrid(x1, x2, indexing='ij')
    r = np.hypot(X, Y)
    psi = np.arctan2(Y, X)
    th = theta(r)
    Mn = np.stack([np.sin(th) * np.cos(psi + helicity), np.sin(th) * np.sin(psi + helicity), np.cos(th)], axis=-1)
    out = (r >= L)
    out[0, :] = out[-1, :] = out[:, 0] = out[:, -1] = True
    Mn[out] = (0.0, 0.0, 1.0)
    return Mn


def belavin_polyakov(R, L, power=2):
    """theta(r) of a Belavin-Polyakov map of radius R, cut off so that theta = 0 for r >= L"""
    def theta(r):
        chi = np.clip(1 - (r / L) ** 2, 0.0, None) ** power
        with np.errstate(divide='ignore'):
            return 2 * np.arctan2(R * chi, r)
    return theta
