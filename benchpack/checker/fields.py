"""Certificate data -> exact polynomial fields on every cell, and the per-cell work of the checker.

Built from benchpack/certificates/FORMAT.md and benchpack/problems/definitions.py only.  Per cell (local coordinates
s, t in [0, 1]) every field is a polynomial with rational coefficients (the stored float64 numbers are taken as the
exact binary rationals they represent).
"""
import json
import os
import sys
from fractions import Fraction as Fr

import numpy as np

import exact as X
import region as RG

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'problems'))
import definitions as D      # noqa: E402  (the only file shared with the builder)

R_HOLE2 = D.P1['hole_radius'] ** 2                    # 1/4
R_PATCH_OUT2 = D.P1['patch']['r_max'] ** 2            # 121/400
assert D.P1['patch']['r_min'] == D.P1['hole_radius'] and D.P1['patch']['theta_half'] == Fr(1, 10)
PATCH_LOWERS = (('line', 'cot(1/10)'), ('circle', R_HOLE2))
PATCH_UPPERS = (('circle', R_PATCH_OUT2),)
# a rational box containing the half patch P_h: x <= (11/20) sin(1/10) < 11/200, y >= (1/2) cos(1/10) > 497/1000
PATCH_BOX = ((Fr(0), Fr(11, 200)), (Fr(497, 1000), Fr(11, 20)))


# ------------------------------------------------------------------------------------------------- loading
def load(stem):
    """read STEM.json and STEM.npz; returns (meta, arrays) with arrays as exact Fraction object arrays"""
    with open(stem + '.json', encoding='utf-8') as f:
        meta = json.load(f)
    npz = np.load(stem + '.npz')
    arrays = {k: X.farr(npz[k]) for k in npz.files}
    raw = {k: np.asarray(npz[k]) for k in npz.files}
    return meta, arrays, raw


class Space:
    def __init__(self, sp):
        self.p = int(sp['degree'])
        self.kx = [Fr(v) for v in sp['knots_x']]
        self.ky = [Fr(v) for v in sp['knots_y']]
        self.nx = len(self.kx) - self.p - 1
        self.ny = len(self.ky) - self.p - 1
        self.sx = X.bspline_pieces(self.kx, self.p)
        self.sy = X.bspline_pieces(self.ky, self.p)
        self.g_exp = [int(e) for e in sp.get('g_exponents', [])]
        self.g_base = Fr(sp['g_base']) if 'g_base' in sp else None

    def cells(self):
        return [(i, j) for i in range(len(self.sx)) for j in range(len(self.sy))]

    def cell(self, i, j):
        kx, x0, x1, Mx = self.sx[i]
        ky, y0, y1, My = self.sy[j]
        return dict(kx=kx, ky=ky, x0=x0, x1=x1, y0=y0, y1=y1, h=x1 - x0, k=y1 - y0, Mx=Mx, My=My)

    def local(self, A, c):
        p = self.p
        return X.tensor(c['Mx'], A[c['kx'] - p:c['kx'] + 1, c['ky'] - p:c['ky'] + 1], c['My'])


class Problem:
    """problem configuration: domain classification, boundary edges, the particular parts of the fields"""

    def __init__(self, meta):
        self.name = meta['problem']
        self.variant = meta.get('variant', 'standard')
        if self.name == 'P1':
            self.box = (Fr(0), Fr(1))
            self.hole = self.variant != 'nohole-test'       # test-only variant: the plate without the hole
            self.frictionless = self.variant in ('frictionless', 'nohole-test')
            self.weight = self.hole                           # stress function multiplied by rho^2
            self.work_factor = Fr(1, 2)                       # int_{grip} t.u = F_q * delta/2
            self.K_factor = D.P1_QUARTER['K_from_quarter_energy']
            self.patch_measure = D.P1_QUARTER['half_patch_area']
        elif self.name == 'P2':
            self.box = (Fr(-1), Fr(1))
            self.hole = False
            self.weight = False
            self.work_factor = Fr(1)
            self.K_factor = Fr(1)
            self.patch_measure = D.P2['patch_area']
        else:
            raise ValueError('unknown problem %r' % self.name)

    # -- domain
    def classify(self, c):
        """'full', 'cut' or 'out' (exact)"""
        x0, x1, y0, y1 = c['x0'], c['x1'], c['y0'], c['y1']
        if self.name == 'P1':
            if not self.hole:
                return 'full'
            assert x0 >= 0 and y0 >= 0
            if x1 * x1 + y1 * y1 <= R_HOLE2:
                return 'out'
            if x0 * x0 + y0 * y0 >= R_HOLE2:
                return 'full'
            return 'cut'
        # P2: Omega = [-1,1]^2 minus (0,1] x [-1,0); 0 must be a knot in both directions
        assert not (x0 < 0 < x1) and not (y0 < 0 < y1), 'a cell straddles x = 0 or y = 0'
        return 'out' if (x0 >= 0 and y1 <= 0) else 'full'

    def edge_in_domain(self, xa, ya, xb, yb):
        """True if the open segment (axis-parallel, xa <= xb, ya <= yb) meets the interior of the domain"""
        if self.name == 'P1':
            if not self.hole:
                return True
            return xb * xb + yb * yb > R_HOLE2
        # P2: the segment is in the closed L; it meets the interior unless it lies in the removed square's closure
        # minus... an edge of the grid is interior to Omega iff both cells beside it are in Omega (decided by caller)
        return True

    # -- particular parts
    def disp_base(self, c):
        """(vx0, vy0): the affine lift of the grip data, as local polynomials"""
        h, x0 = c['h'], c['x0']
        if self.name == 'P1':
            return [[x0 / 2], [h / 2]], [[Fr(0)]]                          # v_x = x/2
        return [[Fr(0)]], [[-(x0 + 1) / 2], [-h / 2]]                       # v_y = -(x+1)/2


def as_arr(lst):
    a = X.zeros((len(lst), len(lst[0])))
    for i, row in enumerate(lst):
        for j, v in enumerate(row):
            a[i, j] = Fr(v)
    return a


def local_poly_from_dict(coefs, x0, h, y0, k):
    """sum c_ij x^i y^j in local coordinates (x = x0 + h s, y = y0 + k t)"""
    di = max(i for i, j in coefs) + 1
    dj = max(j for i, j in coefs) + 1
    xs = [X.ppow([x0, h], i) for i in range(di)]
    ys = [X.ppow([y0, k], j) for j in range(dj)]
    P = X.zeros((di, dj))
    for (i, j), c in coefs.items():
        for a, u in enumerate(xs[i]):
            for b, w in enumerate(ys[j]):
                P[a, b] += c * u * w
    return P


# ------------------------------------------------------------------------------------------------- per-cell fields
def disp_cell(prob, sp, A, B, c):
    """local polynomials vx, vy, exx, eyy, gxy"""
    bx, by = prob.disp_base(c)
    vx = X.b_add(sp.local(A, c), as_arr(bx))
    vy = X.b_add(sp.local(B, c), as_arr(by))
    h, k = c['h'], c['k']
    exx = X.b_ds(vx) / h
    eyy = X.b_dt(vy) / k
    gxy = X.b_add(X.b_dt(vx) / k, X.b_ds(vy) / h)
    return vx, vy, exx, eyy, gxy


def chi_local(y0, k):
    """chi(y) of P2 on a cell (y0 >= 0 or y1 <= 0), as a polynomial in t"""
    if y0 + k <= 0:
        return [Fr(1)]
    assert y0 >= 0
    return X.padd(X.padd([Fr(1)], X.pscale(X.ppow([y0, k], 2), -3)), X.pscale(X.ppow([y0, k], 3), 2))


def stress_cell(prob, sp, C, par, c):
    """local polynomials txx, tyy, txy and the Airy function phi (before the g part)"""
    h, k, x0, y0 = c['h'], c['k'], c['x0'], c['y0']
    Psi = sp.local(C, c)
    if prob.weight:
        rho = local_poly_from_dict(D.RHO, x0, h, y0, k)
        phi = X.b_mul(X.b_mul(rho, rho), Psi)
    else:
        phi = Psi
    if prob.name == 'P2':
        alpha, beta, gamma = par[0], par[1], par[2]
        lin = as_arr([[alpha + beta * x0 + gamma * y0, gamma * k], [beta * h, 0]])
        phi = X.b_add(phi, X.b_conv_t(lin, chi_local(y0, k)))
    txx = X.b_dt(X.b_dt(phi)) / (k * k)
    tyy = X.b_ds(X.b_ds(phi)) / (h * h)
    txy = -X.b_ds(X.b_dt(phi)) / (h * k)
    if prob.name == 'P1' and sp.g_exp:
        base = sp.g_base
        assert not (y0 < base < y0 + k), 'a stress cell straddles y = g_base'
        if y0 >= base:
            g = []
            for e, gam in zip(sp.g_exp, par):
                g = X.padd(g, X.pscale(X.ppow([y0 - base, k], e), gam))
            col = X.zeros((1, len(g)))
            for b, v in enumerate(g):
                col[0, b] = v
            txx = X.b_add(txx, col)
    return txx, tyy, txy, phi


def g_force(sp, par):
    """sum_k gamma_k (1 - base)^{e_k + 1}/(e_k + 1) = int_base^1 g(y) dy (exact)"""
    one = Fr(1)
    return sum((gam * (one - sp.g_base) ** (e + 1) / (e + 1) for e, gam in zip(sp.g_exp, par)), Fr(0))


def divisible_by_monic_t2(Pst, Q):
    """exact test: Q (local array, degree 2 in t with constant t^2 coefficient q2 != 0, no s t^2 terms) divides P.
    Reduces P modulo Q in t over Q[s]; True iff the remainder is 0."""
    q2 = Q[0, 2]
    assert q2 != 0 and all(Q[a, 2] == 0 for a in range(1, Q.shape[0]))
    assert Q.shape[1] == 3
    P = [list(Pst[:, b]) for b in range(Pst.shape[1])]       # P[b] = coefficient polynomial (in s) of t^b
    q1 = list(Q[:, 1])
    q0 = list(Q[:, 0])
    for b in range(len(P) - 1, 1, -1):
        lead = X.ptrim(P[b])
        if not lead:
            continue
        f = X.pscale(lead, 1 / q2)                       # remove f * t^{b-2} * Q
        P[b] = []
        P[b - 1] = X.padd(P[b - 1], X.pscale(X.pmul(f, q1), -1))
        P[b - 2] = X.padd(P[b - 2], X.pscale(X.pmul(f, q0), -1))
    return X.is_zero_poly(P[0]) and (len(P) < 2 or X.is_zero_poly(P[1]))


RHO_FACTORS = ({(2, 0): Fr(1), (0, 2): Fr(1), (0, 0): -R_HOLE2},           # x^2 + y^2 - 1/4
               {(0, 0): Fr(7, 4), (0, 2): Fr(1), (2, 0): Fr(-1)})          # 7/4 + y^2 - x^2


def check_rho_factors():
    """exact: RHO of definitions.py equals the product of the two factors above"""
    P = {}
    for (i1, j1), a in RHO_FACTORS[0].items():
        for (i2, j2), b in RHO_FACTORS[1].items():
            P[(i1 + i2, j1 + j2)] = P.get((i1 + i2, j1 + j2), Fr(0)) + a * b
    P = {key: v for key, v in P.items() if v}
    return P == {key: v for key, v in D.RHO.items() if v}


# ------------------------------------------------------------------------------------------------- worker
_STATE = {}


def init_worker(meta, arrays):
    sys.setrecursionlimit(10000)
    _STATE['meta'] = meta
    _STATE['arrays'] = arrays
    _STATE['prob'] = Problem(meta)
    _STATE['dsp'] = Space(meta['displacement_space'])
    _STATE['ssp'] = Space(meta['stress_space'])


def iv_pair(I):
    return (RG.lo_q(I), RG.hi_q(I))


def work_disp(args):
    """per cell of a displacement field: energy, patch integral, edge traces"""
    name, cells = args
    prob, sp, arr = _STATE['prob'], _STATE['dsp'], _STATE['arrays']
    A, B = arr[name + '_A'], arr[name + '_B']
    Cm = D.C
    out = []
    for (i, j) in cells:
        c = sp.cell(i, j)
        cls = prob.classify(c)
        vx, vy, exx, eyy, gxy = disp_cell(prob, sp, A, B, c)
        rec = dict(cell=(i, j), cls=cls)
        # energy  int eps.C.eps  (eps = (exx, eyy, gxy))
        hk = c['h'] * c['k']
        if cls == 'full':
            w = {(0, 0): Cm[0][0], (1, 1): Cm[1][1], (0, 1): 2 * Cm[0][1], (2, 2): Cm[2][2]}
            rec['energy'] = hk * X.quad_form([exx, eyy, gxy], w)
        elif cls == 'cut':
            dens = X.b_add(X.b_mul(exx, exx) * Cm[0][0], X.b_mul(eyy, eyy) * Cm[1][1],
                           X.b_mul(exx, eyy) * (2 * Cm[0][1]), X.b_mul(gxy, gxy) * Cm[2][2])
            I = RG.integrate_region(dens, c['x0'], c['h'], c['y0'], c['k'], [('circle', R_HOLE2)], [])
            rec['energy_iv'] = iv_pair(I)
        # patch functional integrand (C eps)_xx = C11 exx + C12 eyy
        if cls != 'out':
            rec['patch'] = patch_integral(prob, c, X.b_add(exx * Cm[0][0], eyy * Cm[0][1]))
        rec['tr'] = dict(L=(X.b_at_s(vx, 0), X.b_at_s(vy, 0)), R=(X.b_at_s(vx, 1), X.b_at_s(vy, 1)),
                         B=(X.b_at_t(vx, 0), X.b_at_t(vy, 0)), T=(X.b_at_t(vx, 1), X.b_at_t(vy, 1)))
        out.append(rec)
    return out


def work_stress(args):
    """per cell of a stress field: equilibrium, hole divisibility, energy, patch integral, edge traces"""
    name, cells = args
    prob, sp, arr = _STATE['prob'], _STATE['ssp'], _STATE['arrays']
    Cc = arr[name + '_C']
    par = list(arr[name + '_par']) if (name + '_par') in arr else []
    S = D.S
    out = []
    for (i, j) in cells:
        c = sp.cell(i, j)
        cls = prob.classify(c)
        txx, tyy, txy, phi = stress_cell(prob, sp, Cc, par, c)
        h, k = c['h'], c['k']
        rec = dict(cell=(i, j), cls=cls)
        div1 = X.b_add(X.b_ds(txx) / h, X.b_dt(txy) / k)
        div2 = X.b_add(X.b_ds(txy) / h, X.b_dt(tyy) / k)
        rec['div0'] = all(v == 0 for v in div1.flat) and all(v == 0 for v in div2.flat)
        hk = h * k
        if cls == 'full':
            w = {(0, 0): S[0][0], (1, 1): S[1][1], (0, 1): 2 * S[0][1], (2, 2): S[2][2]}
            rec['energy'] = hk * X.quad_form([txx, tyy, txy], w)
        elif cls == 'cut':
            dens = X.b_add(X.b_mul(txx, txx) * S[0][0], X.b_mul(tyy, tyy) * S[1][1],
                           X.b_mul(txx, tyy) * (2 * S[0][1]), X.b_mul(txy, txy) * S[2][2])
            I = RG.integrate_region(dens, c['x0'], h, c['y0'], k, [('circle', R_HOLE2)], [])
            rec['energy_iv'] = iv_pair(I)
            # zero traction on the hole: tau.grad(rho) divisible by rho (both factors), exactly
            rho = local_poly_from_dict(D.RHO, c['x0'], h, c['y0'], k)
            rx, ry = X.b_ds(rho) / h, X.b_dt(rho) / k
            t1 = X.b_add(X.b_mul(txx, rx), X.b_mul(txy, ry))
            t2 = X.b_add(X.b_mul(txy, rx), X.b_mul(tyy, ry))
            facs = [local_poly_from_dict(f, c['x0'], h, c['y0'], k) for f in RHO_FACTORS]
            rec['hole_div'] = all(divisible_by_monic_t2(t, f) for t in (t1, t2) for f in facs)
        if cls != 'out':
            rec['patch'] = patch_integral(prob, c, txx)
        rec['tr'] = dict(L=(X.b_at_s(txx, 0), X.b_at_s(txy, 0)), R=(X.b_at_s(txx, 1), X.b_at_s(txy, 1)),
                         B=(X.b_at_t(txy, 0), X.b_at_t(tyy, 0)), T=(X.b_at_t(txy, 1), X.b_at_t(tyy, 1)))
        out.append(rec)
    return out


def patch_integral(prob, c, P):
    """int over (patch ∩ cell) of the local polynomial P: exact Fraction (P2) or interval pair (P1); None if the cell
    does not meet the patch"""
    x0, x1, y0, y1, h, k = c['x0'], c['x1'], c['y0'], c['y1'], c['h'], c['k']
    if prob.name == 'P2':
        (px0, px1), (py0, py1) = D.P2['patch']
        a0, a1, b0, b1 = max(x0, px0), min(x1, px1), max(y0, py0), min(y1, py1)
        if a0 >= a1 or b0 >= b1:
            return None
        return h * k * X.integral_rect(P, (a0 - x0) / h, (a1 - x0) / h, (b0 - y0) / k, (b1 - y0) / k)
    (bx0, bx1), (by0, by1) = PATCH_BOX
    if x0 >= bx1 or x1 <= bx0 or y0 >= by1 or y1 <= by0:
        return None
    I = RG.integrate_region(P, x0, h, y0, k, list(PATCH_LOWERS), list(PATCH_UPPERS))
    return iv_pair(I)
