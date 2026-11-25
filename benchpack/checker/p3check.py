"""Independent checker of the P3 certificates (benchpack/certificates/FORMAT.md section 9; octant model, K only).

usage:  python p3check.py <path/to/STEM> [--procs 2] [--log FILE]

Built from FORMAT.md section 9 and definitions.py only (no code of the builder).  On every cell of the uniform grid
(side 1/N) the stress tau = g(y,z) e_x e_x + inc(rho^2 Psi) and the displacement v = (x/2) e_x + spline are rebuilt
as polynomials with rational coefficients (stored floats are exact binary rationals), and in exact arithmetic:

  equilibrium in every cell of the domain; continuity of the traction tau_kj across every face x_k = const that meets
  the domain; on every cell cut by the sphere, each component of tau grad(rho) divisible by rho (zero traction on the
  sphere); zero traction on y = 1 and z = 1; zero shear on the mirror planes and (frictionless) on the grip; the grip
  force F_o integrated from the field and compared with sum (gamma_k + delta_k) (1/2)^(e_k+1)/(e_k+1); the
  displacement continuous across every face, with the mirror and grip values; the coefficient conditions; the knots.

Integrals: superset mode (JSON integration.mode = "superset"): the octree rule of FORMAT.md section 9 is applied
exactly (integer tests), and every kept box is a dyadic sub-box of a cell, so Q_o and a_o are exact rationals.
Bounds: 2 F_o^2/Q_o <= K <= 8 a_o.  Work is done slab by slab (x index) in two processes with a checkpoint after
every slab (cache/), so an interrupted run resumes.
"""
import argparse
import hashlib
import json
import multiprocessing as mp_
import os
import pickle
import sys
import time
from fractions import Fraction as Fr

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '2')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np                       # noqa: E402

import exact as X                        # noqa: E402
import exact3 as E                       # noqa: E402
import fields as FD                      # noqa: E402  (imports definitions.py, the shared problem file)
import region as RG                      # noqa: E402
from mpmath import iv                    # noqa: E402

D = FD.D
COMP = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]     # Psi component order (xx, yy, zz, xy, yz, xz)
HALF = Fr(1, 2)
R2 = D.P3['cavity_radius'] ** 2                            # 1/4


# ------------------------------------------------------------------------------------------------- data
def load(stem):
    with open(stem + '.json', encoding='utf-8') as f:
        meta = json.load(f)
    z = np.load(stem + '.npz')
    arrays = {k: X.farr(z[k]) for k in z.files}
    return meta, arrays


class Axis:
    """B-spline pieces of one axis in the local coordinate of each cell [edge_e, edge_(e+1)]"""

    def __init__(self, knots, p, edges):
        self.p = p
        self.knots = [Fr(v) for v in knots]
        self.ok, self.mult, self.msg = X.knot_report(self.knots, p, Fr(0), Fr(1))
        self.on_grid = set(self.knots) <= set(edges)
        self.n = len(self.knots) - p - 1
        spans = X.bspline_pieces(self.knots, p)
        self.cells = []
        for e in range(len(edges) - 1):
            x0, x1 = edges[e], edges[e + 1]
            sp = [s for s in spans if s[1] <= x0 and x1 <= s[2]]
            assert len(sp) == 1, 'a grid cell is not inside one knot span'
            k, t0, t1, M = sp[0]
            h = t1 - t0
            Mc = X.zeros((p + 1, p + 1))
            for i in range(p + 1):
                row = X.pcompose_lin(list(M[i, :]), (x0 - t0) / h, (x1 - x0) / h)
                for a in range(p + 1):
                    Mc[i, a] = row[a] if a < len(row) else Fr(0)
            self.cells.append((k - p, Mc))


class Cert:
    def __init__(self, meta, arrays):
        self.meta, self.arrays = meta, arrays
        self.N = int(meta['N'])
        self.variant = meta.get('variant')
        self.frictionless = self.variant in ('frictionless', 'nocavity-test')
        self.cavity = not str(self.variant).startswith('nocavity')      # test-only variants: the solid cube
        ds, ss = meta['displacement_space'], meta['stress_space']
        self.pu, self.ps = int(ds['degree']), int(ss['degree'])
        # cells: products of the knot spans (FORMAT.md 9: the knot vectors are authoritative)
        self.edges = [sorted(set(Fr(v) for v in ds['knots'][d]) | set(Fr(v) for v in ss['knots'][d])) for d in range(3)]
        self.n = [len(ed) - 1 for ed in self.edges]
        self.grade = meta.get('grade')
        self.dax = [Axis(k, self.pu, self.edges[d]) for d, k in enumerate(ds['knots'])]
        self.sax = [Axis(k, self.ps, self.edges[d]) for d, k in enumerate(ss['knots'])]
        self.gexp = [int(e) for e in ss.get('g_exponents', [])]
        self.gbase = Fr(ss.get('g_base', '1/2'))
        self.etab = Fr(ss.get('eta_flat_from', '3/4'))
        integ = meta.get('integration', {})
        self.mode = integ.get('mode')
        self.L = int(integ.get('stress_depth', 0))
        self.Lu = int(integ.get('displacement_depth', 0))
        lame = D.P3_LAME
        self.lam, self.mu, self.nu = lame['lam'], lame['mu'], D.P3['nu']
        self.V, self.Psi = arrays['V'], arrays['Psi']
        g = list(arrays['gam']) if 'gam' in arrays else []
        ng = len(self.gexp)
        self.gamma, self.delta = g[:ng], g[ng:]


def cell_class(cert, e):
    """'out' (inside the closed ball), 'full' (outside the open ball) or 'cut', exactly (level 0 of the rule)"""
    if not cert.cavity:
        return 'full'
    x0, h = corner(cert, e)
    rmin = sum(v * v for v in x0)
    rmax = sum((a + b) ** 2 for a, b in zip(x0, h))
    if rmax <= R2:
        return 'out'
    if rmin >= R2:
        return 'full'
    return 'cut'


def kept_boxes(cert, e, L):
    """the boxes of the octree rule of FORMAT.md 9 kept for cell e: (level, jx, jy, jz)"""
    x0, h = corner(cert, e)
    out = []

    def rec(l, j):
        lo = [x0[k] + h[k] * Fr(j[k], 2 ** l) for k in range(3)]
        hi = [lo[k] + h[k] / 2 ** l for k in range(3)]
        rmin = sum(v * v for v in lo)
        rmax = sum(v * v for v in hi)
        if rmax <= R2:
            return
        if rmin >= R2 or l == L:
            out.append((l,) + tuple(j))
            return
        for d in range(8):
            rec(l + 1, tuple(2 * j[k] + ((d >> k) & 1) for k in range(3)))
    if not cert.cavity:
        return [(0, 0, 0, 0)]
    rec(0, (0, 0, 0))
    return out


# ------------------------------------------------------------------------------------------------- fields on a cell
def corner(cert, e):
    """lower corner and side lengths of cell e"""
    ed = cert.edges
    return tuple(ed[d][e[d]] for d in range(3)), tuple(ed[d][e[d] + 1] - ed[d][e[d]] for d in range(3))


def rho_local(cert, e):
    x0, h = corner(cert, e)
    b = cert.etab
    if x0[0] + h[0] <= b:
        coefs = {(2, 0, 0): Fr(1), (0, 2, 0): Fr(1), (0, 0, 2): Fr(1), (0, 0, 0): Fr(-1, 4)}
    else:
        assert x0[0] >= b, 'a stress cell straddles x = eta_flat_from'
        # eta = x^2 - 4 (x - b)^2 = -3 x^2 + 8 b x - 4 b^2
        coefs = {(2, 0, 0): Fr(-3), (1, 0, 0): 8 * b, (0, 2, 0): Fr(1), (0, 0, 2): Fr(1), (0, 0, 0): -4 * b * b - Fr(1, 4)}
    return E.local_from_dict(coefs, x0, h)


def spline3(axes, A, e):
    (ox, Mx), (oy, My), (oz, Mz) = axes[0].cells[e[0]], axes[1].cells[e[1]], axes[2].cells[e[2]]
    p = Mx.shape[0] - 1
    return E.tensor3(Mx, My, Mz, A[ox:ox + p + 1, oy:oy + p + 1, oz:oz + p + 1])


def stress_cell(cert, e, weight=True):
    """tau as dict (i, j) -> local array (physical units), and rho (local)"""
    rho = rho_local(cert, e)
    R2 = E.t_mul(rho, rho) if (weight and cert.cavity) else None
    Phi = {}
    for c, (l, n) in enumerate(COMP):
        P = spline3(cert.sax, cert.Psi[c], e)
        if R2 is not None:
            P = E.t_mul(R2, P)
        Phi[(l, n)] = P
        Phi[(n, l)] = P
    x0, h = corner(cert, e)
    sc = tuple(1 / hh for hh in h)
    tau = E.inc(Phi, sc)
    gx = None
    for ax, coefs in ((1, cert.gamma), (2, cert.delta)):
        lo, hi = x0[ax], x0[ax] + h[ax]
        assert not (lo < cert.gbase < hi), 'a stress cell straddles the base point of g'
        if lo >= cert.gbase and coefs:
            g = []
            for ek, gk in zip(cert.gexp, coefs):
                g = X.padd(g, X.pscale(X.ppow([lo - cert.gbase, h[ax]], ek), gk))
            shp = [1, 1, 1]
            shp[ax] = len(g)
            arr = E.zeros3(tuple(shp))
            for b, v in enumerate(g):
                idx = [0, 0, 0]
                idx[ax] = b
                arr[tuple(idx)] = v
            gx = arr if gx is None else E.t_add(gx, arr)
    if gx is not None:
        tau[(0, 0)] = E.t_add(tau[(0, 0)], gx)
    return tau, rho


def disp_cell(cert, e):
    """v components (local) and strains [exx, eyy, ezz, gxy, gyz, gxz] (physical)"""
    x0, h = corner(cert, e)
    v = [spline3(cert.dax, cert.V[a], e) for a in range(3)]
    lift = E.zeros3((2, 1, 1))
    lift[0, 0, 0], lift[1, 0, 0] = x0[0] / 2, h[0] / 2
    vx = E.t_add(v[0], lift)
    vfull = [vx, v[1], v[2]]
    d = [[E.t_diff(vfull[a], k) / h[k] for k in range(3)] for a in range(3)]
    strains = [d[0][0], d[1][1], d[2][2], E.t_add(d[0][1], d[1][0]), E.t_add(d[1][2], d[2][1]),
               E.t_add(d[0][2], d[2][0])]
    return v, strains


def box_integrals(f, boxes, den=None):
    """exact sum over the boxes (level, j) of int_box f over the unit cell; f is a local Fraction array, or an integer
    array when den is given (the density is then f / den)"""
    if den is None:
        (fi,), Dn = E.to_int(f)
    else:
        fi, Dn = f, den
    tot = Fr(0)
    cache = {}
    for (l, jx, jy, jz) in boxes:
        vecs = []
        den = 1
        for d, j in enumerate((jx, jy, jz)):
            key = (l, j, f.shape[d])
            if key not in cache:
                m = E.moments(Fr(j, 2 ** l), Fr(j + 1, 2 ** l), f.shape[d])
                L = 1
                for v in m:
                    L = X.lcm(L, v.denominator)
                cache[key] = (np.array([v.numerator * (L // v.denominator) for v in m], dtype=object), L)
            vec, L = cache[key]
            vecs.append(vec)
            den *= L
        val = np.tensordot(np.tensordot(np.tensordot(fi, vecs[2], axes=([2], [0])), vecs[1], axes=([1], [0])),
                           vecs[0], axes=([0], [0]))
        tot += Fr(int(val), den)
    return tot / Dn


def cut_exact(cert, e, polys, weights):
    """(lo, hi) enclosing the integral of sum w_pq P_p P_q over (cell e minus the ball): the exact integral over the
    cell minus the closed-form integral over cell cap ball (ball.py)"""
    import ball
    x0, h = corner(cert, e)
    full = E.quad_form3(polys, weights) * (h[0] * h[1] * h[2])
    fi, den = E.density_int(polys, weights)
    g = ball.to_global(fi, x0, h)
    fglob = np.empty(g.shape, dtype=object)
    for idx, v in np.ndenumerate(g):
        fglob[idx] = Fr(v) / den                 # exact (v is an int or a Fraction)
    box = tuple((x0[d], x0[d] + h[d]) for d in range(3))
    I = ball.ball_part(fglob, box, D.P3['cavity_radius'] ** 2)
    return full - RG.hi_q(I), full - RG.lo_q(I)


# ------------------------------------------------------------------------------------------------- slabs
_S = {}


def init(meta, arrays, weight=True):
    FD.assert_independent()
    _S['cert'] = Cert(meta, arrays)
    _S['weight'] = weight


def face_in_domain(cert, ax, val, e):
    """the face x_ax = val of cell e meets the open domain (its farthest corner lies outside the closed ball)"""
    if not cert.cavity:
        return True
    x0, h = corner(cert, e)
    far = [x0[k] + h[k] for k in range(3)]
    far[ax] = val
    return sum(v * v for v in far) > R2


def new_state(a):
    return dict(next=a, Q=Fr(0), A=Fr(0), Qhi=Fr(0), Ahi=Fr(0), volhi=Fr(0), F=Fr(0), cells=0, cut=0, div_bad=0, face=0, face_bad=0, sphere_bad=0,
                bc_bad=0, vface=0, vface_bad=0, vbc_bad=0, fails=[], up=None, vup=None, first_low=None,
                vfirst_low=None, vol=Fr(0))


def note(st, msg):
    if len(st['fails']) < 20:
        st['fails'].append(msg)


def run_range(args):
    a, b, ckpt = args
    cert, weight = _S['cert'], _S['weight']
    nx, ny, nz = cert.n
    ed = cert.edges
    z = E.t_is_zero
    one, zero = Fr(1), Fr(0)
    S_w = {(0, 0): Fr(1), (1, 1): Fr(1), (2, 2): Fr(1), (0, 1): -2 * cert.nu, (0, 2): -2 * cert.nu,
           (1, 2): -2 * cert.nu, (3, 3): 2 * (1 + cert.nu), (4, 4): 2 * (1 + cert.nu), (5, 5): 2 * (1 + cert.nu)}
    C_w = {(0, 0): cert.lam + 2 * cert.mu, (1, 1): cert.lam + 2 * cert.mu, (2, 2): cert.lam + 2 * cert.mu,
           (0, 1): 2 * cert.lam, (0, 2): 2 * cert.lam, (1, 2): 2 * cert.lam,
           (3, 3): cert.mu, (4, 4): cert.mu, (5, 5): cert.mu}
    if ckpt and os.path.exists(ckpt):
        with open(ckpt, 'rb') as f:
            st = pickle.load(f)
    else:
        st = new_state(a)
    for i in range(st['next'], b):
        t0 = time.time()
        low, up, vlow, vup = {}, {}, {}, {}
        pend, vpend = {}, {}
        for j in range(ny):
            for k in range(nz):
                e = (i, j, k)
                cls = cell_class(cert, e)
                if cls == 'out':
                    continue
                st['cells'] += 1
                x0, h = corner(cert, e)
                vol = h[0] * h[1] * h[2]
                tau, rho = stress_cell(cert, e, weight)
                T = [tau[c] for c in COMP]
                # equilibrium
                for r in range(3):                       # physical divergence: d/dx_c = (1/h_c) d/ds_c
                    dv = E.t_add(*[E.t_diff(tau[(r, c)], c) / h[c] for c in range(3)])
                    if not z(dv):
                        st['div_bad'] += 1
                        note(st, 'div tau != 0 in cell %s' % (e,))
                # sphere
                if cls == 'cut':
                    st['cut'] += 1
                    grad = [E.t_diff(rho, c) / h[c] for c in range(3)]
                    for r in range(3):
                        g = E.t_add(*[E.t_mul(tau[(r, c)], grad[c]) for c in range(3)])
                        if not E.divisible_quadratic(g, rho, 0):
                            st['sphere_bad'] += 1
                            note(st, 'tau.grad(rho) not divisible by rho in cell %s' % (e,))
                # boundary faces
                bcs = []
                if j == ny - 1:
                    bcs += [(1, 1, r, 'traction on y = 1') for r in range(3)]
                if k == nz - 1:
                    bcs += [(2, 1, r, 'traction on z = 1') for r in range(3)]
                if i == 0 and face_in_domain(cert, 0, zero, e):
                    bcs += [(0, 0, 1, 'shear on x = 0'), (0, 0, 2, 'shear on x = 0')]
                if j == 0 and face_in_domain(cert, 1, zero, e):
                    bcs += [(1, 0, 0, 'shear on y = 0'), (1, 0, 2, 'shear on y = 0')]
                if k == 0 and face_in_domain(cert, 2, zero, e):
                    bcs += [(2, 0, 0, 'shear on z = 0'), (2, 0, 1, 'shear on z = 0')]
                if i == nx - 1 and cert.frictionless:
                    bcs += [(0, 1, 1, 'shear on the grip'), (0, 1, 2, 'shear on the grip')]
                for ax, side, r, what in bcs:
                    if not z(E.t_face(tau[(ax, r)], ax, side)):
                        st['bc_bad'] += 1
                        note(st, '%s in cell %s' % (what, (e,)))
                # force through the grip
                if i == nx - 1:
                    st['F'] += E.t_integral(E.t_face(tau[(0, 0)], 0, 1)) * (h[1] * h[2])
                # traction across faces: y and z inside the slab, x kept for the neighbouring slab
                for ax in (1, 2):
                    key = (e, ax)
                    if key in pend:
                        mine = [E.t_face(tau[(ax, r)], ax, 0) for r in range(3)]
                        st['face'] += 1
                        if not all(E.t_equal(u, w) for u, w in zip(pend.pop(key), mine)):
                            st['face_bad'] += 1
                            note(st, 'traction jump at the lower %s face of cell %s' % ('xyz'[ax], (e,)))
                    if e[ax] < cert.n[ax] - 1 and face_in_domain(cert, ax, ed[ax][e[ax] + 1], e):
                        nb = list(e)
                        nb[ax] += 1
                        if cell_class(cert, tuple(nb)) != 'out':
                            pend[(tuple(nb), ax)] = [E.t_face(tau[(ax, r)], ax, 1) for r in range(3)]
                low[(j, k)] = [E.t_face(tau[(0, r)], 0, 0) for r in range(3)]
                if i < nx - 1 and face_in_domain(cert, 0, ed[0][i + 1], e):
                    up[(j, k)] = [E.t_face(tau[(0, r)], 0, 1) for r in range(3)]
                # complementary energy: exact on whole cells; cut cells over the superset boxes or exactly
                if cls == 'full':
                    q = E.quad_form3(T, S_w) * vol
                    st['Q'] += q
                    st['Qhi'] += q
                    st['vol'] += vol
                    st['volhi'] += vol
                elif cert.mode == 'exact':
                    lo, hi = cut_exact(cert, e, T, S_w)
                    st['Q'] += lo
                    st['Qhi'] += hi
                    one3 = E.zeros3((1, 1, 1)) + 1
                    vlo, vhi = cut_exact(cert, e, [one3], {(0, 0): Fr(1)})
                    st['vol'] += vlo
                    st['volhi'] += vhi
                else:
                    boxes = kept_boxes(cert, e, cert.L)
                    if boxes:
                        fi, den = E.density_int(T, S_w)
                        q = box_integrals(fi, boxes, den) * vol
                        st['Q'] += q
                        st['Qhi'] += q
                        vb_ = sum((Fr(1, 8 ** bx[0]) for bx in boxes), Fr(0)) * vol
                        st['vol'] += vb_
                        st['volhi'] += vb_
                # displacement
                v, eps = disp_cell(cert, e)
                if cls == 'full':
                    q = E.quad_form3(eps, C_w) * vol
                    st['A'] += q
                    st['Ahi'] += q
                elif cert.mode == 'exact':
                    lo, hi = cut_exact(cert, e, eps, C_w)
                    st['A'] += lo
                    st['Ahi'] += hi
                else:
                    boxes = kept_boxes(cert, e, cert.Lu)
                    if boxes:
                        fi, den = E.density_int(eps, C_w)
                        q = box_integrals(fi, boxes, den) * vol
                        st['A'] += q
                        st['Ahi'] += q
                vb = []
                if i == 0 and face_in_domain(cert, 0, zero, e):
                    vb.append((0, 0, 0, 'v_x = 0 on x = 0'))
                if i == nx - 1:
                    vb.append((0, 0, 1, 'v_x = 1/2 on x = 1'))
                    if not cert.frictionless:
                        vb += [(1, 0, 1, 'v_y = 0 on the clamped grip'), (2, 0, 1, 'v_z = 0 on the clamped grip')]
                if j == 0 and face_in_domain(cert, 1, zero, e):
                    vb.append((1, 1, 0, 'v_y = 0 on y = 0'))
                if k == 0 and face_in_domain(cert, 2, zero, e):
                    vb.append((2, 2, 0, 'v_z = 0 on z = 0'))
                for comp, ax, side, what in vb:          # the spline parts vanish there (the lift x/2 gives 1/2)
                    if not z(E.t_face(v[comp], ax, side)):
                        st['vbc_bad'] += 1
                        note(st, '%s fails in cell %s' % (what, (e,)))
                for ax in (1, 2):
                    key = (e, ax)
                    if key in vpend:
                        mine = [E.t_face(v[c], ax, 0) for c in range(3)]
                        st['vface'] += 1
                        if not all(E.t_equal(u, w) for u, w in zip(vpend.pop(key), mine)):
                            st['vface_bad'] += 1
                            note(st, 'displacement jump at the lower %s face of cell %s' % ('xyz'[ax], (e,)))
                    if e[ax] < cert.n[ax] - 1 and face_in_domain(cert, ax, ed[ax][e[ax] + 1], e):
                        nb = list(e)
                        nb[ax] += 1
                        if cell_class(cert, tuple(nb)) != 'out':
                            vpend[(tuple(nb), ax)] = [E.t_face(v[c], ax, 1) for c in range(3)]
                vlow[(j, k)] = [E.t_face(v[c], 0, 0) for c in range(3)]
                if i < nx - 1 and face_in_domain(cert, 0, ed[0][i + 1], e):
                    vup[(j, k)] = [E.t_face(v[c], 0, 1) for c in range(3)]
        if pend or vpend:
            st['face_bad'] += len(pend)
            st['vface_bad'] += len(vpend)
            note(st, 'unmatched faces inside slab %d' % i)
        # x faces with the previous slab
        if i == a:
            st['first_low'], st['vfirst_low'] = low, vlow
        else:
            match_x(st, st['up'], low, st['vup'], vlow, i)
        st['up'], st['vup'] = up, vup
        st['next'] = i + 1
        if ckpt:
            tmp = ckpt + '.tmp'
            with open(tmp, 'wb') as f:
                pickle.dump(st, f)
            os.replace(tmp, ckpt)
        print('   slab %d of %d done (%.0f s)' % (i + 1, nx, time.time() - t0), flush=True)
    return st


def match_x(st, up, low, vup, vlow, i):
    """faces between slabs i-1 and i: the upper faces of slab i-1 against the lower faces of slab i"""
    for key, faces in (up or {}).items():
        st['face'] += 1
        other = low.get(key)
        if other is None or not all(E.t_equal(u, w) for u, w in zip(faces, other)):
            st['face_bad'] += 1
            note(st, 'traction jump across the x face %d at %s' % (i, key))
    for key, faces in (vup or {}).items():
        st['vface'] += 1
        other = vlow.get(key)
        if other is None or not all(E.t_equal(u, w) for u, w in zip(faces, other)):
            st['vface_bad'] += 1
            note(st, 'displacement jump across the x face %d at %s' % (i, key))


# ------------------------------------------------------------------------------------------------- driver
class Rep:
    def __init__(self, logfile=None):
        self.fail, self.logfile = [], logfile

    def say(self, s):
        print(s, flush=True)
        if self.logfile:
            with open(self.logfile, 'a', encoding='utf-8') as f:
                f.write(s + '\n')

    def check(self, ok, what):
        self.say(('  [pass] ' if ok else '  [FAIL] ') + what)
        if not ok:
            self.fail.append(what)
        return ok


def coefficient_conditions(cert):
    """the stored-coefficient conditions of FORMAT.md 9 (sufficient conditions; identities are checked separately)"""
    V, Psi = cert.V, cert.Psi
    ok_v = (all(x == 0 for x in V[0][0].flat) and all(x == 0 for x in V[0][-1].flat)
            and all(x == 0 for x in V[1][:, 0, :].flat) and all(x == 0 for x in V[2][:, :, 0].flat))
    if not cert.frictionless:
        ok_v = ok_v and all(x == 0 for x in V[1][-1].flat) and all(x == 0 for x in V[2][-1].flat)
    ok_p = True
    for c, (a, b) in enumerate(COMP):
        A = Psi[c]
        for ax in range(3):
            even = ((a == ax) + (b == ax)) % 2 == 0
            first = np.take(A, 0, axis=ax)
            second = np.take(A, 1, axis=ax)
            ok_p &= bool(np.all(first == second)) if even else all(x == 0 for x in first.flat)
            if ax == 0 and cert.frictionless:
                last, prev = np.take(A, -1, axis=0), np.take(A, -2, axis=0)
                ok_p &= bool(np.all(last == prev)) if even else all(x == 0 for x in last.flat)
            if ax in (1, 2):
                ok_p &= all(x == 0 for x in np.take(A, -1, axis=ax).flat) and \
                    all(x == 0 for x in np.take(A, -2, axis=ax).flat)
    return ok_v, ok_p


def check_cert(meta, arrays, rep, procs=2, cache=None, weight=True, claims_required=True):
    FD.assert_independent()
    rep.check(True, 'no module of benchpack/builder is loaded (only definitions.py is shared)')
    rep.check(meta.get('format') == 'benchpack-certificate-1' and meta.get('problem') == 'P3'
              and meta.get('model') == 'octant', 'format, problem P3, octant model')
    cert = Cert(meta, arrays)
    N = cert.N
    rep.check(cert.variant in ('frictionless', 'clamped', 'nocavity-test', 'nocavity-clamped-test'),
              'variant %r' % cert.variant)
    rep.say('grid N = %d, grade %s: %s cells; displacement degree %d, stress degree %d; integration %s'
            % (N, cert.grade, ' x '.join(str(v) for v in cert.n), cert.pu, cert.ps, meta.get('integration')))
    # the cell edges of FORMAT.md 9: k/N, and with "grade": [L, m] also 1 - (1/N) 2^-(i+1) (1 + j/m)
    expect = {Fr(k, N) for k in range(N + 1)}
    if cert.grade:
        Lg, mg = int(cert.grade[0]), int(cert.grade[1])
        expect |= {1 - Fr(1, N) * Fr(1, 2 ** (i + 1)) * (1 + Fr(j, mg)) for i in range(Lg) for j in range(mg)}
    rep.check(all(set(ed) == expect for ed in cert.edges),
              'cell edges on all three axes are those of FORMAT.md 9 for N = %d and grade %s' % (N, cert.grade))
    mult_bad = []
    for nm, axes, need in (('displacement', cert.dax, cert.pu), ('stress', cert.sax, cert.ps - 1)):
        for d, ax in enumerate(axes):
            rep.check(ax.ok and ax.on_grid, '%s knots %s: open on [0, 1], break points among the cell edges %s'
                      % (nm, 'xyz'[d], ax.msg))
            mult_bad.append(not rep.check(ax.mult <= need, '%s knots %s: interior multiplicity %d <= %d'
                                          % (nm, 'xyz'[d], ax.mult, need)))
    rep.check(cert.V.shape == (3,) + tuple(a.n for a in cert.dax), 'V shape %s' % (cert.V.shape,))
    rep.check(cert.Psi.shape == (6,) + tuple(a.n for a in cert.sax), 'Psi shape %s' % (cert.Psi.shape,))
    rep.check(len(cert.gamma) == len(cert.gexp) and len(cert.delta) == len(cert.gexp), 'gam has 2 n_g entries')
    if cert.cavity:
        rep.check(HALF in cert.edges[1] and HALF in cert.edges[2] and Fr(3, 4) in cert.edges[0],
                  '1/2 is a cell edge in y and z, 3/4 in x')
        rep.check(cert.gbase == HALF and cert.etab == Fr(3, 4), 'g_base = 1/2 and eta_flat_from = 3/4')
        rep.check(cert.mode in ('superset', 'exact'), 'integration mode %r is superset or exact' % cert.mode)
    if len(rep.fail) > sum(mult_bad):          # structural errors: the fields cannot be rebuilt
        return False
    okv, okp = coefficient_conditions(cert)
    rep.check(okv, 'displacement: coefficient conditions of FORMAT.md 9')
    rep.check(okp, 'stress potential: parity and free-face conditions of FORMAT.md 9')
    nslab = cert.n[0]
    bounds = [round(r * nslab / procs) for r in range(procs + 1)]
    jobs = [(bounds[r], bounds[r + 1],
             os.path.join(cache, 'range_%d_%d.pkl' % (bounds[r], bounds[r + 1])) if cache else None)
            for r in range(procs) if bounds[r + 1] > bounds[r]]
    t0 = time.time()
    if len(jobs) > 1:
        with mp_.Pool(len(jobs), initializer=init, initargs=(meta, arrays, weight)) as pool:
            states = pool.map(run_range, jobs)
    else:
        init(meta, arrays, weight)
        states = [run_range(j) for j in jobs]
    rep.say('all slabs evaluated (%.0f s)' % (time.time() - t0))
    tot = new_state(0)
    for s in states:
        for key in ('Q', 'A', 'Qhi', 'Ahi', 'volhi', 'F', 'cells', 'cut', 'div_bad', 'face', 'face_bad', 'sphere_bad',
                    'bc_bad', 'vface', 'vface_bad', 'vbc_bad', 'vol'):
            tot[key] += s[key]
        tot['fails'] += s['fails']
    for s1, s2 in zip(states[:-1], states[1:]):
        match_x(tot, s1['up'], s2['first_low'], s1['vup'], s2['vfirst_low'], jobs[states.index(s2)][0])
    for msg in tot['fails'][:10]:
        rep.say('    ' + msg)
    rep.check(tot['div_bad'] == 0, 'div tau = 0 exactly in all %d cells of the domain' % tot['cells'])
    rep.check(tot['face_bad'] == 0, 'traction continuous across all %d faces meeting the domain' % tot['face'])
    if cert.cavity:
        rep.check(tot['sphere_bad'] == 0, 'zero traction on the sphere: tau.grad(rho) divisible by rho in all %d cut '
                                          'cells' % tot['cut'])
    rep.check(tot['bc_bad'] == 0, 'zero traction on y = 1, z = 1; zero shear on the mirror planes%s'
              % (' and the grip' if cert.frictionless else ''))
    Fform = sum(((g + d) * HALF ** (ek + 1) / (ek + 1) if cert.cavity else (g + d) * (1 - cert.gbase) ** (ek + 1) / (ek + 1)
                 for ek, g, d in zip(cert.gexp, cert.gamma, cert.delta)), Fr(0))
    rep.check(tot['F'] == Fform, 'grip force F_o = %.15g exactly, equal to sum (gamma_k + delta_k) (1/2)^(e_k+1)/(e_k+1)'
              % float(tot['F']))
    rep.check(tot['vface_bad'] == 0 and tot['vbc_bad'] == 0,
              'displacement continuous across all %d faces; mirror and grip values exact' % tot['vface'])
    if cert.cavity:
        if cert.mode == 'exact':
            ex = 1 - iv.pi / 48
            rep.check(RG.ivq(tot['vol']).a <= ex.a and ex.b <= RG.ivq(tot['volhi']).b and tot['volhi'] - tot['vol'] < Fr(1, 10 ** 100),
                      'volume of the integration domain [%.15f, %.15f] encloses |Omega_o| = 1 - pi/48 (exact domain)'
                      % (float(tot['vol']), float(tot['volhi'])))
        else:
            rep.check(tot['vol'] >= 1 - Fr(3141593, 48 * 10 ** 6),
                      'kept stress boxes: volume %.9f >= |Omega_o| = 1 - pi/48 = 0.934550 (superset)' % float(tot['vol']))
    rep.energies = (tot['Q'], tot['Qhi'], tot['A'], tot['Ahi'], tot['F'])
    Q, A, F = tot['Qhi'], tot['Ahi'], tot['F']
    if tot['Qhi'] != tot['Q'] or tot['Ahi'] != tot['A']:
        rep.say('energies (enclosures): Q_o in [%.15g, %.15g] (width %.1e), a_o in [%.15g, %.15g] (width %.1e); the bounds '
                'use the upper ends' % (tot['Q'], tot['Qhi'], tot['Qhi'] - tot['Q'], tot['A'], tot['Ahi'],
                                        tot['Ahi'] - tot['A']))
    Klow = 2 * F * F / Q
    Kup = 8 * A
    rep.say('recomputed (exact rationals): F_o = %.15g, Q_o = %.15g, a_o = %.15g' % (F, Q, A))
    rep.say('recomputed K: [%.12g, %.12g]  width %.3e  relative %.3e'
            % (Klow, Kup, Kup - Klow, (Kup - Klow) / ((Kup + Klow) / 2)))
    rep.bounds = (Klow, Kup)
    cl = meta.get('claims', {}).get('Q1')
    if cl is not None:
        clo, chi = Fr(cl['lower']), Fr(cl['upper'])
        rep.check(clo <= Klow and Kup <= chi, 'claim [%s, %s] confirmed: %s <= %.12g and %.12g <= %s'
                  % (cl['lower'], cl['upper'], cl['lower'], Klow, Kup, cl['upper']))
    elif claims_required:
        rep.check(False, 'claim present')
    if meta.get('computed'):
        rep.say("builder's own values (diagnostic): %s" % str(meta['computed'])[:1500])
    return not rep.fail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('stem')
    ap.add_argument('--procs', type=int, default=2)
    ap.add_argument('--log', default=None)
    a = ap.parse_args()
    stem = a.stem[:-5] if a.stem.endswith(('.json')) else (a.stem[:-4] if a.stem.endswith('.npz') else a.stem)
    meta, arrays = load(stem)
    h = hashlib.sha256()
    for p in (stem + '.json', stem + '.npz', os.path.abspath(__file__), os.path.join(HERE, 'exact3.py'),
              os.path.join(HERE, 'exact.py')):
        with open(p, 'rb') as f:
            h.update(f.read())
    cache = os.path.join(HERE, 'cache', '%s_%s' % (os.path.basename(stem), h.hexdigest()[:16]))
    os.makedirs(cache, exist_ok=True)
    rep = Rep(a.log)
    t0 = time.time()
    rep.say('p3check: %s   (checkpoints in %s)' % (os.path.basename(stem), os.path.relpath(cache, HERE)))
    ok = check_cert(meta, arrays, rep, procs=a.procs, cache=cache)
    rep.say('RESULT: %s   (%.0f s)' % ('ALL CHECKS PASSED' if ok else 'FAILED: %d check(s)' % len(rep.fail),
                                       time.time() - t0))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    main()
