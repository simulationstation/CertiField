#!/usr/bin/env python3
"""exact_check.py -- exact rational check that a stored stress potential gives an admissible stress field.

It does not use the evaluation code of curved3d.py: the stress tau = tau_p + inc(rho^2 Psi) is built on each grid cell as a
polynomial with rational coefficients (the stored floats are exact binary fractions), and the conditions the lower bound
rests on are tested as exact identities:

  (a) div tau = 0 in every cell;
  (b) the normal traction is continuous across every face between two cells;
  (c) on every cell cut by the curved surface, tau.grad(rho) is divisible by rho  (zero traction on the surface);
  (d) zero traction on the free faces y = 1, z = 1; zero shear on the mirror planes (and on a frictionless grip);
  (e) the potential part transmits no force through the grip, the particular part transmits exactly area/4 (octant);
  (f) on a few boxes, the exact integral of tau:C^-1:tau (and of the strain energy of v) lies below the interval bound
      of curved3d.py and agrees with its float value;
  (g) the boxes used for the integrals cover the body: every finest-level box that is not covered lies in the closed
      ball (full enumeration where it fits in memory, otherwise two million random points next to the surface);
  (h) if a displacement is given: it is continuous across all cell faces and satisfies the grip and mirror conditions.

Shared with curved3d.py: only the definition of the body (spline shapes, the 1D functions of the geometry, the box list).

usage:  python3 exact_check.py <grip> <file with psi (and v)> [number of sampled cells; 0 or absent: all cells] [--v file]
        python3 exact_check.py <grip> <file with psi and v> total <depth for psi> <depth for v>
            the second form recomputes BOTH bounds of the stored fields in exact rational arithmetic (slow: coarse grids)
"""
import itertools
import random
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import curved3d as c3
from curved3d import Problem, SHAPES, COMPS, CIDX, psi_ties, v_ties, ties_hold


class Poly:
    """polynomial in t1, t2, t3 with Fraction coefficients: {(e1, e2, e3): coefficient}"""
    __slots__ = ('c',)

    def __init__(self, c=None):
        self.c = {} if c is None else c

    @staticmethod
    def const(a):
        a = Fr(a)
        return Poly({(0, 0, 0): a}) if a else Poly()

    @staticmethod
    def uni(coefs, ax):
        d = {}
        for k, a in enumerate(coefs):
            if a:
                e = [0, 0, 0]
                e[ax] = k
                d[tuple(e)] = Fr(a)
        return Poly(d)

    def __add__(self, o):
        d = dict(self.c)
        for k, a in o.c.items():
            v = d.get(k)
            v = a if v is None else v + a
            if v:
                d[k] = v
            else:
                d.pop(k, None)
        return Poly(d)

    def __sub__(self, o):
        d = dict(self.c)
        for k, a in o.c.items():
            v = d.get(k)
            v = -a if v is None else v - a
            if v:
                d[k] = v
            else:
                d.pop(k, None)
        return Poly(d)

    def scale(self, a):
        a = Fr(a)
        return Poly({k: v * a for k, v in self.c.items()}) if a else Poly()

    def __mul__(self, o):
        d = {}
        for (a1, a2, a3), x in self.c.items():
            for (b1, b2, b3), y in o.c.items():
                k = (a1 + b1, a2 + b2, a3 + b3)
                v = d.get(k)
                d[k] = x * y if v is None else v + x * y
        return Poly({k: v for k, v in d.items() if v})

    def diff(self, ax):
        d = {}
        for k, v in self.c.items():
            if k[ax]:
                kk = list(k)
                kk[ax] -= 1
                d[tuple(kk)] = v * k[ax]
        return Poly(d)

    def subs(self, ax, val):
        d = {}
        for k, v in self.c.items():
            kk = list(k)
            p = kk[ax]
            kk[ax] = 0
            kk = tuple(kk)
            w = v * val ** p
            x = d.get(kk)
            d[kk] = w if x is None else x + w
        return Poly({k: v for k, v in d.items() if v})

    def integral(self, lo, hi):
        tot = Fr(0)
        for k, v in self.c.items():
            t = v
            for ax in range(3):
                t *= (hi[ax] ** (k[ax] + 1) - lo[ax] ** (k[ax] + 1)) / (k[ax] + 1)
            tot += t
        return tot

    def deg(self, ax):
        return max((k[ax] for k in self.c), default=-1)

    def is_zero(self):
        return not self.c

    def __eq__(self, o):
        return self.c == o.c


def levi(i, j, k):
    return (i - j) * (j - k) * (k - i) // 2


def spline_poly(arr, e, basis):
    """exact polynomial of sum c_i N_i1(t1) N_i2(t2) N_i3(t3) on cell e"""
    acc = Poly()
    for i1, i2, i3 in itertools.product(range(3), repeat=3):
        c = float(arr[e[0] + i1, e[1] + i2, e[2] + i3])
        if c != 0.0:
            acc = acc + basis[(i1, i2, i3)].scale(Fr(c))
    return acc


def make_basis():
    one = {}
    for ax in range(3):
        for m in range(3):
            one[(ax, m)] = Poly.uni(SHAPES[m], ax)
    return {(i1, i2, i3): one[(0, i1)] * one[(1, i2)] * one[(2, i3)] for i1, i2, i3 in itertools.product(range(3), repeat=3)}


def rho_poly(P, e):
    rho = Poly.const(-P.a2)
    for k in range(3):
        if P.eta[k] is not None:
            rho = rho + Poly.uni(P.funcs[P.eta[k]](e[k]), k)
    return rho


def sep_poly(P, e, names):
    r = Poly.const(1)
    for ax in range(3):
        r = r * Poly.uni(P.funcs[names[ax]](e[ax]), ax)
    return r


def stress_polys(P, psi, e, basis):
    """(potential part, total) stress on cell e as dicts (i, j) -> Poly, physical units (d/dx = N d/dt); and rho"""
    N = P.N
    rho = rho_poly(P, e)
    R = rho * rho
    Phi = [R * spline_poly(psi[p], e, basis) for p in range(6)]
    for (p, coef, nx, ny, nz) in P.pterms:
        Phi[p] = Phi[p] + sep_poly(P, e, (nx, ny, nz)).scale(coef)
    D2 = {}
    for p in range(6):
        for k in range(3):
            dk = Phi[p].diff(k)
            for m in range(k, 3):
                D2[(p, k, m)] = dk.diff(m).scale(N * N)
    pot = {}
    for i in range(3):
        for j in range(i, 3):
            acc = Poly()
            for k, l, m, n in itertools.product(range(3), repeat=4):
                s = levi(i, k, l) * levi(j, m, n)
                if s:
                    term = D2[(CIDX[(l, n)], min(k, m), max(k, m))]
                    acc = acc + term if s > 0 else acc - term
            pot[(i, j)] = acc
            pot[(j, i)] = acc
    tot = dict(pot)
    for (s, coef, nx, ny, nz) in P.sterms:
        i, j = COMPS[s]
        tot[(i, j)] = tot[(i, j)] + sep_poly(P, e, (nx, ny, nz)).scale(coef)
        tot[(j, i)] = tot[(i, j)]
    return pot, tot, rho


def remainder(g, rho, ax):
    """remainder of g modulo rho, as polynomials in t_ax (rho is quadratic in t_ax with a constant leading coefficient)"""
    key = [0, 0, 0]
    key[ax] = 2
    lead = rho.c[tuple(key)]
    r = g
    while not r.is_zero() and r.deg(ax) >= 2:
        d = r.deg(ax)
        q = {}
        for k, v in r.c.items():
            if k[ax] == d:
                kk = list(k)
                kk[ax] -= 2
                q[tuple(kk)] = v / lead
        r = r - Poly(q) * rho
    return r


def cell_class(P, e):
    N, a2 = P.N, P.a2
    rmin = sum(P.kappa[k] * e[k] ** 2 for k in range(3))
    rmax = sum(P.kappa[k] * (e[k] + 1) ** 2 for k in range(3))
    thr = a2 * N * N
    if P.solid_inside:
        return 0 if rmax <= thr else (2 if rmin >= thr else 1)
    return 0 if rmin >= thr else (2 if rmax <= thr else 1)


def q_dual(P, tot):
    nu = P.nu
    diag = [tot[(i, i)] for i in range(3)]
    off = [tot[(0, 1)], tot[(1, 2)], tot[(0, 2)]]
    acc = Poly()
    for t in diag:
        acc = acc + t * t
    for t in off:
        acc = acc + (t * t).scale(2)
    tr = diag[0] + diag[1] + diag[2]
    return acc.scale(1 + nu) - (tr * tr).scale(nu)


def q_primal(P, v, e, basis):
    N = P.N
    vp = [spline_poly(v[a], e, basis) for a in range(3)]
    dv = [[vp[a].diff(k).scale(N) for k in range(3)] for a in range(3)]
    exx = dv[0][0] + Poly.const(Fr(1, 2))
    eyy, ezz = dv[1][1], dv[2][2]
    gxy, gyz, gxz = dv[0][1] + dv[1][0], dv[1][2] + dv[2][1], dv[0][2] + dv[2][0]
    tr = exx + eyy + ezz
    acc = (tr * tr).scale(P.lam)
    for t in (exx, eyy, ezz):
        acc = acc + (t * t).scale(2 * P.mu)
    for t in (gxy, gyz, gxz):
        acc = acc + (t * t).scale(P.mu)
    return acc


def coverage_check(P, L, seed=3, nsamples=2_000_000):
    """(g) every point of the octant body lies in a kept box of depth L"""
    boxes = P.boxes(L)
    N, kap = P.N, np.array(P.kappa, dtype=np.int64)
    n = N * 2 ** L
    num, den = P.a2.numerator, P.a2.denominator
    if n ** 3 <= 6e7:
        cover = np.zeros((n, n, n), dtype=bool)
        for b in boxes:
            s_ = 2 ** (L - int(b[3]))
            x0 = [(int(b[k]) * 2 ** int(b[3]) + int(b[4 + k])) * s_ for k in range(3)]
            cover[x0[0]:x0[0] + s_, x0[1]:x0[1] + s_, x0[2]:x0[2] + s_] = True
        idx = np.argwhere(~cover).astype(np.int64)
        if P.solid_inside:
            ok = ((idx ** 2) * kap).sum(axis=1) * den >= num * n * n
        else:
            ok = (((idx + 1) ** 2) * kap).sum(axis=1) * den <= num * n * n
        assert bool(np.all(ok)), "an uncovered box reaches into the body"
        return f"depth {L}: all {n}^3 finest boxes enumerated, {len(idx)} uncovered, all inside the void"
    rng = np.random.default_rng(seed)
    f = 8                                                   # sample points are centres of a lattice 8 times finer
    u = np.abs(rng.standard_normal((nsamples, 3)))
    u /= np.linalg.norm(u, axis=1)[:, None]
    r = float(P.a) + rng.uniform(-0.5, 3.0, size=(nsamples, 1)) / N
    c = np.floor(u * r * n * f).astype(np.int64)
    c = c[np.all((c >= 0) & (c < n * f), axis=1)]
    r2 = (((2 * c + 1) ** 2) * kap).sum(axis=1)             # (2 x n f)^2 at the point x = (c + 1/2)/(n f)
    inside_body = (r2 * den < num * (2 * n * f) ** 2) if P.solid_inside else (r2 * den > num * (2 * n * f) ** 2)
    c = c[inside_body]
    covered = np.zeros(len(c), dtype=bool)
    for lev in range(L + 1):
        m = N * 2 ** lev
        bl = boxes[boxes[:, 3] == lev]
        if not len(bl):
            continue
        key = ((bl[:, 0] * 2 ** lev + bl[:, 4]) * m + (bl[:, 1] * 2 ** lev + bl[:, 5])) * m + (bl[:, 2] * 2 ** lev + bl[:, 6])
        q = c // (f * 2 ** (L - lev))
        covered |= np.isin((q[:, 0] * m + q[:, 1]) * m + q[:, 2], key)
    assert bool(np.all(covered)), "a point of the body is not covered by the boxes"
    return f"depth {L}: {len(c)} random points of the body next to the surface, all covered"


def check_v(P, v, basis, grip):
    """(h) exact admissibility of the displacement  (x/2) e_x + spline"""
    N = P.N
    Pm, own = v_ties(N, grip, owners=True)
    assert ties_hold(Pm, own, v), "displacement: ties violated"
    faces, nface = {}, 0
    one, zero = Fr(1), Fr(0)
    for e in itertools.product(range(N), repeat=3):
        vp = [spline_poly(v[a], e, basis) for a in range(3)]
        for k in range(3):
            if (e, k) in faces:
                assert faces.pop((e, k)) == [vp[a].subs(k, zero) for a in range(3)], ("displacement jump", e, k)
                nface += 1
            if e[k] < N - 1:
                nb = list(e)
                nb[k] += 1
                faces[(tuple(nb), k)] = [vp[a].subs(k, one) for a in range(3)]
            if e[k] == 0:
                assert vp[k].subs(k, zero).is_zero(), ("normal displacement on a mirror plane", e, k)
        if e[0] == N - 1:
            assert vp[0].subs(0, one).is_zero(), ("grip: u_x", e)
            if grip == 'clamped':
                assert vp[1].subs(0, one).is_zero() and vp[2].subs(0, one).is_zero(), ("grip: u_y, u_z", e)
    assert not faces
    return nface


def mono_tables(L, maxdeg):
    """exact integrals of t^k over every dyadic sub-interval of [0,1] up to depth L"""
    tab = {}
    for lev in range(L + 1):
        s = Fr(1, 2 ** lev)
        for j in range(2 ** lev):
            lo, hi = j * s, (j + 1) * s
            tab[(lev, j)] = [(hi ** (k + 1) - lo ** (k + 1)) / (k + 1) for k in range(maxdeg + 1)]
    return tab


def exact_total(P, boxes, qfun, maxdeg, label, t0):
    """exact sum over the boxes of the integral of the cell polynomial qfun(e) (octant, physical units)"""
    tab = mono_tables(int(boxes[:, 3].max()), maxdeg)
    order = np.lexsort((boxes[:, 2], boxes[:, 1], boxes[:, 0]))
    total, cur, q = Fr(0), None, None
    for n, b in enumerate(boxes[order]):
        e = (int(b[0]), int(b[1]), int(b[2]))
        if e != cur:
            cur, q = e, list(qfun(e).c.items())
        lev = int(b[3])
        I0, I1, I2 = tab[(lev, int(b[4]))], tab[(lev, int(b[5]))], tab[(lev, int(b[6]))]
        acc = Fr(0)
        for (k0, k1, k2), cf in q:
            acc += cf * I0[k0] * I1[k1] * I2[k2]
        total += acc
        if (n + 1) % 20000 == 0:
            print(f"   {label}: {n + 1}/{len(boxes)} boxes   {time.time() - t0:.0f}s", flush=True)
    return total / P.N ** 3


def main_total():
    grip, path, Lp, Lv = sys.argv[1], sys.argv[2], int(sys.argv[4]), int(sys.argv[5])
    data = np.load(path)
    psi, v = data['psi'], data['v']
    N = psi.shape[1] - 2
    assert v.shape[1] == psi.shape[1]
    P = Problem(N, grip=grip, particular='step', column=Fr(1, 2), kflat=(3 * N // 4 if grip == 'frictionless' else None),
                L_max=max(Lp, Lv))
    basis = make_basis()
    t0 = time.time()
    print(f"grid {N}, {grip} grips: both bounds in exact rational arithmetic (depths {Lp} and {Lv})", flush=True)
    bp, bv = P.boxes(Lp), P.boxes(Lv)
    Q = 8 * exact_total(P, bp, lambda e: q_dual(P, stress_polys(P, psi, e, basis)[1]), 12, "stress", t0)
    K_low = P.area[0] ** 2 / Q
    W = 8 * exact_total(P, bv, lambda e: q_primal(P, v, e, basis), 4, "strain", t0)
    Q_up = float(np.nextafter(8 * c3.dual_compliance(P, bp, psi, rigorous=True), np.inf))
    alo = c3.frac_enclosure(P.area[0])[0]
    c_low = float(c3.dn(c3.dn(alo * alo) / Q_up))
    c_up = float(np.nextafter(8 * c3.primal_energy(P, bv, v, rigorous=True), np.inf))
    print(f"   exact lower bound  {float(K_low):.15f}     interval result {c_low:.15f}    exact - interval = {float(K_low - Fr(c_low)):.2e}")
    print(f"   exact upper bound  {float(W):.15f}     interval result {c_up:.15f}    interval - exact = {float(Fr(c_up) - W):.2e}")
    assert Fr(c_low) <= K_low and W <= Fr(c_up), "an interval result is on the wrong side of the exact value"
    assert K_low - Fr(c_low) < Fr(1, 10 ** 9) and Fr(c_up) - W < Fr(1, 10 ** 9)
    assert K_low <= W
    print(f"EXACT TOTALS PASSED   {time.time() - t0:.0f}s")


def main():
    if len(sys.argv) > 3 and sys.argv[3] == 'total':
        return main_total()
    argv = list(sys.argv)
    vfile = None
    if '--v' in argv:
        i = argv.index('--v')
        vfile = argv[i + 1]
        del argv[i:i + 2]
    grip, path = argv[1], argv[2]
    nsample = int(argv[3]) if len(argv) > 3 and int(argv[3]) > 0 else None
    data = np.load(path)
    psi = data['psi']
    v = np.load(vfile)['v'] if vfile else (data['v'] if 'v' in data else None)
    N = psi.shape[1] - 2
    P = Problem(N, grip=grip, particular='step', column=Fr(1, 2), kflat=(3 * N // 4 if grip == 'frictionless' else None), L_max=3)
    Pm, own = psi_ties(N, grip, owners=True)
    assert ties_hold(Pm, own, psi), "ties violated"
    basis = make_basis()
    t0 = time.time()
    cells = list(itertools.product(range(N), repeat=3))
    if nsample is not None and nsample < len(cells):
        rnd = random.Random(1)
        cut = [e for e in cells if cell_class(P, e) == 1]
        grip_cells = [e for e in cells if e[0] == N - 1]
        cells = rnd.sample(cut, min(len(cut), nsample // 3)) + rnd.sample(grip_cells, min(len(grip_cells), nsample // 6)) \
            + rnd.sample(cells, nsample - nsample // 3 - nsample // 6)
        cells = sorted(set(cells))
        sampled = True
    else:
        sampled = False
    print(f"grid {N}, {grip} grips: exact check of {len(cells)} cells" + (" (sample)" if sampled else " (all)"), flush=True)
    faces = {}
    cellset = set(cells)
    count = dict(div=0, sphere=0, free=0, mirror=0, grip=0, face=0)
    res_pot, res_tot = Fr(0), Fr(0)
    one = Fr(1)
    for idx, e in enumerate(cells):
        pot, tot, rho = stress_polys(P, psi, e, basis)
        # (a) equilibrium
        for i in range(3):
            d = Poly()
            for j in range(3):
                d = d + tot[(i, j)].diff(j)
            assert d.is_zero(), ("div tau != 0", e, i)
        count['div'] += 1
        # (c) curved surface
        if cell_class(P, e) == 1:
            ax = next(k for k in range(3) if P.eta[k] is not None)
            for i in range(3):
                g = Poly()
                for j in range(3):
                    g = g + tot[(i, j)] * rho.diff(j)
                assert remainder(g, rho, ax).is_zero(), ("traction on the curved surface", e, i)
            count['sphere'] += 1
        # (d) flat boundary faces
        for k in (1, 2):
            if e[k] == N - 1:
                for i in range(3):
                    assert tot[(i, k)].subs(k, one).is_zero(), ("traction on a free face", e, k, i)
                count['free'] += 1
        for k in range(3):
            if e[k] == 0:
                for i in range(3):
                    if i != k:
                        assert tot[(i, k)].subs(k, Fr(0)).is_zero(), ("shear on a mirror plane", e, k, i)
                count['mirror'] += 1
        if e[0] == N - 1:
            if grip == 'frictionless':
                for i in (1, 2):
                    assert tot[(i, 0)].subs(0, one).is_zero(), ("shear on the grip", e, i)
            count['grip'] += 1
            lo, hi = (Fr(0),) * 3, (Fr(1),) * 3
            res_pot += pot[(0, 0)].subs(0, one).integral(lo, hi) / (N * N)
            res_tot += tot[(0, 0)].subs(0, one).integral(lo, hi) / (N * N)
        # (b) compare with the faces left by the lower neighbours, leave the upper faces for the next cells
        for k in range(3):
            key = (e, k)
            if key in faces:
                assert faces.pop(key) == [tot[(i, k)].subs(k, Fr(0)) for i in range(3)], ("traction jump across a cell face", e, k)
                count['face'] += 1
            if e[k] < N - 1:
                nb = list(e)
                nb[k] += 1
                nb = tuple(nb)
                if nb in cellset:
                    faces[(nb, k)] = [tot[(i, k)].subs(k, one) for i in range(3)]
        if (idx + 1) % 500 == 0:
            print(f"   {idx + 1} cells   {time.time() - t0:.0f}s", flush=True)
    assert not faces, "unmatched faces"
    print(f"   (a) div tau = 0 exactly in {count['div']} cells")
    print(f"   (b) normal traction continuous across {count['face']} interior faces")
    print(f"   (c) zero traction on the curved surface, as a polynomial identity, in {count['sphere']} cut cells")
    print(f"   (d) zero traction on {count['free']} free-face cells; zero shear on {count['mirror']} mirror-plane cells"
          + (f" and on the grip" if grip == 'frictionless' else ""))
    if not sampled:
        assert res_pot == 0 and res_tot == P.area[0] / 4, (res_pot, res_tot)
        print(f"   (e) force through the grip (octant): potential part {res_pot}, total {res_tot} (expected 0 and {P.area[0] / 4})")
    # (f) exact integrals on a few boxes against curved3d
    bx_all = P.boxes(3)
    rnd = random.Random(2)
    pick = rnd.sample([b for b in bx_all if b[7] == 1], 3) + rnd.sample([b for b in bx_all if b[3] == 0], 3)
    worst = 0.0
    for b in pick:
        e = tuple(int(x) for x in b[:3])
        s = Fr(1, 2 ** int(b[3]))
        lo = tuple(int(b[4 + k]) * s for k in range(3))
        hi = tuple(x + s for x in lo)
        _, tot, _ = stress_polys(P, psi, e, basis)
        exact = q_dual(P, tot).integral(lo, hi) / N ** 3
        b1 = np.array([b])
        up = c3.dual_compliance(P, b1, psi, rigorous=True)
        fl = c3.dual_compliance(P, b1, psi, rigorous=False)
        assert exact >= 0 and Fr(up) >= exact, ("interval bound below the exact integral", b)
        rel_up, rel_fl = float((Fr(up) - exact) / exact), float(abs(Fr(fl) - exact) / exact)
        worst = max(worst, rel_up, rel_fl)
        assert rel_up < 1e-8 and rel_fl < 1e-8
        if v is not None and v.shape[1] == psi.shape[1]:
            exv = q_primal(P, v, e, basis).integral(lo, hi) / N ** 3
            upv = c3.primal_energy(P, b1, v, rigorous=True)
            assert Fr(upv) >= exv and float((Fr(upv) - exv) / exv) < 1e-9, ("strain energy", b)
    print(f"   (f) exact integrals on {len(pick)} boxes: the interval bounds are above them, by at most {worst:.1e} (relative)")
    for L in (3, 4, 5, 6):
        print(f"   (g) {coverage_check(P, L)}   {time.time() - t0:.0f}s", flush=True)
    if v is not None:
        nf = check_v(P, v, basis, grip)
        print(f"   (h) displacement: continuous across {nf} faces, grip and mirror conditions exact")
    print(f"EXACT CHECK PASSED   {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
