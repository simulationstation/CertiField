"""Tests of the independent checker.

  python tests.py            (about two minutes, one process)

1. B-splines: exact reproduction of linear functions, smoothness C^{p-m} across a knot of multiplicity m.
2. Curved integrals: the quarter domain and the half patch P_h against closed forms (pi, sin(1/10), cos(1/10)), on
   grids with grid points exactly on the circles; one cut cell against mpmath quadrature (diagnostic agreement).
3. Known answers: the plate without the hole, K = 1 and Q2 = 1/2 (the frictionless grips give sigma_xx = 1/2).
4. Green identity with random admissible fields: int tau:eps(v) = (work factor) * F(tau), exactly for P2 and to
   1e-100 for P1 (this needs equilibrium, every boundary condition and the cut-cell integrals to be right).
5. Falsification: random admissible fields never give a lower bound above an upper bound (K, and m(k) for Q2);
   near-optimal fields (ritz.py, a small float Ritz solver) give close bounds that do not cross, also after random
   admissible perturbations of relative size 1e-8 to 1e-2.
6. Negative tests: broken certificates are rejected (shear on a mirror line, a C^0 stress function, a wrong grip
   value, a discontinuous displacement, a stress field not traction-free on the hole, claims tighter than the bounds,
   a free edge of P2 with traction).
7. P3: the solid cube (K = 2 exactly); the Green identity on the solid cube (exact); octree box integrals against
   direct integration; random fields with the cavity; negative tests (mirror shear, no weight, C^0 potential, claim).
8. P3 exact mode: the closed forms of ball.py (exact identities); moments of the octant of the ball against the
   Gamma-function values; the lift alone, a_o = (lambda + 2 mu)(1 - pi/48)/4; the Green identity over the octant
   minus the ball for random admissible fields (frictionless and clamped).
"""
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from mpmath import iv, mp

import check
import exact as X
import fields as FD
import region as RG
import ritz
import synth

FAILS = []


def ok(cond, what):
    print(('  ok    ' if cond else '  FAIL  ') + what, flush=True)
    if not cond:
        FAILS.append(what)


def encl(lo, hi, I):
    """[lo, hi] (Fractions) encloses the interval I"""
    return RG.ivq(lo).a <= I.a and I.b <= RG.ivq(hi).b


# ------------------------------------------------------------------------------------------------- 1
def test_bspline():
    print('1. B-splines')
    for p, kn in ((2, [0, 0, 0, Fr(1, 3), Fr(1, 2), Fr(1, 2), 1, 1, 1]),
                  (5, [0] * 6 + [Fr(1, 7), Fr(2, 7), Fr(2, 7), Fr(1, 2), Fr(1, 2), Fr(1, 2), Fr(5, 7)] + [1] * 6)):
        kn = [Fr(v) for v in kn]
        g = synth.greville(kn, p)
        spans = X.bspline_pieces(kn, p)
        good = True
        for (k, t0, t1, M) in spans:
            poly = X.ptrim(list(np.array(g[k - p:k + 1], dtype=object).dot(M)))
            good &= poly == X.ptrim([t0, t1 - t0])          # sum xi_j B_j(x) = x = t0 + h s
        ok(good, 'degree %d: sum of Greville abscissae times B-splines equals x exactly on every span' % p)
        # smoothness across each interior knot: derivatives 0..p-m agree
        good = True
        for (k, t0, t1, M), (k2, u0, u1, M2) in zip(spans[:-1], spans[1:]):
            m = kn.count(t1)
            h1, h2 = t1 - t0, u1 - u0
            for i in range(k - p, k + 1):
                if not (k2 - p <= i <= k2):
                    continue
                a = list(M[i - (k - p), :])
                b = list(M2[i - (k2 - p), :])
                for d in range(p - m + 1):
                    va = X.peval(a, Fr(1)) / h1 ** d
                    vb = X.peval(b, Fr(0)) / h2 ** d
                    good &= va == vb
                    a, b = X.pder(a), X.pder(b)
        ok(good, 'degree %d: B-splines are C^(p-m) across knots of multiplicity m (exact)' % p)


# ------------------------------------------------------------------------------------------------- 2
def grid_cells(n):
    h = Fr(1, n)
    return [(i * h, h, j * h, h) for i in range(n) for j in range(n)]


def test_region():
    print('2. curved integrals')
    one = X.zeros((1, 1)) + 1
    for n in (10, 16):
        lo = hi = Fr(0)
        lo2 = hi2 = Fr(0)
        for (x0, h, y0, k) in grid_cells(n):
            x1, y1 = x0 + h, y0 + k
            if x1 * x1 + y1 * y1 <= Fr(1, 4):
                continue
            P2 = FD.local_poly_from_dict({(2, 0): Fr(1)}, x0, h, y0, k)
            if x0 * x0 + y0 * y0 >= Fr(1, 4):
                a = h * k
                b = h * k * X.integral(P2)
                lo, hi, lo2, hi2 = lo + a, hi + a, lo2 + b, hi2 + b
            else:
                I = RG.integrate_region(one, x0, h, y0, k, [('circle', Fr(1, 4))], [])
                J = RG.integrate_region(P2, x0, h, y0, k, [('circle', Fr(1, 4))], [])
                lo, hi = lo + RG.lo_q(I), hi + RG.hi_q(I)
                lo2, hi2 = lo2 + RG.lo_q(J), hi2 + RG.hi_q(J)
        ok(encl(lo, hi, 1 - iv.pi / 16) and hi - lo < Fr(1, 10 ** 300),
           'grid 1/%d: area of the quarter minus the disc encloses 1 - pi/16 (width %.1e)' % (n, float(hi - lo)))
        ok(encl(lo2, hi2, RG.ivq(Fr(1, 3)) - iv.pi / 256) and hi2 - lo2 < Fr(1, 10 ** 300),
           'grid 1/%d: int x^2 over the quarter minus the disc encloses 1/3 - pi/256' % n)
    # half patch: closed forms in polar coordinates
    r1, r2 = Fr(1, 2), Fr(11, 20)
    s, c = iv.sin(RG.ivq(Fr(1, 10))), iv.cos(RG.ivq(Fr(1, 10)))
    exact = {'1': RG.ivq(Fr(21, 8000)), 'x': RG.ivq((r2 ** 3 - r1 ** 3) / 3) * (1 - c),
             'y': RG.ivq((r2 ** 3 - r1 ** 3) / 3) * s, 'x2+y2': RG.ivq((r2 ** 4 - r1 ** 4) / 4 / 10),
             'x^3 y^5': None}
    polys = {'1': {(0, 0): Fr(1)}, 'x': {(1, 0): Fr(1)}, 'y': {(0, 1): Fr(1)},
             'x2+y2': {(2, 0): Fr(1), (0, 2): Fr(1)}, 'x^3 y^5': {(3, 5): Fr(1)}}
    # x^3 y^5 in polar: int r^9 dr * int cos^3 sin^5 dt, by mpmath quadrature (diagnostic)
    with mp.workdps(60):
        ref = ((mp.mpf(11) / 20) ** 10 - (mp.mpf(1) / 2) ** 10) / 10 * mp.quad(
            lambda t: mp.cos(t) ** 3 * mp.sin(t) ** 5, [mp.pi / 2 - mp.mpf(1) / 10, mp.pi / 2])
    for n in (40, 100, 200):          # 200: grid lines through (33/100, 44/100) on r = 11/20 are not met; fine
        tot = {key: [Fr(0), Fr(0)] for key in polys}
        for (x0, h, y0, k) in grid_cells(n):
            if x0 >= Fr(11, 200) or y0 + k <= Fr(497, 1000) or y0 >= Fr(11, 20):
                continue
            for key, pd in polys.items():
                P = FD.local_poly_from_dict(pd, x0, h, y0, k)
                I = RG.integrate_region(P, x0, h, y0, k, list(FD.PATCH_LOWERS), list(FD.PATCH_UPPERS))
                tot[key][0] += RG.lo_q(I)
                tot[key][1] += RG.hi_q(I)
        for key in polys:
            lo, hi = tot[key]
            if exact[key] is not None:
                ok(encl(lo, hi, exact[key]) and hi - lo < Fr(1, 10 ** 200),
                   'grid 1/%d: int_{P_h} %s encloses the closed form (width %.1e)' % (n, key, float(hi - lo)))
            else:
                with mp.workdps(60):
                    err = abs(mp.mpf(lo.numerator) / lo.denominator - ref)
                ok(err < mp.mpf(10) ** -40 and hi - lo < Fr(1, 10 ** 200),
                   'grid 1/%d: int_{P_h} x^3 y^5 agrees with mpmath quadrature to 1e-40 (diagnostic)' % n)
    # a cut cell whose corner (3/10, 2/5) lies on the circle, with a random polynomial, against quadrature
    rng = np.random.default_rng(5)
    P = X.farr(rng.standard_normal((5, 6)))
    x0, h, y0, k = Fr(3, 10), Fr(1, 10), Fr(3, 10), Fr(1, 10)
    I = RG.integrate_region(P, x0, h, y0, k, [('circle', Fr(1, 4))], [])

    def f(x, y):
        s_, t_ = (x - mp.mpf(3) / 10) * 10, (y - mp.mpf(3) / 10) * 10
        return sum(mp.mpf(P[a, b].numerator) / P[a, b].denominator * s_ ** a * t_ ** b
                   for a in range(5) for b in range(6))
    # the cell [3/10, 2/5]^2 minus the disc: for x in [3/10, 2/5], y from sqrt(1/4 - x^2) (>= 3/10) to 2/5
    with mp.workdps(50):
        ref = mp.quad(lambda x: mp.quad(lambda y: f(x, y), [mp.sqrt(mp.mpf(1) / 4 - x * x), mp.mpf(4) / 10]),
                      [mp.mpf(3) / 10, mp.mpf(4) / 10])
        mid = mp.mpf(RG.lo_q(I).numerator) / RG.lo_q(I).denominator
        err = abs(mid - ref)
    ok(err < mp.mpf(10) ** -30 and RG.hi_q(I) - RG.lo_q(I) < Fr(1, 10 ** 200),
       'cut cell with a corner on the circle: closed form agrees with mpmath quadrature (%.1e, diagnostic)' % float(err))


# ------------------------------------------------------------------------------------------------- 3
def run(meta, arrays):
    rep = check.Report()
    import io
    import contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        good = check.check_cert(meta, arrays, rep, procs=1, claims_required=False)
    return good, rep


def test_known_answer():
    print('3. known answers: plate without the hole')
    meta, arrays = synth.make_nohole()
    good, rep = run(meta, arrays)
    ok(good, 'all checks pass on the exact fields of the plate without the hole')
    ok(rep.bounds['Q1'] == (1, 1), 'K: recomputed bounds are [1, 1] exactly (known answer K = 1)')
    n = rep.numbers
    ok(n[('l', 'v')][0] <= Fr(1, 2) <= n[('l', 'v')][1] and n[('l', 'tau')][0] <= Fr(1, 2) <= n[('l', 'tau')][1]
       and n[('l', 'v')][1] - n[('l', 'v')][0] < Fr(1, 10 ** 200),
       'Q2: the patch average of sigma_xx of the exact fields encloses 1/2 (width < 1e-200)')
    lo, hi = rep.bounds['Q2']
    ok(lo <= Fr(1, 2) <= hi, 'Q2: the duality bounds [%.6f, %.6f] contain 1/2' % (float(lo), float(hi)))
    # for these fields (no adjoint part) m_up(+-k) = 1/8 -+ k/2 and m_low(+-k) = 1/8 -+ k/2 - k^2 c/2, c = C11/|P_h|,
    # so the bracket is [1/2 - k c/4, 1/2 + k c/4]
    c = FD.D.C[0][0] / FD.D.P1_QUARTER['half_patch_area']
    w = Fr(1, 1000) * c / 4
    ok(abs(lo - (Fr(1, 2) - w)) < Fr(1, 10 ** 100) and abs(hi - (Fr(1, 2) + w)) < Fr(1, 10 ** 100),
       'Q2: bracket equals [1/2 - kappa c/4, 1/2 + kappa c/4], c = C11/|P_h|, as derived by hand')


# ------------------------------------------------------------------------------------------------- 4
def cross_integral(meta, arrays, vname, tname):
    """int over the domain of tau : eps(v) (same cell grid for both spaces); interval (lo, hi)"""
    prob = FD.Problem(meta)
    dsp, ssp = FD.Space(meta['displacement_space']), FD.Space(meta['stress_space'])
    assert [s[1] for s in dsp.sx] == [s[1] for s in ssp.sx]
    A, B = arrays[vname + '_A'], arrays[vname + '_B']
    C = arrays[tname + '_C']
    par = list(arrays[tname + '_par'])
    lo = hi = Fr(0)
    for (i, j) in dsp.cells():
        cd, cs = dsp.cell(i, j), ssp.cell(i, j)
        cls = prob.classify(cd)
        if cls == 'out':
            continue
        vx, vy, exx, eyy, gxy = FD.disp_cell(prob, dsp, A, B, cd)
        txx, tyy, txy, _ = FD.stress_cell(prob, ssp, C, par, cs)
        if cls == 'full':
            v = cd['h'] * cd['k'] * X.quad_form([txx, tyy, txy, exx, eyy, gxy], {(0, 3): 1, (1, 4): 1, (2, 5): 1})
            lo, hi = lo + v, hi + v
        else:
            dens = X.b_add(X.b_mul(txx, exx), X.b_mul(tyy, eyy), X.b_mul(txy, gxy))
            I = RG.integrate_region(dens, cd['x0'], cd['h'], cd['y0'], cd['k'], [('circle', FD.R_HOLE2)], [])
            lo, hi = lo + RG.lo_q(I), hi + RG.hi_q(I)
    return lo, hi


def test_green():
    print('4. Green identity with random admissible fields')
    rng = np.random.default_rng(11)
    br = [Fr(i, 8) for i in range(9)]
    for variant in ('frictionless', 'clamped'):
        meta, arrays = synth.make_P1(rng, br, br, variant=variant)
        good, rep = run(meta, arrays)
        ok(good, 'P1 %s, random fields: all exact checks pass' % variant)
        for vn, tn in (('v', 'tau'), ('vp', 'tm')):
            lo, hi = cross_integral(meta, arrays, vn, tn)
            F = rep.numbers[('F', tn)]
            ok(lo <= F / 2 <= hi and hi - lo < Fr(1, 10 ** 100),
               'P1 %s: int tau:eps(v) in [%.15g, %.15g] encloses F/2 = %.15g' % (variant, float(lo), float(hi),
                                                                                float(F / 2)))
    br2 = [Fr(i, 4) for i in range(-4, 5)]
    meta, arrays = synth.make_P2(rng, br2)
    good, rep = run(meta, arrays)
    ok(good, 'P2, random fields: all exact checks pass')
    for vn, tn in (('v', 'tau'), ('vm', 'tp')):
        lo, hi = cross_integral(meta, arrays, vn, tn)
        F = rep.numbers[('F', tn)]
        ok(lo == hi == F, 'P2: int tau:eps(v) = F(tau) = %.15g exactly' % float(F))


# ------------------------------------------------------------------------------------------------- 5
def test_falsification():
    print('5. falsification: random admissible fields')
    rng = np.random.default_rng(23)
    br = [Fr(i, 6) for i in range(7)]
    worst = []
    for trial in range(3):
        for prob in ('P1', 'P2'):
            if prob == 'P1':
                meta, arrays = synth.make_P1(rng, br, [Fr(i, 6) for i in range(7)], scale=10.0 ** rng.uniform(-3, 0))
            else:
                meta, arrays = synth.make_P2(rng, [Fr(i, 3) for i in range(-3, 4)], scale=10.0 ** rng.uniform(-3, 0))
            good, rep = run(meta, arrays)
            n = rep.numbers
            pr = FD.Problem(meta)
            Ks = []
            for t in ('tau', 'tp', 'tm'):
                for v in ('v', 'vp', 'vm'):
                    Klow = pr.K_factor * pr.work_factor ** 2 * n[('F', t)] ** 2 / n[('Q', t)][1]
                    Kup = pr.K_factor * n[('a', v)][1]
                    Ks.append(Klow <= Kup)
                    for kap in (Fr(1, 10), Fr(-3, 100), Fr(2)):
                        mlow = (pr.work_factor * n[('F', t)] - n[('Q', t)][1] / 2 - kap * n[('l', t)][1 if kap > 0 else 0]
                                - kap * kap / 2 * FD.D.C[0][0] / pr.patch_measure)
                        mup = n[('a', v)][1] / 2 - kap * n[('l', v)][0 if kap > 0 else 1]
                        Ks.append(mlow <= mup)
            lo, hi = rep.bounds['Q2']
            ok(good and all(Ks) and lo <= hi, '%s trial %d: all K and m(k) lower bounds below all upper bounds '
                                              '(%d pairs), Q2 lower <= upper' % (prob, trial, len(Ks)))


def test_near_optimal():
    """Random fields give bounds far apart, which cannot reveal an error that makes bounds too tight.  Here the fields
    are near-optimal (a small float Ritz solver, ritz.py), so such an error would make the bounds cross."""
    print('5b. falsification with near-optimal fields and random admissible perturbations of them')
    rng = np.random.default_rng(31)
    cases = (('P1', 'frictionless', [Fr(i, 8) for i in range(9)], Fr(1, 1000)),
             ('P1', 'clamped', [Fr(i, 8) for i in range(9)], Fr(1, 1000)),
             ('P2', 'standard', [Fr(i, 4) for i in range(-4, 5)], Fr(1, 100)))
    for name, variant, br, gap in cases:
        meta, arrays = ritz.certificate(name, variant, br, br, 3, 4, Fr(1, 100))
        good, rep = run(meta, arrays)
        (klo, khi), (qlo, qhi) = rep.bounds['Q1'], rep.bounds['Q2']
        rel = (khi - klo) / khi
        ok(good and klo <= khi and qlo <= qhi and (rel < gap if name == 'P1' else rel < Fr(1, 10)),
           '%s %s: near-optimal fields, K in [%.9f, %.9f] (relative gap %.1e), Q2 in [%.6f, %.6f]'
           % (name, variant, float(klo), float(khi), float(rel), float(qlo), float(qhi)))
        for eps in (Fr(1, 10 ** 8), Fr(1, 10 ** 5), Fr(1, 100)):
            if name == 'P1':
                _, R = synth.make_P1(rng, br, br, pu=3, ps=4, gexp=(0, 1, 2, 3), variant=variant)
            else:
                _, R = synth.make_P2(rng, br, pu=3, ps=4)
            a2 = {}
            for key, val in arrays.items():
                scale = max(abs(v) for v in val.flat) or Fr(1)
                a2[key] = val + R[key] * (eps * scale)
            good2, rep2 = run(meta, a2)
            cross = []
            for q in ('Q1', 'Q2'):
                (l1, h1), (l2, h2) = rep.bounds[q], rep2.bounds[q]
                cross.append(max(l1, l2) <= min(h1, h2))
            ok(good2 and all(cross), '%s %s: fields moved by %.0e (relative): all checks pass, K and Q2 lower bounds '
                                     'stay below the upper bounds of the unmoved fields and vice versa'
               % (name, variant, float(eps)))


# ------------------------------------------------------------------------------------------------- 6
def expect_fail(meta, arrays, needle, what):
    good, rep = run(meta, arrays)
    hit = [f for f in rep.fail if needle in f]
    ok((not good) and bool(hit), 'rejected: %s  (failed check: %s)' % (what, hit[0][:90] if hit else rep.fail[:1]))


def test_negative():
    print('6. negative tests: broken certificates are rejected')
    rng = np.random.default_rng(3)
    br = [Fr(i, 4) for i in range(5)]
    meta, arrays = synth.make_P1(rng, br, br)
    a2 = dict(arrays)
    C = arrays['tau_C'].copy()
    C[0, 2] += Fr(1, 8)
    a2['tau_C'] = C
    expect_fail(meta, a2, 'zero shear on x = 0', 'P1 stress with shear on the mirror line x = 0')
    a2 = dict(arrays)
    A = arrays['v_A'].copy()
    A[-1, 1] = Fr(1, 64)
    a2['v_A'] = A
    expect_fail(meta, a2, 'v_x = 1/2 on x = 1', 'P1 displacement with a wrong grip value')
    # C^0 stress function: knot 3/4 of multiplicity ps in x
    m2 = synth.make_P1(rng, br, br, mult_s={Fr(3, 4): 3})
    expect_fail(*m2, 'normal traction continuous', 'P1 stress function only C^0 across x = 3/4')
    # discontinuous displacement: knot 1/4 of multiplicity pu+1
    meta3, arr3 = synth.make_P1(rng, br, br)
    sp = dict(meta3['displacement_space'])
    kx = [Fr(v) for v in sp['knots_x']]
    kx.insert(kx.index(Fr(1, 4)), Fr(1, 4))
    kx.insert(kx.index(Fr(1, 4)), Fr(1, 4))
    sp['knots_x'] = [str(v) for v in kx]
    meta3 = dict(meta3, displacement_space=sp)
    nx = len(kx) - sp['degree'] - 1
    for nm in ('v', 'vp', 'vm'):
        for c in 'AB':
            M = synth.rnd(rng, (nx, arr3[nm + '_' + c].shape[1]), 0.1)
            if c == 'A':
                M[0, :] = M[-1, :] = Fr(0)
            else:
                M[:, 0] = Fr(0)
            arr3[nm + '_' + c] = M
    expect_fail(meta3, arr3, 'continuous across', 'P1 displacement discontinuous across x = 1/4')
    # stress not traction-free on the hole: drop the weight rho^2 (checker variant, same data)
    orig = FD.Problem.__init__

    def no_weight(self, m):
        orig(self, m)
        self.weight = False
    FD.Problem.__init__ = no_weight
    try:
        expect_fail(meta, arrays, 'divisible', 'P1 stress function without the weight rho^2 (traction on the hole)')
    finally:
        FD.Problem.__init__ = orig
    # claims tighter than the recomputed bounds
    meta4, arr4 = synth.make_nohole()
    meta4 = dict(meta4, claims={'Q1': {'lower': '1.0000001', 'upper': '1.1'}, 'Q2': {'lower': '0', 'upper': '1'}})
    expect_fail(meta4, arr4, 'Q1 claim', 'claimed lower bound above the recomputed one')
    meta4 = dict(meta4, claims={'Q1': {'lower': '0.99', 'upper': '1.01'}, 'Q2': {'lower': '0.4', 'upper': '0.6'}})
    expect_fail(meta4, arr4, 'Q2 claim', 'claimed Q2 interval narrower than the recomputed one')
    # P2: a stress coefficient that puts traction on the free edge x = 0, y <= 0
    meta5, arr5 = synth.make_P2(rng, [Fr(i, 4) for i in range(-4, 5)])
    kxs = [Fr(v) for v in meta5['stress_space']['knots_x']]
    Lx, Rx, Sx = synth.p2_classes(kxs, meta5['stress_space']['degree'])
    Ly, Ry, Sy = synth.p2_classes(kxs, meta5['stress_space']['degree'])
    C = arr5['tau_C'].copy()
    C[Sx[0], Ly[-1]] = Fr(1, 4)
    arr5['tau_C'] = C
    expect_fail(meta5, arr5, 'traction on y = 1, y = -1', 'P2 stress with traction on the free edge x = 0')


# ------------------------------------------------------------------------------------------------- 7 (P3)
def run3(meta, arrays, weight=True):
    import io
    import contextlib
    import p3check
    rep = p3check.Rep()
    with contextlib.redirect_stdout(io.StringIO()):
        try:
            good = p3check.check_cert(meta, arrays, rep, procs=1, weight=weight, claims_required=False)
        except AssertionError as exc:
            rep.fail.append('assertion: %s' % (exc,))
            good = False
    return good, rep


def cross3(meta, arrays):
    """int over the unit octant cube of tau:eps(v) (no cavity; every cell whole), exact"""
    import exact3 as E
    import p3check
    cert = p3check.Cert(meta, arrays)
    N = cert.N
    tot = Fr(0)
    w = {(p, p + 6): Fr(1) for p in range(6)}
    for i in range(N):
        for j in range(N):
            for k in range(N):
                tau, _ = p3check.stress_cell(cert, (i, j, k))
                _, eps = p3check.disp_cell(cert, (i, j, k))
                tot += E.quad_form3([tau[c] for c in p3check.COMP] + eps, w) / N ** 3
    return tot


def test_p3():
    print('7. P3 (3D, octant): known answer, Green identity, box integrals, falsification, negative tests')
    import exact3 as E
    import p3check
    meta, arrays = synth.make_P3_solid()
    good, rep = run3(meta, arrays)
    ok(good and rep.bounds == (2, 2), 'solid cube, frictionless grips: all checks pass, K = [2, 2] exactly (known answer)')
    rng = np.random.default_rng(41)
    for variant in ('frictionless', 'clamped'):
        meta, arrays = synth.make_P3(rng, variant=variant, cavity=False)
        good, rep = run3(meta, arrays)
        cert = p3check.Cert(meta, arrays)
        F = sum(((g + d) * (1 - cert.gbase) ** (e + 1) / (e + 1) for e, g, d in zip(cert.gexp, cert.gamma, cert.delta)),
                Fr(0))
        val = cross3(meta, arrays)
        ok(good and val == F / 2, 'solid cube, %s grips, random admissible fields: all checks pass and int tau:eps(v) = '
                                  'F_o/2 = %.12g exactly' % (variant, float(F / 2)))
    # superset box integrals against direct integration on each box
    meta, arrays = synth.make_P3(rng)
    cert = p3check.Cert(meta, arrays)
    e = next(c for c in ((0, 0, 1), (0, 1, 1), (1, 1, 0)) if p3check.cell_class(cert, c) == 'cut')
    tau, _ = p3check.stress_cell(cert, e)
    T = [tau[c] for c in p3check.COMP]
    f = E.t_add(*[E.t_mul(T[p], T[q]) for p in range(6) for q in range(6)])
    boxes = p3check.kept_boxes(cert, e, 3)
    direct = sum((E.t_integral(f, tuple((Fr(j, 2 ** l), Fr(j + 1, 2 ** l)) for j in jj))
                  for (l, *jj) in boxes), Fr(0))
    fi, den = E.density_int(T, {(p, q): Fr(1) for p in range(6) for q in range(6)})
    ok(p3check.box_integrals(f, boxes) == direct == p3check.box_integrals(fi, boxes, den) and len(boxes) > 8,
       'octree boxes of a cut cell (%d boxes, depth 3): integer box integrals (Fraction and integer densities) equal '
       'direct exact integrals' % len(boxes))
    for trial in range(2):
        meta, arrays = synth.make_P3(rng, variant=('frictionless', 'clamped')[trial], scale=10.0 ** rng.uniform(-3, 0))
        good, rep = run3(meta, arrays)
        ok(good and rep.bounds[0] <= rep.bounds[1], 'random admissible fields with the cavity (%s): all checks pass, '
                                                    'lower bound below upper bound' % meta['variant'])
    # negative tests
    meta, arrays = synth.make_P3(rng)
    a2 = dict(arrays)
    P = arrays['Psi'].copy()
    P[3][0, 2, 2] = Fr(1, 8)               # xy component is odd in x: must vanish at index 0 (shear on x = 0)
    a2['Psi'] = P
    good, rep = run3(meta, a2)
    ok(not good and any('mirror' in f for f in rep.fail), 'rejected: P3 potential that puts shear on the mirror x = 0')
    good, rep = run3(meta, arrays, weight=False)
    ok(not good and any('sphere' in f for f in rep.fail), 'rejected: P3 stress without the weight rho^2 (traction on '
                                                         'the sphere)')
    meta3, arr3 = synth.make_P3(np.random.default_rng(5))
    ss = dict(meta3['stress_space'])
    kx = [Fr(v) for v in ss['knots'][0]]
    kx.insert(kx.index(Fr(1, 4)), Fr(1, 4))
    kx.insert(kx.index(Fr(1, 4)), Fr(1, 4))
    ss['knots'] = [[str(v) for v in kx]] + ss['knots'][1:]
    meta3 = dict(meta3, stress_space=ss)
    P = synth.p3_conditions(synth.rnd(np.random.default_rng(6), (6, len(kx) - 4) + arr3['Psi'].shape[2:], 0.3), True)
    arr3 = dict(arr3, Psi=P)
    good, rep = run3(meta3, arr3)
    ok(not good and any('traction continuous' in f for f in rep.fail), 'rejected: P3 potential only C^0 across x = 1/4')
    meta4, arr4 = synth.make_P3_solid()
    meta4 = dict(meta4, claims={'Q1': {'lower': '2.000001', 'upper': '2.1'}})
    good, rep = run3(meta4, arr4)
    ok(not good and any('claim' in f for f in rep.fail), 'rejected: claimed lower bound above the recomputed one')


def gamma_half(m2):
    """Gamma(m2/2) = rational * sqrt(pi)^k"""
    from math import factorial
    if m2 % 2 == 0:
        return Fr(factorial(m2 // 2 - 1)), 0
    n = (m2 - 1) // 2
    return Fr(factorial(2 * n), 4 ** n * factorial(n)), 1


def octant_ball_moment(a, b, c, R):
    """int over {x, y, z >= 0, |x| <= R} of x^a y^b z^c  =  R^(s)/s * Gamma(A)Gamma(B)Gamma(C) / (4 Gamma(A+B+C)),
    s = a+b+c+3, A = (a+1)/2, ...  (an interval)"""
    num, p = Fr(1), 0
    for m in (a + 1, b + 1, c + 1):
        g, q = gamma_half(m)
        num, p = num * g, p + q
    gd, qd = gamma_half(a + b + c + 3)
    s = a + b + c + 3
    return RG.ivq(num / gd / 4 * R ** s / s) * iv.sqrt(iv.pi) ** (p - qd)


def test_p3_exact():
    print('8. P3 exact integration over cells cut by the sphere (ball.py)')
    import ball
    import exact3 as E
    import p3check
    ok(all(ball.check_Pn(n) for n in range(24)), 'closed forms int y^n sqrt(c^2 - y^2) dy: exact identities, n < 24')
    R2 = Fr(1, 4)
    for (a, b, c) in ((0, 0, 0), (2, 0, 0), (1, 1, 1), (3, 2, 0), (2, 3, 5)):
        f = E.zeros3((a + 1, b + 1, c + 1))
        f[a, b, c] = Fr(1)
        N = 4
        lo = hi = Fr(0)
        for e in [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]:
            rmin, rmax = sum(v * v for v in e), sum((v + 1) ** 2 for v in e)
            if rmin * 4 >= N * N:
                continue
            box = tuple((Fr(v, N), Fr(v + 1, N)) for v in e)
            if rmax * 4 <= N * N:
                v = E.t_integral(f, box)
                lo, hi = lo + v, hi + v
            else:
                I = ball.ball_part(f, box, R2)
                lo, hi = lo + RG.lo_q(I), hi + RG.hi_q(I)
        ex = octant_ball_moment(a, b, c, Fr(1, 2))
        ok(RG.ivq(lo).a <= ex.a and ex.b <= RG.ivq(hi).b and hi - lo < Fr(1, 10 ** 150),
           'int x^%d y^%d z^%d over the octant of the ball (grid 4, cut cells in closed form) encloses the Gamma-function '
           'value (width %.0e)' % (a, b, c, float(hi - lo)))
    # known answer: the lift alone, v = (x/2) e_x, has a_o = (lambda + 2 mu) |Omega_o| / 4, |Omega_o| = 1 - pi/48
    rng = np.random.default_rng(43)
    meta, arrays = synth.make_P3(rng, N=4, pu=2, ps=3)
    meta = dict(meta, integration={'mode': 'exact'})
    arr0 = dict(arrays, V=arrays['V'] * 0)
    good, rep = run3(meta, arr0)
    lam, mu = FD.D.P3_LAME['lam'], FD.D.P3_LAME['mu']
    exact_a = RG.ivq((lam + 2 * mu) / 4) * (1 - iv.pi / 48)
    Alo, Ahi = rep.energies[2], rep.energies[3]
    ok(good and RG.ivq(Alo).a <= exact_a.a and exact_a.b <= RG.ivq(Ahi).b and Ahi - Alo < Fr(1, 10 ** 100),
       'exact mode, the lift alone: a_o encloses (lambda + 2 mu)(1 - pi/48)/4 (width %.0e)' % float(Ahi - Alo))
    # Green identity over the octant minus the ball, random admissible fields, exact mode
    for variant in ('frictionless', 'clamped'):
        meta, arrays = synth.make_P3(rng, N=4, pu=2, ps=3, variant=variant)
        meta = dict(meta, integration={'mode': 'exact'})
        good, rep = run3(meta, arrays)
        cert = p3check.Cert(meta, arrays)
        lo = hi = Fr(0)
        w = {(p, p + 6): Fr(1) for p in range(6)}
        N = cert.N
        for e in [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]:
            cls = p3check.cell_class(cert, e)
            if cls == 'out':
                continue
            tau, _ = p3check.stress_cell(cert, e)
            _, eps = p3check.disp_cell(cert, e)
            polys = [tau[c] for c in p3check.COMP] + eps
            if cls == 'full':
                v = E.quad_form3(polys, w) / N ** 3
                lo, hi = lo + v, hi + v
            else:
                a_, b_ = p3check.cut_exact(cert, e, polys, w)
                lo, hi = lo + a_, hi + b_
        F = rep.energies[4]
        ok(good and lo <= F / 2 <= hi and hi - lo < Fr(1, 10 ** 100),
           'octant minus the ball, %s grips, random admissible fields: all checks pass and int tau:eps(v) encloses '
           'F_o/2 = %.12g (width %.0e)' % (variant, float(F / 2), float(hi - lo)))


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    t0 = time.time()
    test_bspline()
    test_region()
    test_known_answer()
    test_green()
    test_falsification()
    test_near_optimal()
    test_negative()
    test_p3()
    test_p3_exact()
    print('%s   (%.0f s)' % ('ALL TESTS PASSED' if not FAILS else '%d TEST(S) FAILED: %s' % (len(FAILS), FAILS),
                            time.time() - t0))
    sys.exit(1 if FAILS else 0)
