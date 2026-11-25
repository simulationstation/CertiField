"""Independent checker of the benchmark-pack certificates P1 and P2 (format benchpack-certificate-1).

usage:  python check.py <path/to/STEM> [--procs N] [--log FILE]
        (STEM.json and STEM.npz are read; the extension may be given or not)

What it does (details in README.md):
 1. rebuilds every displacement and stress field as one polynomial with rational coefficients per cell;
 2. checks in exact rational arithmetic: equilibrium in every cell; continuity of the normal traction across every
    interior cell edge; zero traction on the hole (tau.grad(rho) divisible by both factors of rho) and on the free
    edges; zero shear on the mirror lines and frictionless grips; the grip force; continuity of the displacement
    and its grip and mirror values; the coefficient conditions of FORMAT.md; that the cells cover the domain;
 3. recomputes a, Q, F and the patch functional l of every field: exact rationals on every cell inside the domain,
    closed forms evaluated with mpmath.iv (outward rounding) on cells cut by the hole and on the patch;
 4. forms the bounds of FORMAT.md 3.4 / 4.4 and confirms lower_claim <= recomputed lower and recomputed upper <=
    upper_claim.
Exit code 0 if every check passes.
"""
import argparse
import multiprocessing as mp_
import os
import sys
import time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np            # noqa: E402
from mpmath import iv        # noqa: E402

import exact as X            # noqa: E402
import fields as FD          # noqa: E402
import region as RG          # noqa: E402

D = FD.D


class Report:
    def __init__(self, logfile=None):
        self.lines = []
        self.fail = []
        self.logfile = logfile
        self.t0 = time.time()

    def say(self, s=''):
        line = s
        print(line, flush=True)
        self.lines.append(line)
        if self.logfile:
            with open(self.logfile, 'a', encoding='utf-8') as f:
                f.write(line + '\n')

    def check(self, ok, what):
        self.say(('  [pass] ' if ok else '  [FAIL] ') + what)
        if not ok:
            self.fail.append(what)
        return ok


def fmt(q, n=12):
    return ('%.' + str(n) + 'g') % float(q) if not isinstance(q, str) else q


def ivfmt(I, n=13):
    return '[%s, %s]' % (fmt(RG.lo_q(I), n), fmt(RG.hi_q(I), n))


def chunks(cells, nchunk):
    return [cells[i::nchunk] for i in range(nchunk)]


def run_field(pool, kind, name, cells, nchunk):
    fn = FD.work_disp if kind == 'disp' else FD.work_stress
    jobs = [(name, ch) for ch in chunks(cells, nchunk) if ch]
    out = {}
    if pool is None:
        res = map(fn, jobs)
    else:
        res = pool.imap_unordered(fn, jobs)
    for lst in res:
        for rec in lst:
            out[rec['cell']] = rec
    return out


def total(recs):
    """interval of the sum of 'energy' (exact) and 'energy_iv' (pairs) over the records"""
    ex = Fr(0)
    lo = hi = Fr(0)
    for r in recs.values():
        if 'energy' in r:
            ex += r['energy']
        if 'energy_iv' in r:
            lo += r['energy_iv'][0]
            hi += r['energy_iv'][1]
    return ex, (ex + lo, ex + hi)


def patch_total(recs):
    lo = hi = Fr(0)
    for r in recs.values():
        p = r.get('patch')
        if p is None:
            continue
        if isinstance(p, tuple):
            lo += p[0]
            hi += p[1]
        else:
            lo += p
            hi += p
    return lo, hi


def ivp(pair):
    lo, hi = pair
    return iv.mpf([RG.ivq(lo).a, RG.ivq(hi).b])


def coefficient_conditions(prob, meta, arrays, dnames, snames):
    """the stored-coefficient conditions of FORMAT.md (sufficient for the identities, which are checked separately)"""
    res = []
    for n in dnames:
        A, B = arrays[n + '_A'], arrays[n + '_B']
        if prob.name == 'P1':
            ok = all(v == 0 for v in A[0, :]) and all(v == 0 for v in A[-1, :]) and all(v == 0 for v in B[:, 0])
            if not prob.frictionless:
                ok = ok and all(v == 0 for v in B[-1, :])
        else:
            ok = all(v == 0 for M in (A, B) for v in list(M[0, :]) + list(M[-1, :]))
        res.append(('displacement %s: grip/mirror coefficient conditions' % n, ok))
    for n in snames:
        C = arrays[n + '_C']
        if prob.name == 'P1':
            ok = (all(C[0, :] == C[1, :]) and all(C[:, 0] == C[:, 1]) and all(v == 0 for v in C[:, -1])
                  and all(v == 0 for v in C[:, -2]))
            if prob.frictionless:
                ok = ok and all(C[-1, :] == C[-2, :])
        else:
            ok = True
        res.append(('stress %s: FORMAT coefficient conditions' % n, bool(ok)))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('stem')
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--log', default=None)
    a = ap.parse_args()
    stem = a.stem
    for ext in ('.json', '.npz'):
        if stem.endswith(ext):
            stem = stem[:-len(ext)]
    meta, arrays, raw = FD.load(stem)
    rep = Report(a.log)
    rep.say('checker: %s' % os.path.basename(stem))
    ok = check_cert(meta, arrays, rep, procs=a.procs)
    rep.say('RESULT: %s   (%.0f s)' % ('ALL CHECKS PASSED' if ok else 'FAILED: %d check(s)' % len(rep.fail),
                                       time.time() - rep.t0))
    sys.exit(0 if ok else 1)


def check_cert(meta, arrays, rep, procs=8, claims_required=True):
    """run every check on one certificate (meta: the JSON dict; arrays: exact Fraction arrays). Returns True if all
    pass; the recomputed numbers are stored in rep.numbers."""
    rep.numbers = {}
    prob = FD.Problem(meta)
    rep.say('problem %s, variant %s, model %s, label %r' % (prob.name, prob.variant, meta.get('model'),
                                                            meta.get('label')))
    try:
        FD.assert_independent()
        indep = True
    except AssertionError:
        indep = False
    rep.check(indep, 'no module of benchpack/builder is loaded (only definitions.py is shared)')
    rep.check(meta.get('format') == 'benchpack-certificate-1', 'format string benchpack-certificate-1')
    if prob.name == 'P1':
        rep.check(meta.get('model') == 'quarter', 'P1 model is the quarter model')
        rep.check(FD.check_rho_factors(), 'rho of definitions.py = (x^2+y^2-1/4)(7/4+y^2-x^2) exactly')
    dsp = FD.Space(meta['displacement_space'])
    ssp = FD.Space(meta['stress_space'])
    lo_b, hi_b = prob.box
    # ---- knots
    rep.say('spaces: displacement degree %d, %d x %d cells; stress degree %d, %d x %d cells'
            % (dsp.p, len(dsp.sx), len(dsp.sy), ssp.p, len(ssp.sx), len(ssp.sy)))
    mult_ok = []
    for nm, sp, need in (('displacement', dsp, dsp.p), ('stress', ssp, ssp.p - 1)):
        for ax, kn in (('x', sp.kx), ('y', sp.ky)):
            ok, mult, msg = X.knot_report(kn, sp.p, lo_b, hi_b)
            rep.check(ok, '%s knots_%s open on [%s, %s] %s' % (nm, ax, lo_b, hi_b, msg))
            mult_ok.append(rep.check(mult <= need, '%s knots_%s: interior multiplicity %d <= %d (C^%d splines)'
                                     % (nm, ax, mult, need, sp.p - need)))
            if prob.name == 'P2':
                rep.check(Fr(0) in kn, '%s knots_%s contain 0' % (nm, ax))
    if prob.name == 'P1' and ssp.g_exp:
        rep.check(ssp.g_base in ssp.ky, 'g_base %s is a knot of the stress knots_y' % ssp.g_base)
        if prob.hole:
            rep.check(ssp.g_base >= FD.D.P1['hole_radius'],
                      'g vanishes below y = g_base >= 1/2, so g is zero on cells meeting the hole')
    # ---- fields
    flds = meta['fields']
    q1 = flds['Q1']
    q2 = flds.get('Q2')
    dnames = [q1['v']] + ([q2['v_plus'], q2['v_minus']] if q2 else [])
    snames = [q1['tau']] + ([q2['tau_plus'], q2['tau_minus']] if q2 else [])
    dnames = list(dict.fromkeys(dnames))
    snames = list(dict.fromkeys(snames))
    for n in dnames:
        rep.check(arrays[n + '_A'].shape == (dsp.nx, dsp.ny) and arrays[n + '_B'].shape == (dsp.nx, dsp.ny),
                  'displacement %s: array shapes (%d, %d)' % (n, dsp.nx, dsp.ny))
    for n in snames:
        rep.check(arrays[n + '_C'].shape == (ssp.nx, ssp.ny), 'stress %s: array shape (%d, %d)' % (n, ssp.nx, ssp.ny))
        if prob.name == 'P1':
            npar = len(arrays.get(n + '_par', []))
            rep.check(npar == len(ssp.g_exp), 'stress %s: %d g coefficients for %d exponents' % (n, npar, len(ssp.g_exp)))
        else:
            rep.check(len(arrays[n + '_par']) == 3, 'stress %s: (alpha, beta, gamma) given' % n)
    if len(rep.fail) > mult_ok.count(False):     # structural errors: the fields cannot be rebuilt
        return False
    for what, ok in coefficient_conditions(prob, meta, arrays, dnames, snames):
        rep.check(ok, what)     # sufficient conditions only; the identities themselves are checked below

    pool = mp_.Pool(procs, initializer=FD.init_worker, initargs=(meta, arrays)) if procs > 1 else None
    if pool is None:
        FD.init_worker(meta, arrays)
    nchunk = max(1, 4 * procs)
    try:
        dres = {}
        sres = {}
        for n in dnames:
            t = time.time()
            dres[n] = run_field(pool, 'disp', n, dsp.cells(), nchunk)
            rep.say('displacement %s: %d cells evaluated (%.0f s)' % (n, len(dres[n]), time.time() - t))
        for n in snames:
            t = time.time()
            sres[n] = run_field(pool, 'stress', n, ssp.cells(), nchunk)
            rep.say('stress %s: %d cells evaluated (%.0f s)' % (n, len(sres[n]), time.time() - t))
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    # ---- exact identities: displacement
    for n in dnames:
        disp_checks(rep, prob, dsp, n, dres[n])
    # ---- exact identities: stress
    forces = {}
    for n in snames:
        par = list(arrays[n + '_par']) if (n + '_par') in arrays else []
        forces[n] = stress_checks(rep, prob, ssp, n, sres[n], par)
    # ---- coverage
    coverage(rep, prob, ssp, dsp)

    # ---- bounds
    meas = prob.patch_measure
    c11 = D.C[0][0]
    num = {}
    for n in dnames:
        ex, (lo, hi) = total(dres[n])
        num[('a', n)] = (lo, hi)
        num[('l', n)] = tuple(v / meas for v in patch_total(dres[n]))
    for n in snames:
        ex, (lo, hi) = total(sres[n])
        num[('Q', n)] = (lo, hi)
        num[('l', n)] = tuple(v / meas for v in patch_total(sres[n]))
        num[('F', n)] = forces[n]
    rep.numbers = num
    rep.say('recomputed integrals (exact on cells inside the domain; intervals where curved):')
    for key in sorted(num, key=str):
        v = num[key]
        if isinstance(v, tuple):
            rep.say('  %s(%s) in [%.15g, %.15g]  width %.2e' % (key[0], key[1], float(v[0]), float(v[1]),
                                                               float(v[1] - v[0])))
        else:
            rep.say('  %s(%s) = %.15g (exact)' % (key[0], key[1], float(v)))
    K = prob.K_factor
    wf = prob.work_factor
    F = forces[q1['tau']]
    Qi = ivp(num[('Q', q1['tau'])])
    ai = ivp(num[('a', q1['v'])])
    Kup = (RG.ivq(K) * ai).b
    # K >= K_factor w^2 F^2/Q:  max_l [l w F - l^2 Q/2] = w^2 F^2/(2Q) <= (1/2) a(u,u) = K/(2 K_factor)
    # (w = work factor: P1 quarter w = 1/2, so 4 (1/4) F^2/Q = F^2/Q; P2 w = 1)
    Klow = (RG.ivq(K * wf * wf) * RG.ivq(F) ** 2 / Qi).a
    bounds = {'Q1': (Klow, Kup)}
    if q2:
        kap = Fr(q2['kappa'])
        rep.check(kap > 0, 'kappa = %s > 0' % q2['kappa'])
        kv = RG.ivq(kap)
        cst = RG.ivq(kap * kap / 2 * c11 / meas)
        vp, vm, tp, tm = q2['v_plus'], q2['v_minus'], q2['tau_plus'], q2['tau_minus']
        m_up_p = ivp(num[('a', vp)]) / 2 - kv * ivp(num[('l', vp)])
        m_up_m = ivp(num[('a', vm)]) / 2 + kv * ivp(num[('l', vm)])
        m_lo_p = RG.ivq(wf * forces[tp]) - ivp(num[('Q', tp)]) / 2 - kv * ivp(num[('l', tp)]) - cst
        m_lo_m = RG.ivq(wf * forces[tm]) - ivp(num[('Q', tm)]) / 2 + kv * ivp(num[('l', tm)]) - cst
        rep.say('  m_up(+k) <= %.15g, m_up(-k) <= %.15g, m_low(+k) >= %.15g, m_low(-k) >= %.15g'
                % (float(RG.hi_q(m_up_p)), float(RG.hi_q(m_up_m)), float(RG.lo_q(m_lo_p)), float(RG.lo_q(m_lo_m))))
        q2lo = ((m_lo_m - m_up_p) / (2 * kv)).a
        q2hi = ((m_up_m - m_lo_p) / (2 * kv)).b
        bounds['Q2'] = (q2lo, q2hi)
    rep.bounds = {k: (RG.lo_q(v[0]), RG.hi_q(v[1])) for k, v in bounds.items()}
    claims = meta.get('claims', {})
    for qn, (lo, hi) in rep.bounds.items():
        rel = (hi - lo) / abs((hi + lo) / 2) if hi + lo else Fr(0)
        rep.say('recomputed %s: [%.12g, %.12g]  width %.3e  relative %.3e' % (qn, float(lo), float(hi),
                                                                            float(hi - lo), float(rel)))
        cl = claims.get(qn)
        if cl is None:
            if claims_required:
                rep.check(False, 'claim for %s present' % qn)
            continue
        clo, chi = Fr(cl['lower']), Fr(cl['upper'])
        rep.check(clo <= lo and hi <= chi, '%s claim [%s, %s] confirmed: %s <= %.12g and %.12g <= %s'
                  % (qn, cl['lower'], cl['upper'], cl['lower'], float(lo), float(hi), cl['upper']))
    comp = meta.get('computed')
    if comp:
        rep.say('builder\'s own enclosures (diagnostic, not used): %s' % short(comp))
    return not rep.fail


def short(obj, n=1500):
    s = str(obj)
    return s if len(s) <= n else s[:n] + ' ...'


def neighbours(sp):
    nx, ny = len(sp.sx), len(sp.sy)
    for i in range(nx):
        for j in range(ny):
            if i + 1 < nx:
                yield (i, j), (i + 1, j), 'v'
            if j + 1 < ny:
                yield (i, j), (i, j + 1), 'h'


def edge_counts(prob, sp, recs, a, b, kind):
    """True if the edge between cells a and b is part of the domain's interior"""
    ra, rb = recs[a], recs[b]
    if prob.name == 'P2':
        return ra['cls'] != 'out' and rb['cls'] != 'out'
    ca = sp.cell(*a)
    if kind == 'v':
        return prob.edge_in_domain(ca['x1'], ca['y0'], ca['x1'], ca['y1'])
    return prob.edge_in_domain(ca['x0'], ca['y1'], ca['x1'], ca['y1'])


def disp_checks(rep, prob, sp, n, recs):
    nx, ny = len(sp.sx), len(sp.sy)
    bad = 0
    cnt = 0
    for a, b, kind in neighbours(sp):
        if not edge_counts(prob, sp, recs, a, b, kind):
            continue
        cnt += 1
        ta, tb = recs[a]['tr'], recs[b]['tr']
        pa, pb = (ta['R'], tb['L']) if kind == 'v' else (ta['T'], tb['B'])
        if not all(X.same_poly(u, w) for u, w in zip(pa, pb)):
            bad += 1
    rep.check(bad == 0, 'displacement %s: continuous across all %d interior edges (exact)' % (n, cnt))
    half = Fr(1, 2)
    errs = []
    for (i, j), r in recs.items():
        c = sp.cell(i, j)
        tr = r['tr']
        if prob.name == 'P1':
            if i == 0 and prob.edge_in_domain(c['x0'], c['y0'], c['x0'], c['y1']):
                if not X.is_zero_poly(tr['L'][0]):
                    errs.append(('v_x = 0 on x = 0', (i, j)))
            if i == nx - 1:
                if not X.same_poly(tr['R'][0], [half]):
                    errs.append(('v_x = 1/2 on x = 1', (i, j)))
                if not prob.frictionless and not X.is_zero_poly(tr['R'][1]):
                    errs.append(('v_y = 0 on x = 1 (clamped)', (i, j)))
            if j == 0 and prob.edge_in_domain(c['x0'], c['y0'], c['x1'], c['y0']):
                if not X.is_zero_poly(tr['B'][1]):
                    errs.append(('v_y = 0 on y = 0', (i, j)))
        else:
            if i == 0 and not (X.is_zero_poly(tr['L'][0]) and X.is_zero_poly(tr['L'][1])):
                errs.append(('v = 0 on x = -1', (i, j)))
            if i == nx - 1 and c['y0'] >= 0:
                if not (X.is_zero_poly(tr['R'][0]) and X.same_poly(tr['R'][1], [Fr(-1)])):
                    errs.append(('v = (0, -1) on the grip x = 1', (i, j)))
    what = ('v_x = 0 on x = 0, v_x = 1/2 on x = 1%s, v_y = 0 on y = 0'
            % ('' if prob.frictionless else ', v_y = 0 on x = 1')) if prob.name == 'P1' else \
        'v = 0 on x = -1 and v = (0, -1) on x = 1, 0 <= y <= 1'
    rep.check(not errs, 'displacement %s: %s, exactly%s' % (n, what, '' if not errs else ' -- %s' % errs[:5]))


def stress_checks(rep, prob, sp, n, recs, par):
    nx, ny = len(sp.sx), len(sp.sy)
    dom = [c for c, r in recs.items() if r['cls'] != 'out']
    rep.check(all(recs[c]['div0'] for c in dom), 'stress %s: div tau = 0 exactly in all %d cells of the domain'
              % (n, len(dom)))
    bad, cnt = [], 0
    for a, b, kind in neighbours(sp):
        if not edge_counts(prob, sp, recs, a, b, kind):
            continue
        cnt += 1
        ta, tb = recs[a]['tr'], recs[b]['tr']
        pa, pb = (ta['R'], tb['L']) if kind == 'v' else (ta['T'], tb['B'])
        if not all(X.same_poly(u, w) for u, w in zip(pa, pb)):
            bad.append((a, b))
    rep.check(not bad, 'stress %s: normal traction continuous across all %d interior edges (exact)%s'
              % (n, cnt, '' if not bad else ' -- fails at %s' % bad[:5]))
    if prob.name == 'P1' and prob.hole:
        cut = [c for c, r in recs.items() if r['cls'] == 'cut']
        okh = all(recs[c]['hole_div'] for c in cut)
        rep.check(okh, 'stress %s: on all %d cells cut by the hole, both components of tau.grad(rho) are divisible '
                       'by x^2+y^2-1/4 and by 7/4+y^2-x^2, i.e. by rho (zero traction on the hole, exact)'
                  % (n, len(cut)))
    errs = []
    for (i, j), r in recs.items():
        if r['cls'] == 'out':
            continue
        c = sp.cell(i, j)
        tr = r['tr']
        z = X.is_zero_poly
        if prob.name == 'P1':
            if j == ny - 1 and not (z(tr['T'][0]) and z(tr['T'][1])):
                errs.append(('traction on y = 1', (i, j)))
            if i == nx - 1 and prob.frictionless and not z(tr['R'][1]):
                errs.append(('shear on the grip x = 1', (i, j)))
            if i == 0 and prob.edge_in_domain(c['x0'], c['y0'], c['x0'], c['y1']) and not z(tr['L'][1]):
                errs.append(('shear on x = 0', (i, j)))
            if j == 0 and prob.edge_in_domain(c['x0'], c['y0'], c['x1'], c['y0']) and not z(tr['B'][0]):
                errs.append(('shear on y = 0', (i, j)))
        else:
            if j == ny - 1 and not (z(tr['T'][0]) and z(tr['T'][1])):
                errs.append(('traction on y = 1', (i, j)))
            if j == 0 and c['x1'] <= 0 and not (z(tr['B'][0]) and z(tr['B'][1])):
                errs.append(('traction on y = -1', (i, j)))
            if c['x1'] == 0 and c['y1'] <= 0 and not (z(tr['R'][0]) and z(tr['R'][1])):
                errs.append(('traction on x = 0, y <= 0', (i, j)))
            if c['y0'] == 0 and c['x0'] >= 0 and not (z(tr['B'][0]) and z(tr['B'][1])):
                errs.append(('traction on y = 0, x >= 0', (i, j)))
    what = ('zero traction on y = 1; zero shear on x = 0, y = 0%s' % (' and the grip x = 1' if prob.frictionless else '')
            if prob.name == 'P1' else 'zero traction on y = 1, y = -1 (x <= 0), x = 0 (y <= 0), y = 0 (x >= 0)')
    rep.check(not errs, 'stress %s: %s, exactly%s' % (n, what, '' if not errs else ' -- %s' % errs[:5]))
    # force through the grip, exact
    Fsum = Fr(0)
    for (i, j), r in recs.items():
        if i != nx - 1:
            continue
        c = sp.cell(i, j)
        if prob.name == 'P1':
            Fsum += c['k'] * sum((v / (b + 1) for b, v in enumerate(r['tr']['R'][0])), Fr(0))
        elif c['y0'] >= 0:
            Fsum -= c['k'] * sum((v / (b + 1) for b, v in enumerate(r['tr']['R'][1])), Fr(0))
    if prob.name == 'P1':
        Fformula = FD.g_force(sp, par) if sp.g_exp else Fr(0)
        rep.check(Fsum == Fformula, 'stress %s: grip force int_0^1 tau_xx(1,y) dy = %.15g exactly, equal to sum_k '
                                    'gamma_k (1-g_base)^(e_k+1)/(e_k+1)' % (n, float(Fsum)))
    else:
        rep.check(Fsum == -par[1], 'stress %s: grip force -int_0^1 tau_xy(1,y) dy = %.15g exactly, equal to -beta'
                  % (n, float(Fsum)))
    return Fsum


def coverage(rep, prob, ssp, dsp):
    """the cells tile the box (knot vectors span it); classification is exact; area check of the domain and patch"""
    for nm, sp in (('stress', ssp), ('displacement', dsp)):
        full = cut = 0
        area_lo = area_hi = Fr(0)
        for (i, j) in sp.cells():
            c = sp.cell(i, j)
            cls = prob.classify(c)
            if cls == 'full':
                full += 1
                area_lo += c['h'] * c['k']
                area_hi += c['h'] * c['k']
            elif cls == 'cut':
                cut += 1
                I = RG.integrate_region(X.zeros((1, 1)) + 1, c['x0'], c['h'], c['y0'], c['k'],
                                        [('circle', FD.R_HOLE2)], [])
                area_lo += RG.lo_q(I)
                area_hi += RG.hi_q(I)
        if prob.name == 'P1' and prob.hole:
            exact_area = 1 - iv.pi / 16
            ok = RG.ivq(area_lo).a <= exact_area.a and exact_area.b <= RG.ivq(area_hi).b
            rep.check(ok and area_hi - area_lo < Fr(1, 10 ** 100),
                      '%s cells: %d inside, %d cut by the hole; area of the covered domain [%.15g, %.15g] '
                      'encloses |Omega_q| = 1 - pi/16 (width %.1e)' % (nm, full, cut, float(area_lo), float(area_hi),
                                                                     float(area_hi - area_lo)))
        elif prob.name == 'P1':
            rep.check(area_lo == 1, '%s cells: %d, total area 1 exactly' % (nm, full))
        else:
            rep.check(area_lo == 3, '%s cells: %d in Omega, total area %s (exact; |Omega| = 3)' % (nm, full, area_lo))
    # patch area
    sp = ssp
    lo = hi = Fr(0)
    for (i, j) in sp.cells():
        c = sp.cell(i, j)
        if prob.classify(c) == 'out':
            continue
        p = FD.patch_integral(prob, c, X.zeros((1, 1)) + 1)
        if p is None:
            continue
        if isinstance(p, tuple):
            lo, hi = lo + p[0], hi + p[1]
        else:
            lo, hi = lo + p, hi + p
    rep.check(lo <= prob.patch_measure <= hi and hi - lo < Fr(1, 10 ** 100),
              'patch: sum over cells of |patch ∩ cell| in [%.15g, %.15g] encloses %s' % (float(lo), float(hi),
                                                                                       prob.patch_measure))


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    main()
