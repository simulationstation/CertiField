"""Resumable two-process run of the separate exact checker of P3, Reference/Curved3D/exact_check.py.

    python p3_exact.py <grip> <psi file> --v <v file> [--procs 2]

It performs the checks (a)-(h) of exact_check.main() with exact_check's own functions (stress_polys, remainder,
cell_class, spline_poly, coverage_check, ...), but slab by slab (the layers of cells with the same x index): two
processes take two halves of the slabs, each saves a checkpoint after every slab (cache/), and the faces between the
two halves are matched at the end.  An interrupted run resumes from the last saved slab.  Nothing in
Reference/Curved3D is modified.  The problem parameters are those of exact_check.main (step particular field, column
1/2, flat part 3N/4 for frictionless grips).
"""
import argparse
import hashlib
import itertools
import multiprocessing as mp_
import os
import pickle
import random
import sys
import time
from fractions import Fraction as Fr

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '2')

HERE = os.path.dirname(os.path.abspath(__file__))
C3DIR = os.path.normpath(os.path.join(HERE, '..', '..', 'Reference', 'Curved3D'))
sys.path.insert(0, C3DIR)

import numpy as np                       # noqa: E402

import exact_check as EC                 # noqa: E402
from exact_check import Poly             # noqa: E402

_S = {}


def problem(N, grip):
    return EC.Problem(N, grip=grip, particular='step', column=Fr(1, 2),
                      kflat=(3 * N // 4 if grip == 'frictionless' else None), L_max=3)


def init(grip, psi_path, v_path):
    psi = np.load(psi_path)['psi']
    v = np.load(v_path)['v'] if v_path else None
    N = psi.shape[1] - 2
    _S.update(grip=grip, psi=psi, v=v, N=N, P=problem(N, grip), basis=EC.make_basis())


def save(path, state):
    tmp = path + '.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(state, f)
    os.replace(tmp, path)


def run_range(args):
    """slabs a..b-1 in order; checkpoint after each slab"""
    a, b, ckpt = args
    grip, psi, v, N, P, basis = (_S[k] for k in ('grip', 'psi', 'v', 'N', 'P', 'basis'))
    one, zero = Fr(1), Fr(0)
    if os.path.exists(ckpt):
        with open(ckpt, 'rb') as f:
            st = pickle.load(f)
    else:
        st = dict(next=a, count=dict(div=0, sphere=0, free=0, mirror=0, grip=0, face=0, vface=0), res_pot=Fr(0),
                  res_tot=Fr(0), up=None, vup=None, first_low=None, vfirst_low=None)
    for i in range(st['next'], b):
        t0 = time.time()
        cnt = st['count']
        low, up, vlow, vup = {}, {}, {}, {}
        faces, vfaces = {}, {}
        for j, k in itertools.product(range(N), repeat=2):
            e = (i, j, k)
            pot, tot, rho = EC.stress_polys(P, psi, e, basis)
            for r in range(3):                                            # (a)
                d = Poly()
                for c in range(3):
                    d = d + tot[(r, c)].diff(c)
                assert d.is_zero(), ('div tau != 0', e, r)
            cnt['div'] += 1
            if EC.cell_class(P, e) == 1:                                  # (c)
                ax = next(m for m in range(3) if P.eta[m] is not None)
                for r in range(3):
                    g = Poly()
                    for c in range(3):
                        g = g + tot[(r, c)] * rho.diff(c)
                    assert EC.remainder(g, rho, ax).is_zero(), ('traction on the curved surface', e, r)
                cnt['sphere'] += 1
            for m in (1, 2):                                              # (d)
                if e[m] == N - 1:
                    for r in range(3):
                        assert tot[(r, m)].subs(m, one).is_zero(), ('traction on a free face', e, m, r)
                    cnt['free'] += 1
            for m in range(3):
                if e[m] == 0:
                    for r in range(3):
                        if r != m:
                            assert tot[(r, m)].subs(m, zero).is_zero(), ('shear on a mirror plane', e, m, r)
                    cnt['mirror'] += 1
            if i == N - 1:
                if grip == 'frictionless':
                    for r in (1, 2):
                        assert tot[(r, 0)].subs(0, one).is_zero(), ('shear on the grip', e, r)
                cnt['grip'] += 1
                lo, hi = (Fr(0),) * 3, (Fr(1),) * 3
                st['res_pot'] += pot[(0, 0)].subs(0, one).integral(lo, hi) / (N * N)
                st['res_tot'] += tot[(0, 0)].subs(0, one).integral(lo, hi) / (N * N)
            # (b) faces: y and z inside the slab; x faces kept for the neighbouring slabs
            for m in (1, 2):
                key = (e, m)
                if key in faces:
                    assert faces.pop(key) == [tot[(r, m)].subs(m, zero) for r in range(3)], ('traction jump', e, m)
                    cnt['face'] += 1
                if e[m] < N - 1:
                    nb = list(e)
                    nb[m] += 1
                    faces[(tuple(nb), m)] = [tot[(r, m)].subs(m, one) for r in range(3)]
            low[(j, k)] = [tot[(r, 0)].subs(0, zero) for r in range(3)]
            if i < N - 1:
                up[(j, k)] = [tot[(r, 0)].subs(0, one) for r in range(3)]
            # (h) displacement
            if v is not None:
                vp = [EC.spline_poly(v[c], e, basis) for c in range(3)]
                for m in (1, 2):
                    key = (e, m)
                    if key in vfaces:
                        assert vfaces.pop(key) == [vp[c].subs(m, zero) for c in range(3)], ('displacement jump', e, m)
                        cnt['vface'] += 1
                    if e[m] < N - 1:
                        nb = list(e)
                        nb[m] += 1
                        vfaces[(tuple(nb), m)] = [vp[c].subs(m, one) for c in range(3)]
                for m in range(3):
                    if e[m] == 0:
                        assert vp[m].subs(m, zero).is_zero(), ('normal displacement on a mirror plane', e, m)
                if i == N - 1:
                    assert vp[0].subs(0, one).is_zero(), ('grip: u_x', e)
                    if grip == 'clamped':
                        assert vp[1].subs(0, one).is_zero() and vp[2].subs(0, one).is_zero(), ('grip: u_y, u_z', e)
                vlow[(j, k)] = [vp[c].subs(0, zero) for c in range(3)]
                if i < N - 1:
                    vup[(j, k)] = [vp[c].subs(0, one) for c in range(3)]
        assert not faces and not vfaces, 'unmatched faces inside a slab'
        if i == a:
            st['first_low'], st['vfirst_low'] = low, vlow
        else:
            assert st['up'] == low, ('traction jump across x faces', i)
            cnt['face'] += len(low)
            if v is not None:
                assert st['vup'] == vlow, ('displacement jump across x faces', i)
                cnt['vface'] += len(vlow)
        st['up'], st['vup'] = up, vup
        st['next'] = i + 1
        save(ckpt, st)
        print('   slab %d done (%d cells, %.0f s)' % (i, N * N, time.time() - t0), flush=True)
    return st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('grip', choices=('clamped', 'frictionless'))
    ap.add_argument('psi')
    ap.add_argument('--v', default=None)
    ap.add_argument('--procs', type=int, default=2)
    a = ap.parse_args()
    t0 = time.time()
    init(a.grip, a.psi, a.v)
    N, P, psi, v, basis = _S['N'], _S['P'], _S['psi'], _S['v'], _S['basis']
    Pm, own = EC.psi_ties(N, a.grip, owners=True)
    assert EC.ties_hold(Pm, own, psi), 'ties violated'
    if v is not None:
        Pv, ownv = EC.v_ties(N, a.grip, owners=True)
        assert EC.ties_hold(Pv, ownv, v), 'displacement: ties violated'
    h = hashlib.sha256()
    for p in [a.psi] + ([a.v] if a.v else []) + [os.path.join(C3DIR, 'exact_check.py'), os.path.join(C3DIR, 'curved3d.py'),
                                                  os.path.abspath(__file__)]:
        with open(p, 'rb') as f:
            h.update(f.read())
    cdir = os.path.join(HERE, 'cache', 'p3_%s_N%d_%s' % (a.grip, N, h.hexdigest()[:16]))
    os.makedirs(cdir, exist_ok=True)
    print('grid %d, %s grips: exact check of all %d cells, %d slabs in %d ranges; checkpoints in %s'
          % (N, a.grip, N ** 3, N, a.procs, os.path.relpath(cdir, HERE)), flush=True)
    bounds = [round(r * N / a.procs) for r in range(a.procs + 1)]
    jobs = [(bounds[r], bounds[r + 1], os.path.join(cdir, 'range_%d_%d.pkl' % (bounds[r], bounds[r + 1])))
            for r in range(a.procs)]
    if a.procs > 1:
        with mp_.Pool(a.procs, initializer=init, initargs=(a.grip, a.psi, a.v)) as pool:
            states = pool.map(run_range, jobs)
    else:
        states = [run_range(j) for j in jobs]
    cnt = {k: sum(s['count'][k] for s in states) for k in states[0]['count']}
    for s1, s2 in zip(states[:-1], states[1:]):                     # faces between the ranges
        assert s1['up'] == s2['first_low'], 'traction jump between ranges'
        cnt['face'] += len(s2['first_low'])
        if v is not None:
            assert s1['vup'] == s2['vfirst_low'], 'displacement jump between ranges'
            cnt['vface'] += len(s2['vfirst_low'])
    res_pot = sum((s['res_pot'] for s in states), Fr(0))
    res_tot = sum((s['res_tot'] for s in states), Fr(0))
    nf = 3 * N * N * (N - 1)
    assert cnt['face'] == nf and (v is None or cnt['vface'] == nf), (cnt, nf)
    print('   (a) div tau = 0 exactly in %d cells' % cnt['div'])
    print('   (b) normal traction continuous across %d interior faces' % cnt['face'])
    print('   (c) zero traction on the curved surface, as a polynomial identity, in %d cut cells' % cnt['sphere'])
    print('   (d) zero traction on %d free-face cells; zero shear on %d mirror-plane cells%s'
          % (cnt['free'], cnt['mirror'], ' and on the grip' if a.grip == 'frictionless' else ''))
    assert res_pot == 0 and res_tot == P.area[0] / 4, (res_pot, res_tot)
    print('   (e) force through the grip (octant): potential part %s, total %s (expected 0 and %s)'
          % (res_pot, res_tot, P.area[0] / 4))
    # (f) as in exact_check.main: exact integrals on six boxes against curved3d's interval bounds
    bx_all = P.boxes(3)
    rnd = random.Random(2)
    pick = rnd.sample([b for b in bx_all if b[7] == 1], 3) + rnd.sample([b for b in bx_all if b[3] == 0], 3)
    worst = 0.0
    for b in pick:
        e = tuple(int(x) for x in b[:3])
        s = Fr(1, 2 ** int(b[3]))
        lo = tuple(int(b[4 + k]) * s for k in range(3))
        hi = tuple(x + s for x in lo)
        _, tot, _ = EC.stress_polys(P, psi, e, basis)
        ex = EC.q_dual(P, tot).integral(lo, hi) / N ** 3
        b1 = np.array([b])
        up = EC.c3.dual_compliance(P, b1, psi, rigorous=True)
        fl = EC.c3.dual_compliance(P, b1, psi, rigorous=False)
        assert ex >= 0 and Fr(up) >= ex, ('interval bound below the exact integral', b)
        rel_up, rel_fl = float((Fr(up) - ex) / ex), float(abs(Fr(fl) - ex) / ex)
        worst = max(worst, rel_up, rel_fl)
        assert rel_up < 1e-8 and rel_fl < 1e-8
        if v is not None and v.shape[1] == psi.shape[1]:
            exv = EC.q_primal(P, v, e, basis).integral(lo, hi) / N ** 3
            upv = EC.c3.primal_energy(P, b1, v, rigorous=True)
            assert Fr(upv) >= exv and float((Fr(upv) - exv) / exv) < 1e-9, ('strain energy', b)
    print('   (f) exact integrals on %d boxes: the interval bounds are above them, by at most %.1e (relative)'
          % (len(pick), worst))
    for L in (3, 4, 5, 6):
        print('   (g) %s   %.0fs' % (EC.coverage_check(P, L), time.time() - t0), flush=True)
    if v is not None:
        print('   (h) displacement: continuous across %d faces, grip and mirror conditions exact' % cnt['vface'])
    print('EXACT CHECK PASSED   %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    main()
