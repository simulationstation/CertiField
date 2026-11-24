#!/usr/bin/env python3
"""M1: certified intervals for the infimum of the Zeeman model (A = 1/2, kappa = 1, K = 0) in degree -1 on a square
[-L, L]^2 with m = e3 on the boundary.  Lower bound: the stored certificate results/zeeman_h<H>.json (valid on every
domain), re-checked exactly.  Upper bound: the certified energy of the normalised interpolant of the axisymmetric
profile on a graded grid, with certified degree.
usage:  python run_m1.py H L N BETA [H L N BETA ...]        e.g.  python run_m1.py 2 4 160 4 10 2 200 6
Writes results/m1_h<H>.json (the interval)."""
import json
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import magnet as M
import competitor as C

sys.stdout.reconfigure(encoding="utf-8")
args = sys.argv[1:]
for k in range(0, len(args), 4):
    h, L, n, beta = Fr(args[k]), Fr(args[k + 1]), int(args[k + 2]), float(args[k + 3])
    tag = str(h).replace('/', 'o')
    t0 = time.time()
    d = json.load(open(f'results/zeeman_h{tag}.json'))
    cert = M.Certificate.from_json(d)
    assert cert.model.h == h and cert.model.kappa == 1 and cert.model.A == Fr(1, 2) and cert.model.K == 0
    ok = cert.check()
    lb = M.floor_float(cert.lower(-1, area=4 * L * L))
    t1 = time.time()
    prof = C.skyrmion_profile(cert.model, float(L))
    x = C.graded_axis(float(L), n, beta)
    assert x[0] == -float(L) and x[-1] == float(L)
    Mn = C.radial_field(x, x, prof, float(L))
    ub, parts = C.upper_bound(cert.model, x, x, Mn)
    Q, enc = C.degree(x, x, Mn)
    t2 = time.time()
    good = ok and Q == -1 and lb <= ub
    gap = (ub - lb) / abs(ub)
    print(f"h = {h}, square [-{L}, {L}]^2, grid {2 * n}x{2 * n}: {lb:.8f} <= inf E <= {ub:.8f}   gap {100 * gap:.3f} %   degree {Q}   "
          f"{'OK' if good else 'FAILED'}   (lower {t1 - t0:.1f} s, upper {t2 - t1:.1f} s)", flush=True)
    json.dump(dict(problem=dict(A='1/2', kappa='1', K='0', h=str(h), degree=-1, domain=[[str(-L), str(L)], [str(-L), str(L)]]),
                   lower=repr(lb), upper=repr(ub), gap=gap, grid=[2 * n, 2 * n], parts=parts, degree_enclosure=list(enc)),
              open(f'results/m1_h{tag}.json', 'w'), indent=1)
