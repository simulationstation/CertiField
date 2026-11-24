#!/usr/bin/env python3
"""Recompute both bounds of every stored certificate in results/ from its exact data (problem, candidate, stress
field) and compare with the stored claims (exact decimals, see elasticity.claim_strings).
usage: python verify_all.py"""
import glob
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

import elasticity as E

bad = 0
here = os.path.dirname(os.path.abspath(__file__))
for js in sorted(glob.glob(os.path.join(here, 'results', '*.json'))):
    stem = js[:-5]
    t0 = time.time()
    meta = json.load(open(js))
    if meta.get('format') != 'CertiField-elasticity-certificate-1':
        continue
    ok, lb, ub = E.verify_certificate(stem, verbose=False)
    bad += (not ok)
    exact = (meta.get('info') or {}).get('exact')
    inside = ''
    if exact is not None:
        good = Fr(lb) <= Fr(exact) <= Fr(ub)
        bad += (not good)
        inside = '   exact %s %s' % (exact, 'inside' if good else 'OUTSIDE')
    P = meta['problem']
    lo, hi = float(meta['lower']), float(meta['upper'])
    print(f"{os.path.basename(stem):13s} N={'x'.join(map(str, P['N'])):6s} p={P['degree']}  claim {meta['lower']} <= inf E <= {meta['upper']}"
          f"   gap {math.ceil(1e6 * (hi - lo) / abs(hi)) / 1e4:7.4f} %   {'OK' if ok else 'FAILED'}  ({time.time() - t0:.0f} s){inside}", flush=True)
print('ALL CERTIFICATES VERIFIED' if bad == 0 else f'{bad} PROBLEM(S)')
sys.exit(1 if bad else 0)
