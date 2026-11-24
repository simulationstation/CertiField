#!/usr/bin/env python3
"""The numbers of version 0.3.0 must not move.  For every test vector in testvectors/ (written by the reference of
version 0.3.0) this recomputes the competitor energy, the boundary flux and the lower bound from the stored inputs
and compares them with the stored outputs bit for bit.
usage: python check_unchanged.py        (seconds)"""
import json
import os
import sys
from fractions import Fraction as Fr

import numpy as np

import scalarnd as S

here = os.path.dirname(os.path.abspath(__file__))
tv = os.path.join(here, 'testvectors')
index = json.load(open(os.path.join(tv, 'index.json')))
bad = 0
for name, entry in index.items():
    d = os.path.join(tv, name)
    if not os.path.exists(os.path.join(d, 'field.json')):
        continue
    meta = json.load(open(os.path.join(d, 'field.json')))
    shape = tuple(meta['NodeDims'])
    gval = np.fromfile(os.path.join(d, 'gval.f64'), dtype='<f8').reshape(shape)
    U = np.fromfile(os.path.join(d, 'U.f64'), dtype='<f8').reshape(shape)
    P = S.Problem.from_spec(meta['problem'], gval)
    K = meta['K']
    cal = S.Calibration(P, Fr(meta['TRange'][0]), Fr(meta['TRange'][1]), K, affine=meta['Affine'])
    for i in range(P.n):
        fs = tuple(meta['FaceDims'][i])
        cal.S0[i] = np.fromfile(os.path.join(d, f'S0_{i}.f64'), dtype='<f8').reshape(fs)
        if K > 0:
            cal.C[i] = np.fromfile(os.path.join(d, f'C{i}.f64'), dtype='<f8').reshape(fs + (K + 2,))
    ex = entry['expected']
    ub = S.upper_bound(P, U)
    fl = S.flux_lower(P, cal)
    lb = S.lower_bound(P, cal, Qs=ex['Subintervals'])
    same = (list(ub) == ex['upper_enclosure'] and fl == ex['flux_lower'] and lb == ex['lower'])
    bad += not same
    print(f"{name:10s} upper enclosure {'identical' if list(ub) == ex['upper_enclosure'] else 'DIFFERENT'}   "
          f"flux {'identical' if fl == ex['flux_lower'] else 'DIFFERENT'}   "
          f"lower bound {'identical' if lb == ex['lower'] else 'DIFFERENT (%r, stored %r)' % (lb, ex['lower'])}")
print('ALL NUMBERS OF VERSION 0.3.0 REPRODUCED BIT FOR BIT' if bad == 0 else f'{bad} TEST VECTOR(S) DIFFER')
sys.exit(1 if bad else 0)
