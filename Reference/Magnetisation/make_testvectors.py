#!/usr/bin/env python3
"""Test vector for the Wolfram Language model: a small perturbed field, its certified energy and degree enclosure.
Writes testvectors/small/{x1,x2,Mn}.f64 (little-endian doubles, C order) and testvectors/small/expected.json."""
import json
import os
from fractions import Fraction as Fr

import numpy as np

import magnet as M
import competitor as C

rng = np.random.default_rng(11)
d = 'testvectors/small'
os.makedirs(d, exist_ok=True)
x = C.graded_axis(2.0, 8, 1.5)
Mn = C.radial_field(x, x, C.belavin_polyakov(0.6, 2.0), 2.0)
Mn[1:-1, 1:-1] += 0.02 * rng.standard_normal(Mn[1:-1, 1:-1].shape)
out = {}
for name, mod in [('zeeman', M.Model(A=Fr(1, 2), kappa=1, h=2)), ('all', M.Model(A=Fr(3, 4), kappa=Fr(1, 2), K=Fr(-1, 4), h=Fr(3, 2))),
                  ('positiveK', M.Model(A=Fr(1, 2), kappa=Fr(1, 4), K=Fr(1, 2), h=0))]:
    ub, parts = C.upper_bound(mod, x, x, Mn)
    out[name] = dict(model=mod.spec(), upper=ub, parts=parts)
Q, enc = C.degree(x, x, Mn)
out['degree'] = Q
out['degree_enclosure'] = list(enc)
for nm, arr in (('x1', x), ('x2', x), ('Mn', Mn)):
    np.asarray(arr, dtype='<f8').tofile(os.path.join(d, nm + '.f64'))
json.dump(dict(dims=list(Mn.shape), **out), open(os.path.join(d, 'expected.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
