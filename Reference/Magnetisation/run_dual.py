#!/usr/bin/env python3
"""Search for a certificate of Proposition 5 for the Zeeman model (A = 1/2, kappa = 1, K = 0, degree -1), check it in
exact arithmetic, test the pointwise inequality independently, and store it.
usage:  python run_dual.py H [H ...]       e.g.  python run_dual.py 2 5 10
Writes results/zeeman_h<H>.json (the certificate and the exact bound)."""
import json
import os
import sys
import time
from fractions import Fraction as Fr

import magnet as M
import dual_opt as O
from check_dual import pointwise_ok

sys.stdout.reconfigure(encoding="utf-8")
os.makedirs('results', exist_ok=True)
for arg in sys.argv[1:]:
    h = Fr(arg)
    model = M.Model(A=Fr(1, 2), kappa=1, h=h)
    t0 = time.time()
    cert, info = O.build_certificate(model)
    t1 = time.time()
    ok = cert.check()
    t2 = time.time()
    if not ok:
        print(f"h = {h}: no valid certificate (cells without a valid u: {info['failed_cells']})", flush=True)
        continue
    okp, worst = pointwise_ok(cert, 300)
    c, d = cert.bound(-1)
    low = cert.lower(-1)
    out = cert.to_json()
    out['degree'] = -1
    out['bound_coefficient_of_pi'] = str(c)
    out['lower'] = repr(M.floor_float(low))
    out['info'] = dict(search_seconds=round(t1 - t0, 1), check_seconds=round(t2 - t1, 1), damping=info['damping'])
    name = f"results/zeeman_h{str(h).replace('/', 'o')}.json"
    json.dump(out, open(name, 'w'))
    print(f"h = {h}: VALID (exact check {t2 - t1:.1f} s, search {t1 - t0:.0f} s, damping {info['damping']:g}); "
          f"E >= {M.floor_float(low):.9f} = 4 pi * {float(c) / 4:.6f};  simple bound 4 pi * {float(1 - 1 / h):.6f};  "
          f"pointwise test {'ok' if okp else 'FAILED'}", flush=True)
