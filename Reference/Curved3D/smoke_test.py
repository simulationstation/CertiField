#!/usr/bin/env python3
"""Small end-to-end test: certify the block with a spherical cavity on grid 8 for both kinds of grips (about 40 s
each) and compare with the stored grid-8 certificates.  usage: python smoke_test.py"""
import sys
from fractions import Fraction as Fr
import numpy as np
from curved3d import Problem, certify

ok = True
for grip in ("clamped", "frictionless"):
    P = Problem(8, grip=grip, particular="step", column=Fr(1, 2), kflat=(6 if grip == "frictionless" else None), L_max=5)
    r = certify(P, L_solve=2, L_cert=3, Lp_solve=3, Lp_cert=5)
    z = np.load(f"certificates/cert_{grip}_N8_L2335.npz")
    lo, up, slo, sup = float(r["K_low"]), float(r["K_up"]), float(z["K_low"]), float(z["K_up"])
    good = lo <= up and abs(lo - slo) < 1e-6 and abs(up - sup) < 1e-6
    ok &= good
    print(f"{grip:13s} {lo:.6f} <= K <= {up:.6f}   stored {slo:.6f}, {sup:.6f}   {'ok' if good else 'MISMATCH'}", flush=True)
print("SMOKE TEST PASSED" if ok else "SMOKE TEST FAILED")
sys.exit(0 if ok else 1)
