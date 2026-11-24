#!/usr/bin/env python3
"""verify.py -- independent re-check of a stored certificate.

A certificate is a pair of coefficient arrays (psi: stress potential, v: displacement) and the description of the body.
Nothing from the solver is trusted: the script checks that the arrays satisfy the admissibility ties exactly and
re-evaluates both bounds with outward-rounded interval arithmetic.

usage:  python3 verify.py <N> <grip> <psi.npz> <v.npz> [L_psi] [L_v]
"""
import sys, time
from fractions import Fraction as Fr
import numpy as np
from curved3d import *

def main():
    N, grip = int(sys.argv[1]), sys.argv[2]
    psi = np.load(sys.argv[3])['psi']
    v = np.load(sys.argv[4])['v']
    Lp = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    Lv = int(sys.argv[6]) if len(sys.argv) > 6 else 6
    t0 = time.time()
    Np, Nv = psi.shape[1] - 2, v.shape[1] - 2
    out = {}
    # ---- lower bound
    P = Problem(Np, grip=grip, particular='step', column=Fr(1, 2),
                kflat=(3 * Np // 4 if grip == 'frictionless' else None), L_max=Lp)
    Pm, own = psi_ties(Np, grip, owners=True)
    assert psi.dtype == np.float64 and np.all(np.isfinite(psi)) and ties_hold(Pm, own, psi), "stress potential: ties violated"
    bc = P.boxes(Lp)
    Q_up = float(np.nextafter(8 * dual_compliance(P, bc, psi, rigorous=True), np.inf))
    alo = frac_enclosure(P.area[0])[0]
    K_low = float(dn(dn(alo * alo) / Q_up))
    print(f"lower bound (grid {Np}, depth {Lp}, {len(bc)} boxes):  K >= {K_low:.8f}    {time.time() - t0:.0f}s", flush=True)
    # ---- upper bound
    P = Problem(Nv, grip=grip, particular='step', column=Fr(1, 2),
                kflat=(3 * Nv // 4 if grip == 'frictionless' else None), L_max=Lv)
    Pm, own = v_ties(Nv, grip, owners=True)
    assert v.dtype == np.float64 and np.all(np.isfinite(v)) and ties_hold(Pm, own, v), "displacement: ties violated"
    bc = P.boxes(Lv)
    K_up = float(np.nextafter(8 * primal_energy(P, bc, v, rigorous=True), np.inf))
    print(f"upper bound (grid {Nv}, depth {Lv}, {len(bc)} boxes):  K <= {K_up:.8f}    {time.time() - t0:.0f}s", flush=True)
    assert K_low <= K_up
    print(f"VERIFIED  {K_low:.6f} <= K <= {K_up:.6f}   gap {100 * (K_up - K_low) / K_low:.3f} %")

if __name__ == '__main__':
    main()
