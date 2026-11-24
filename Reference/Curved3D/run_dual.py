"""lower bound only, on a finer grid (the upper bound of a coarser grid stays valid for the same body)"""
import numpy as np, time, sys
from fractions import Fraction as Fr
from curved3d import *
t0 = time.time()
N = int(sys.argv[1]); grip = sys.argv[2]; Ls = int(sys.argv[3]); Lc = int(sys.argv[4])
P = Problem(N, grip=grip, particular='step', column=Fr(1, 2), kflat=(3 * N // 4 if grip == 'frictionless' else None),
            L_max=max(Ls, Lc))
bs = P.boxes(Ls)
print(f"N={N} {grip}: {len(bs)} boxes to solve (L={Ls})", flush=True)
A, b, a0 = assemble_dual(P, bs, verbose=False)
print("assembled", A.shape, A.nnz, round(time.time() - t0), "s", flush=True)
Pp = psi_ties(N, grip)
y = solve_reduced(A, -b, Pp, verbose=True, maxiter=3000)
del A
cp = (Pp @ y).reshape(6, P.M, P.M, P.M)
np.savez(f"psi_{grip}_N{N}.npz", psi=cp)
bc = P.boxes(Lc)
v_all, v_in = box_volume_bounds(bc, N)
print(f"{len(bc)} boxes to certify (L={Lc}); octant volume between {float(v_in):.6f} and {float(v_all):.6f}", flush=True)
Qf = 8 * dual_compliance(P, bc, cp, rigorous=False)
print(f"float value K_low = {16 / Qf:.8f}   {time.time() - t0:.0f}s", flush=True)
Q_up = float(np.nextafter(8 * dual_compliance(P, bc, cp, rigorous=True, verbose=False), np.inf))
alo = frac_enclosure(P.area[0])[0]
Klow = float(dn(dn(alo * alo) / Q_up))
print(f"certified lower bound  K >= {Klow:.8f}    total {time.time() - t0:.0f}s", flush=True)
