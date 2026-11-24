"""upper bound only"""
import numpy as np, time, sys
from fractions import Fraction as Fr
from curved3d import *
t0 = time.time()
N = int(sys.argv[1]); grip = sys.argv[2]; Ls = int(sys.argv[3]); Lc = int(sys.argv[4])
P = Problem(N, grip=grip, particular='step', column=Fr(1, 2), kflat=(3 * N // 4 if grip == 'frictionless' else None),
            L_max=max(Ls, Lc))
bs = P.boxes(Ls)
print(f"N={N} {grip}: {len(bs)} boxes to solve (L={Ls})   {time.time() - t0:.0f}s", flush=True)
K, f, k0 = assemble_primal(P, bs)
Pv = v_ties(N, grip)
y = solve_reduced(K, -f, Pv, verbose=True, maxiter=3000)
cv = (Pv @ y).reshape(3, P.M, P.M, P.M)
np.savez(f"v_{grip}_N{N}.npz", v=cv)
bc = P.boxes(Lc)
v_all, v_in = box_volume_bounds(bc, N)
print(f"{len(bc)} boxes to certify (L={Lc}); octant volume between {float(v_in):.6f} and {float(v_all):.6f}   {time.time() - t0:.0f}s", flush=True)
Kf = 8 * primal_energy(P, bc, cv, rigorous=False)
Kup = float(np.nextafter(8 * primal_energy(P, bc, cv, rigorous=True), np.inf))
print(f"certified upper bound  K <= {Kup:.8f}   (float value {Kf:.8f})   total {time.time() - t0:.0f}s", flush=True)
