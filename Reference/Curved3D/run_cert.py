import numpy as np, time, sys
from fractions import Fraction as Fr
from curved3d import *
t0 = time.time()
N = int(sys.argv[1]); grip = sys.argv[2]
Ls, Lc, Lps, Lpc = [int(x) for x in sys.argv[3:7]]
P = Problem(N, grip=grip, particular='step', column=Fr(1, 2), kflat=(3 * N // 4 if grip == 'frictionless' else None),
            L_max=max(Ls, Lc, Lps, Lpc))
r = certify(P, L_solve=Ls, L_cert=Lc, Lp_solve=Lps, Lp_cert=Lpc)
print("total time", round(time.time() - t0, 1), "s")
np.savez(f"cert_{grip}_N{N}_L{Ls}{Lc}{Lps}{Lpc}.npz", psi=r['psi'], v=r['v'], K_low=r['K_low'], K_up=r['K_up'])
