#!/usr/bin/env python3
"""Run the reference examples and store a certificate (problem, candidate, field, claimed bounds) for each.
usage:  python run_examples.py NAME [NAME ...]        names: see EXAMPLES below;  'fast' = everything under a minute
Run it from this directory.  It overwrites results/NAME.json and results/NAME.npz."""
import json, sys, time
from fractions import Fraction as Fr
import numpy as np
import scalarnd as S
from ref1d import wall_energy
from refseries import screened_min

W4 = [Fr(1, 4), 0, Fr(-1, 2), 0, Fr(1, 4)]                       # (u^2 - 1)^2 / 4


def dw(lam):                                                     # lam (u^2 - 1)^2 / 4
    return [Fr(lam, 4), 0, Fr(-lam, 2), 0, Fr(lam, 4)]


def finish(name, P, U, cal, rows, ref=None, note=''):
    r = rows[-1]
    info = dict(example=name, note=note, levels=rows, reference_float=ref)
    S.save_certificate(f'results/{name}', P, U, cal, r['lower'], r['upper'], extra=info)
    ok, lb, ub = S.verify_certificate(f'results/{name}', verbose=False)
    line = (f"{name:14s} n={P.n} N={tuple(P.N)} K={cal.K:2d}  {lb:+.8f} <= inf E <= {ub:+.8f}   gap {100 * (ub - lb) / abs(ub):.3f} %"
            f"   {rows[-1]['seconds']:.0f} s   stored-certificate check: {'OK' if ok else 'FAILED'}"
            + (f"   reference {ref:+.8f} {'inside' if lb <= ref <= ub else 'OUTSIDE'}" if ref is not None else ''))
    print(line, flush=True)
    open('results/summary.txt', 'a').write(line + '\n')


def convex(name, n, N, W, g, ref=None, box=None, note=''):
    P = S.Problem(box=box or [(0, 1)] * n, N=(N,) * n, k=Fr(1, 2), W=W, bc=[('D', 'D')] * n, g=g, name=name)
    out, U, cal = S.certify_convex(P)
    finish(name, P, U, cal, [out], ref, note)


def calib(name, mk, levels, trange, U0, ref=None, scale=1.0, note='', **kw):
    rows, U, cal, P = S.certify_calibration(mk, levels, trange, U0=U0, scale=scale, log=f'results/{name}_levels.txt', **kw)
    finish(name, P, U, cal, rows, ref, note)


EXAMPLES = {}
zero = lambda *X: 0 * X[0]
# ---- convex: W = u^2 - u (screened Poisson with unit source), exact minimum from the sine series ----
EXAMPLES['screened2d'] = lambda: convex('screened2d', 2, 64, [0, -1, 1], zero, screened_min(2, 0.5, 1.0, 1.0), note='W = u^2 - u, u = 0 on the boundary of the unit square')
EXAMPLES['screened2d_128'] = lambda: convex('screened2d_128', 2, 128, [0, -1, 1], zero, screened_min(2, 0.5, 1.0, 1.0))
EXAMPLES['screened3d'] = lambda: convex('screened3d', 3, 32, [0, -1, 1], zero, screened_min(3, 0.5, 1.0, 1.0), note='W = u^2 - u, u = 0 on the boundary of the unit cube')
EXAMPLES['screened3d_48'] = lambda: convex('screened3d_48', 3, 48, [0, -1, 1], zero, screened_min(3, 0.5, 1.0, 1.0))
# ---- convex, not quadratic:  -Laplace u + u^3 + u = 1  on the unit cube (no closed form) ----
EXAMPLES['quartic3d'] = lambda: convex('quartic3d', 3, 32, [0, -1, Fr(1, 2), 0, Fr(1, 4)], zero, note='W = u^4/4 + u^2/2 - u (convex, not quadratic), u = 0 on the boundary of the unit cube')
EXAMPLES['quartic2d'] = lambda: convex('quartic2d', 2, 64, [0, -1, Fr(1, 2), 0, Fr(1, 4)], zero, note='W = u^4/4 + u^2/2 - u, unit square')
# ---- non-convex with known answer: a domain wall; exact minimum = cross-section x 1D value ----
wallg = lambda *X: np.where(X[0] < 0, -1.0, 1.0)
wallU = lambda *X: np.tanh(X[0] / np.sqrt(2))
EXAMPLES['wall2d'] = lambda: calib('wall2d', lambda N: S.Problem(box=[(-3, 3), (0, 1)], N=N, k=Fr(1, 2), W=W4, bc=[('D', 'D'), ('N', 'N')], g=wallg, name='wall2d'),
                                   [((24, 1), 8), ((48, 2), 8), ((48, 2), 16), ((48, 2), 32), ((96, 4), 32)], (-1, 1), wallU, wall_energy(3)[0],
                                   note='double well (u^2-1)^2/4, u = -1 | +1 on the ends of (-3,3) x (0,1), free on the sides')
EXAMPLES['wall3d'] = lambda: calib('wall3d', lambda N: S.Problem(box=[(-3, 3), (0, 1), (0, 1)], N=N, k=Fr(1, 2), W=W4, bc=[('D', 'D'), ('N', 'N'), ('N', 'N')], g=wallg, name='wall3d'),
                                   [((24, 1, 1), 8), ((48, 2, 2), 8), ((48, 2, 2), 16), ((48, 2, 2), 32)], (-1, 1), wallU, wall_energy(3)[0],
                                   note='the same wall in the box (-3,3) x (0,1)^2')
# ---- non-convex, genuinely multi-dimensional, no closed form ----
EXAMPLES['ramp2d'] = lambda: calib('ramp2d', lambda N: S.Problem(box=[(-1, 1)] * 2, N=N, k=Fr(1, 2), W=dw(16), bc=[('D', 'D')] * 2, g=lambda X, Y: X + 0 * Y, name='ramp2d'),
                                   [((8, 8), 8), ((16, 16), 8), ((16, 16), 16), ((32, 32), 16), ((32, 32), 32), ((64, 64), 32)], (-1, 1), lambda X, Y: X + 0 * Y, scale=4,
                                   note='double well 4 (u^2-1)^2, u = x on the boundary of (-1,1)^2')
EXAMPLES['saddle2d'] = lambda: calib('saddle2d', lambda N: S.Problem(box=[(-1, 1)] * 2, N=N, k=Fr(1, 2), W=dw(16), bc=[('D', 'D')] * 2, g=lambda X, Y: X * Y, name='saddle2d'),
                                     [((8, 8), 8), ((16, 16), 8), ((16, 16), 16), ((32, 32), 16), ((32, 32), 32), ((64, 64), 32)], (-1, 1),
                                     lambda X, Y: X * Y + 0.3 * (1 - X * X) * (1 - Y * Y), scale=4,
                                     note='double well 4 (u^2-1)^2, u = x y on the boundary of (-1,1)^2 (two symmetric minimisers)')
EXAMPLES['ramp3d'] = lambda: calib('ramp3d', lambda N: S.Problem(box=[(-1, 1)] * 3, N=N, k=Fr(1, 2), W=dw(4), bc=[('D', 'D')] * 3, g=lambda X, Y, Z: X + 0 * Y, name='ramp3d'),
                                   [((8, 8, 8), 8), ((16, 16, 16), 8), ((16, 16, 16), 16)], (-1, 1), lambda X, Y, Z: X + 0 * Y,
                                   note='double well (u^2-1)^2, u = x on the boundary of the cube (-1,1)^3')
EXAMPLES['ramp3d_32'] = lambda: calib('ramp3d_32', lambda N: S.Problem(box=[(-1, 1)] * 3, N=N, k=Fr(1, 2), W=dw(4), bc=[('D', 'D')] * 3, g=lambda X, Y, Z: X + 0 * Y, name='ramp3d_32'),
                                      [((8, 8, 8), 8), ((16, 16, 16), 8), ((16, 16, 16), 16), ((32, 32, 32), 16)], (-1, 1), lambda X, Y, Z: X + 0 * Y,
                                      last=((1e-4, 200), (3e-5, 400)),
                                      note='double well (u^2-1)^2, u = x on the boundary of the cube (-1,1)^3')
FAST = ['screened2d', 'screened2d_128', 'screened3d', 'screened3d_48', 'quartic2d', 'quartic3d']

if __name__ == '__main__':
    names = sys.argv[1:] or FAST
    if names == ['fast']:
        names = FAST
    for nm in names:
        EXAMPLES[nm]()
