#!/usr/bin/env python3
"""Run the three-dimensional reference examples and store a certificate (problem, candidate, stress field, claimed
bounds) for each.
usage:  python run_examples3d.py [NAME ...]      without names: all of them (a few minutes)
        python run_examples3d.py --table [N ...] the block of E-T4 on N^3 cells (default 2 3 4 5 6), nothing stored
Run it from this directory.  It overwrites results3d/NAME.json and results3d/NAME.npz."""
import os
import sys
import time
from fractions import Fraction as Fr

import elasticity3d as E

sys.stdout.reconfigure(encoding="utf-8")
lam, mu = Fr(3, 2), Fr(1)                    # Lame constants of the block examples (Poisson ratio 0.3)
DD, DN, NN = ('D', 'D'), ('D', 'N'), ('N', 'N')
RING = (DD, NN, NN)
S1 = E.Problem3D(box=[(0, 1)] * 3, N=(1, 1, 1), lam=Fr(2), mu=Fr(1)).S          # dyadic compliance for the exact example


def block(N, p=2, load='pull', **kw):
    """E-T4: unit cube, clamped on the face x1 = 0, displaced on the face x1 = 1, no traction on the four others"""
    g = {'pull': lambda X, Y, Z: (0.1 * X, 0 * Y, 0 * Z),
         'shear': lambda X, Y, Z: (0 * X, 0.1 * X, 0 * Z),
         'twist': lambda X, Y, Z: (0 * X, -0.1 * X * (Z - 0.5), 0.1 * X * (Y - 0.5))}[load]
    return E.Problem3D(box=[(0, 1)] * 3, N=(N, N, N), lam=lam, mu=mu, degree=p, bc=RING, g=g, name=f'block3d_{load}_p{p}_{N}', **kw)


EXAMPLES = {
    'block3d_4': (lambda: block(4), None, 'cube clamped on one face, pulled on the opposite face, four free faces'),
    'block3d_6': (lambda: block(6), None, 'the same'),
    'shear3d_4': (lambda: block(4, load='shear'), None, 'the same cube, end face displaced sideways'),
    'twist3d_4': (lambda: block(4, load='twist'), None, 'the same cube, end face rotated about the axis'),
    'bending3d': (lambda: E.Problem3D(box=[(0, 2), (Fr(-1, 2), Fr(1, 2)), (Fr(-1, 2), Fr(1, 2))], N=(2, 2, 2), lam=Fr(2), mu=Fr(1), degree=2, bc=RING,
                                      g=lambda X, Y, Z: (float(S1[0][0]) * X * Y, -float(S1[0][0]) * X * X / 2 + float(S1[1][0]) * (Y * Y - Z * Z) / 2,
                                                         float(S1[1][0]) * Y * Z), name='bending3d'),
                  S1[0][0] / 12, 'pure bending of a bar with four free faces (exact minimum 1/32)'),
    'gravity3d': (lambda: E.Problem3D(box=[(0, 1), (0, 1), (0, 2)], N=(3, 3, 6), lam=lam, mu=mu, degree=2, bc=(NN, NN, DN), b=(Fr(0), Fr(0), Fr(-1)),
                                      g=lambda X, Y, Z: (0 * X, 0 * Y, 0 * Z), name='gravity3d'),
                  None, 'column clamped at the bottom under its own weight, five free faces'),
    'plate3d': (lambda: E.Problem3D(box=[(0, 1), (0, Fr(1, 2)), (0, 1)], N=(4, 2, 4), lam=lam, mu=mu, degree=2, bc=(DD, NN, DD), b=(Fr(0), Fr(-1), Fr(0)),
                                    g=lambda X, Y, Z: (0 * X, 0 * Y, 0 * Z), name='plate3d'),
                None, 'plate clamped on its four edges under its own weight, two free faces'),
}


def run(name):
    mk, exact, note = EXAMPLES[name]
    P = mk()
    out, U, st = E.certify(P)
    info = dict(example=name, note=note, seconds=out['seconds'], exact=(str(exact) if exact is not None else None), solver=st.info)
    E.save_certificate(f'results3d/{name}', P, U, st, out['lower'], out['upper'], extra=info)
    ok, lb, ub = E.verify_certificate(f'results3d/{name}', verbose=False)
    inside = '' if exact is None else ('   exact %.10f %s' % (float(exact), 'inside' if Fr(lb) <= exact <= Fr(ub) else 'OUTSIDE'))
    print(f"{name:11s} N={'x'.join(map(str, P.N)):6s} p={P.p}  {lb:+.10f} <= inf E <= {ub:+.10f}   gap {100 * (ub - lb) / abs(ub):.4f} %   "
          f"{out['seconds']:.1f} s   stored-certificate check: {'OK' if ok else 'FAILED'}{inside}", flush=True)


def table(Ns, load='pull'):
    print(f"block of E-T4 ({load}), six stress functions, p = 2")
    for N in Ns:
        P = block(N, load=load)
        t0 = time.time()
        out, U, st = E.certify(P)
        print(f"  {N}x{N}x{N}:  {out['lower']:.8f} <= inf E <= {out['upper']:.8f}   gap {100 * out['rel_gap']:.3f} %   {time.time() - t0:.1f} s   "
              f"(candidate {out['seconds_candidate']:.1f} s, stress field {out['seconds_stress']:.1f} s with {st.info}, bounds {out['seconds_bounds']:.1f} s)",
              flush=True)


if __name__ == '__main__':
    args = sys.argv[1:]
    if args and args[0] == '--table':
        load = 'pull'
        rest = args[1:]
        if rest and rest[0] in ('pull', 'shear', 'twist'):
            load, rest = rest[0], rest[1:]
        table([int(x) for x in rest] or [2, 3, 4, 5, 6], load)
    else:
        os.makedirs('results3d', exist_ok=True)
        for nm in (args or list(EXAMPLES)):
            run(nm)
