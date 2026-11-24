#!/usr/bin/env python3
"""Run the reference examples and store a certificate (problem, candidate, stress field, claimed bounds) for each.
usage:  python run_examples.py [NAME ...]        without names: all of them (a few seconds)
Run it from this directory.  It overwrites results/NAME.json and results/NAME.npz."""
import os
import sys
from fractions import Fraction as Fr

import elasticity as E

lam, mu = Fr(3, 2), Fr(1)                    # Lame constants of the block examples (Poisson ratio 0.3)
S1 = E.Problem(box=[(0, 1), (0, 1)], N=(1, 1), lam=Fr(1), mu=Fr(1)).S


def block(N, p):
    """E-T4: unit square, clamped on the left, displaced by (0.1, 0) on the right, no traction on top and bottom"""
    return E.Problem(box=[(0, 1), (0, 1)], N=(N, N), lam=lam, mu=mu, degree=p, bc=(('D', 'D'), ('N', 'N')),
                     g=lambda X, Y: (0.1 * X, 0 * Y), name=f'block_p{p}_{N}')


EXAMPLES = {
    'block_p1_16': (lambda: block(16, 1), None, 'block clamped on the left, pulled on the right, free sides'),
    'block_p1_32': (lambda: block(32, 1), None, 'the same'),
    'block_p2_16': (lambda: block(16, 2), None, 'the same'),
    'block_p2_32': (lambda: block(32, 2), None, 'the same'),
    'bending': (lambda: E.Problem(box=[(0, 2), (Fr(-1, 2), Fr(1, 2))], N=(4, 2), lam=Fr(1), mu=Fr(1), degree=2, bc=(('D', 'D'), ('N', 'N')),
                                  g=lambda X, Y: (float(S1[0][0]) * X * Y, float(S1[1][0]) * Y * Y / 2 - float(S1[0][0]) * X * X / 2), name='bending'),
                S1[0][0] / 12, 'pure bending with free sides (exact minimum 1/32)'),
    'gravity': (lambda: E.Problem(box=[(0, 1), (0, 2)], N=(16, 32), lam=lam, mu=mu, degree=2, bc=(('N', 'N'), ('D', 'N')), b=(Fr(0), Fr(-1)),
                                  g=lambda X, Y: (0 * X, 0 * Y), name='gravity'),
                None, 'column clamped at the bottom under its own weight, three free faces'),
}


def run(name):
    mk, exact, note = EXAMPLES[name]
    P = mk()
    out, U, st = E.certify(P)
    info = dict(example=name, note=note, seconds=out['seconds'], exact=(str(exact) if exact is not None else None))
    E.save_certificate(f'results/{name}', P, U, st, out['lower'], out['upper'], extra=info)
    ok, lb, ub = E.verify_certificate(f'results/{name}', verbose=False)
    inside = '' if exact is None else ('   exact %.10f %s' % (float(exact), 'inside' if Fr(lb) <= exact <= Fr(ub) else 'OUTSIDE'))
    print(f"{name:13s} N={P.N} p={P.p}  {lb:+.10f} <= inf E <= {ub:+.10f}   gap {100 * (ub - lb) / abs(ub):.4f} %   "
          f"{out['seconds']:.1f} s   stored-certificate check: {'OK' if ok else 'FAILED'}{inside}", flush=True)


if __name__ == '__main__':
    os.makedirs('results', exist_ok=True)
    for nm in (sys.argv[1:] or list(EXAMPLES)):
        run(nm)
