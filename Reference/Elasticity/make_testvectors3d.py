#!/usr/bin/env python3
"""Small test vectors of the three-dimensional reference for a port: inputs as raw doubles, expected outputs.
Each case stores gval.f64 (3 x nodes), U.f64 (3 x nodes), X.f64 (functions x 8 x grid nodes), par.f64, field.json and
states.json (the polynomials Psi_b, Psi_1 ... Psi_r as lists of [column, "numerator/denominator"], column =
64 f + 16 g1 + 4 g2 + g3 for the coefficient of x1^g1 x2^g2 x3^g3 in the f-th function in use);
index.json holds a suggested model expression for the Wolfram Language and the expected numbers."""
import json
import os
import shutil
import sys
from fractions import Fraction as Fr

import numpy as np

import elasticity3d as E

sys.stdout.reconfigure(encoding="utf-8")
rng = np.random.default_rng(17)
OUT = 'testvectors3d'
shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(OUT)
index = {}
DD, DN, ND, NN = ('D', 'D'), ('D', 'N'), ('N', 'D'), ('N', 'N')
lam, mu = Fr(3, 2), Fr(1)
Daniso = [[(lam + 2 * mu if i == j else lam) if (i < 3 and j < 3) else (mu if i == j else Fr(0)) for j in range(6)] for i in range(6)]
for (i, j), v in {(0, 3): Fr(1, 4), (0, 5): Fr(-1, 8), (1, 4): Fr(1, 5), (2, 5): Fr(1, 3), (3, 4): Fr(1, 8), (4, 5): Fr(-1, 6), (0, 1): Fr(1, 4)}.items():
    Daniso[i][j] += v
    Daniso[j][i] += v
gq = lambda X, Y, Z: (0.1 * X - 0.15 * Y + 0.3 * X * Z, 0.25 * X + 0.2 * Y * Z - 0.1, 0.05 * Z - 0.2 * X * Y + 0.1 * Y)
GQ = 'Function[{x, y, z}, {0.1 x - 0.15 y + 0.3 x z, 0.25 x + 0.2 y z - 0.1, 0.05 z - 0.2 x y + 0.1 y}]'
DWL = '{' + ', '.join('{' + ', '.join(str(x) for x in row) + '}' for row in Daniso) + '}'


def record(name, P, wl_model, randomise=False):
    out, U, st = E.certify(P)
    if randomise:                                # any nodal values and parameters are admissible
        U = U + 0.2 * rng.standard_normal(U.shape)
        U[:, P.fixed] = P.gval[:, P.fixed]
        st = E.Stress(P, st.X + 0.2 * rng.standard_normal(st.X.shape), st.par + 0.2 * rng.standard_normal(st.par.shape))
    lb, parts = E.lower_bound(P, st, return_parts=True)
    ub = E.upper_bound(P, U)
    mask, Tb, T = P.prescribed()
    psi_b, psis = P.states()
    exp = dict(lower=lb, boundary=list(parts['boundary']), intQ=list(parts['intQ']), upper_enclosure=list(ub), parameters=len(T),
               prescribed=int(mask.sum()))
    E.export_dir(f'{OUT}/{name}', P, U, st, expected=exp)
    vec = lambda v: [[int(c), str(x)] for c, x in sorted(v.items())]
    json.dump(dict(particular=vec(psi_b), states=[vec(v) for v in psis]), open(f'{OUT}/{name}/states.json', 'w', newline='\n'), indent=0)
    index[name] = dict(wolfram_model=wl_model, expected=exp)
    print(f"{name:14s} lower {lb:+.12f}   upper enclosure [{ub[0]:+.12f}, {ub[1]:+.12f}]   parameters {len(T)}   prescribed {int(mask.sum())}")


record('dirichlet_p1', E.Problem3D(box=[(0, 2), (0, 1), (0, 1)], N=(2, 1, 1), lam=lam, mu=mu, degree=1, g=gq),
       '<|"Lame" -> {3/2, 1}, "Domain" -> {{0, 2}, {0, 1}, {0, 1}}, "Cells" -> {2, 1, 1}, "Degree" -> 1, "BoundaryData" -> ' + GQ + '|>')
record('bodyforce_p2', E.Problem3D(box=[(0, 2), (0, 1), (Fr(1, 3), 1)], N=(1, 2, 1), lam=lam, mu=mu, degree=2, g=gq, b=(Fr(1, 2), Fr(-1, 3), Fr(1, 5))),
       '<|"Lame" -> {3/2, 1}, "Domain" -> {{0, 2}, {0, 1}, {1/3, 1}}, "Cells" -> {1, 2, 1}, "Degree" -> 2, "BodyForce" -> {1/2, -1/3, 1/5}, '
       '"BoundaryData" -> ' + GQ + '|>')
record('block_p1', E.Problem3D(box=[(0, 1)] * 3, N=(2, 2, 2), lam=lam, mu=mu, degree=1, bc=(DD, NN, NN), g=lambda X, Y, Z: (0.1 * X, 0 * Y, 0 * Z)),
       '<|"Lame" -> {3/2, 1}, "Domain" -> {{0, 1}, {0, 1}, {0, 1}}, "Cells" -> {2, 2, 2}, "Degree" -> 1, '
       '"TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}, {3, "Lower"}, {3, "Upper"}}, "BoundaryData" -> Function[{x, y, z}, {0.1 x, 0, 0}]|>')
record('ring_aniso', E.Problem3D(box=[(-1, 1), (Fr(1, 2), 2), (0, 1)], N=(1, 1, 2), D=Daniso, degree=2, bc=(NN, DD, NN), g=gq, b=(Fr(1, 3), Fr(1, 4), Fr(-1, 2))),
       '<|"StiffnessMatrix" -> ' + DWL + ', "Domain" -> {{-1, 1}, {1/2, 2}, {0, 1}}, "Cells" -> {1, 1, 2}, "Degree" -> 2, "BodyForce" -> {1/3, 1/4, -1/2}, '
       '"TractionFreeFaces" -> {{1, "Lower"}, {1, "Upper"}, {3, "Lower"}, {3, "Upper"}}, "BoundaryData" -> ' + GQ + '|>', randomise=True)
record('corner', E.Problem3D(box=[(0, 1)] * 3, N=(1, 1, 2), lam=Fr(1), mu=Fr(1), degree=1, bc=(ND, ND, ND), g=gq, b=(Fr(1), Fr(1), Fr(1))),
       '<|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {0, 1}, {0, 1}}, "Cells" -> {1, 1, 2}, "Degree" -> 1, "BodyForce" -> {1, 1, 1}, '
       '"TractionFreeFaces" -> {{1, "Lower"}, {2, "Lower"}, {3, "Lower"}}, "BoundaryData" -> ' + GQ + '|>')
record('cup_random', E.Problem3D(box=[(0, 1), (0, 1), (0, 2)], N=(1, 1, 2), lam=lam, mu=mu, degree=2, bc=(NN, NN, DN), g=gq, b=(Fr(0), Fr(0), Fr(-1))),
       '<|"Lame" -> {3/2, 1}, "Domain" -> {{0, 1}, {0, 1}, {0, 2}}, "Cells" -> {1, 1, 2}, "Degree" -> 2, "BodyForce" -> {0, 0, -1}, '
       '"TractionFreeFaces" -> {{1, "Lower"}, {1, "Upper"}, {2, "Lower"}, {2, "Upper"}, {3, "Upper"}}, "BoundaryData" -> ' + GQ + '|>', randomise=True)
json.dump(index, open(f'{OUT}/index.json', 'w', newline='\n'), indent=1)
print('wrote', OUT)
