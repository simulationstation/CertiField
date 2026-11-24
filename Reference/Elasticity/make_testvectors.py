#!/usr/bin/env python3
"""Small test vectors for the Wolfram Language model: inputs as raw doubles, expected outputs of the reference.
Each case stores gval.f64 (2 x nodes), U.f64 (2 x nodes), X.f64 (4 x Airy nodes), par.f64 and field.json;
index.json holds the model as a Wolfram Language expression and the expected numbers."""
import json
import os
import shutil
from fractions import Fraction as Fr

import numpy as np

import elasticity as E

rng = np.random.default_rng(17)
OUT = 'testvectors'
shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(OUT)
index = {}
Daniso = [[Fr(4), Fr(1), Fr(1, 2)], [Fr(1), Fr(3), Fr(-1, 4)], [Fr(1, 2), Fr(-1, 4), Fr(2)]]
gq = lambda X, Y: (0.1 * X - 0.15 * Y + 0.3 * X * Y, 0.25 * X + 0.2 * Y * Y - 0.1)
GQ = 'Function[{x, y}, {0.1 x - 0.15 y + 0.3 x y, 0.25 x + 0.2 y y - 0.1}]'


def record(name, P, wl_model, randomise=False):
    out, U, st = E.certify(P)
    if randomise:                                # any nodal values and parameters are admissible
        U = U + 0.2 * rng.standard_normal(U.shape)
        U[:, P.fixed] = P.gval[:, P.fixed]
        st = E.Stress(P, st.X + 0.2 * rng.standard_normal(st.X.shape), st.par + 0.2 * rng.standard_normal(st.par.shape))
    lb, parts = E.lower_bound(P, st, return_parts=True)
    ub = E.upper_bound(P, U)
    nodes, Tb, T = P.traces()
    exp = dict(lower=lb, boundary=list(parts['boundary']), intQ=list(parts['intQ']), upper_enclosure=list(ub), parameters=len(T),
               traction_free_nodes=len(nodes))
    E.export_dir(f'{OUT}/{name}', P, U, st, expected=exp)
    index[name] = dict(wolfram_model=wl_model, expected=exp)
    print(f"{name:12s} lower {lb:+.12f}   upper enclosure [{ub[0]:+.12f}, {ub[1]:+.12f}]   parameters {len(T)}")


record('dirichlet_p1', E.Problem(box=[(0, 2), (0, 1)], N=(3, 2), lam=Fr(3, 2), mu=Fr(1), degree=1, g=gq),
       '<|"Lame" -> {3/2, 1}, "Domain" -> {{0, 2}, {0, 1}}, "Cells" -> {3, 2}, "Degree" -> 1, "BoundaryData" -> ' + GQ + '|>')
record('bodyforce_p2', E.Problem(box=[(0, 2), (0, 1)], N=(2, 2), lam=Fr(3, 2), mu=Fr(1), degree=2, g=gq, b=(Fr(1, 2), Fr(-1, 3))),
       '<|"Lame" -> {3/2, 1}, "Domain" -> {{0, 2}, {0, 1}}, "Cells" -> {2, 2}, "Degree" -> 2, "BodyForce" -> {1/2, -1/3}, "BoundaryData" -> ' + GQ + '|>')
record('block_p1', E.Problem(box=[(0, 1), (0, 1)], N=(4, 4), lam=Fr(3, 2), mu=Fr(1), degree=1, bc=(('D', 'D'), ('N', 'N')), g=lambda X, Y: (0.1 * X, 0 * Y)),
       '<|"Lame" -> {3/2, 1}, "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> {4, 4}, "Degree" -> 1, "TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}}, '
       '"BoundaryData" -> Function[{x, y}, {0.1 x, 0}]|>')
record('threefree', E.Problem(box=[(Fr(1, 3), 1), (0, 1)], N=(3, 2), D=Daniso, degree=2, bc=(('D', 'N'), ('N', 'N')), g=gq, b=(Fr(1, 5), Fr(-2, 3))),
       '<|"StiffnessMatrix" -> {{4, 1, 1/2}, {1, 3, -1/4}, {1/2, -1/4, 2}}, "Domain" -> {{1/3, 1}, {0, 1}}, "Cells" -> {3, 2}, "Degree" -> 2, '
       '"BodyForce" -> {1/5, -2/3}, "TractionFreeFaces" -> {{1, "Upper"}, {2, "Lower"}, {2, "Upper"}}, "BoundaryData" -> ' + GQ + '|>')
record('random', E.Problem(box=[(-1, 1), (Fr(1, 2), 2)], N=(2, 3), lam=Fr(3, 2), mu=Fr(1), degree=2, bc=(('D', 'D'), ('N', 'N')), g=gq, b=(Fr(1, 3), Fr(1, 4))),
       '<|"Lame" -> {3/2, 1}, "Domain" -> {{-1, 1}, {1/2, 2}}, "Cells" -> {2, 3}, "Degree" -> 2, "BodyForce" -> {1/3, 1/4}, '
       '"TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}}, "BoundaryData" -> ' + GQ + '|>', randomise=True)
record('corner', E.Problem(box=[(0, 1), (0, 1)], N=(2, 2), lam=Fr(1), mu=Fr(1), degree=1, bc=(('N', 'D'), ('N', 'D')), g=gq, b=(Fr(1), Fr(1))),
       '<|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> {2, 2}, "Degree" -> 1, "BodyForce" -> {1, 1}, '
       '"TractionFreeFaces" -> {{1, "Lower"}, {2, "Lower"}}, "BoundaryData" -> ' + GQ + '|>')
json.dump(index, open(f'{OUT}/index.json', 'w'), indent=1)
print('wrote', OUT)
