#!/usr/bin/env python3
"""Small test vectors for the Wolfram Language port: inputs as raw doubles, expected outputs of the reference."""
import json, os, shutil
from fractions import Fraction as Fr
import numpy as np
import scalarnd as S
from ia import I, fr

rng = np.random.default_rng(11)
OUT = 'testvectors'
shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(OUT)
W4 = [Fr(1, 4), 0, Fr(-1, 2), 0, Fr(1, 4)]
WQ = [0, -1, Fr(1, 2), 0, Fr(1, 4)]                      # u^4/4 + u^2/2 - u  (convex)
index = {}


def record(name, P, U, cal, wl_model, extra=None, Qs=16):
    lb, parts = S.lower_bound(P, cal, Qs=Qs, return_parts=True)
    ub = S.upper_bound(P, U)
    exp = dict(lower=lb, flux_lower=S.flux_lower(P, cal), upper_enclosure=list(ub), Subintervals=Qs)
    if extra:
        exp.update(extra)
    S.export_field_dir(f'{OUT}/{name}', P, U, cal, expected=exp)
    index[name] = dict(wolfram_model=wl_model, expected=exp)
    print(f"{name:14s} lower {lb:+.12f}   upper enclosure [{ub[0]:+.12f}, {ub[1]:+.12f}]")


# 1. convex mode, 2D, rectangular grid, non-constant multilinear datum
P = S.Problem(box=[(0, 1), (0, 2)], N=(6, 5), k=Fr(1, 2), W=WQ, bc=[('D', 'D'), ('D', 'D')], g=lambda X, Y: 0.25 * X * (1 + Y))
out, U, cal = S.certify_convex(P)
record('convex2d', P, U, cal, '<|"Potential" -> Function[u, u^4/4 + u^2/2 - u], "Domain" -> {{0, 1}, {0, 2}}, "Cells" -> {6, 5}, "BoundaryData" -> Function[{x, y}, 0.25 x (1 + y)]|>')

# 2. convex mode, 3D, one natural face
P = S.Problem(box=[(0, 1), (0, 1), (0, 2)], N=(4, 3, 2), k=Fr(1, 2), W=WQ, bc=[('D', 'D'), ('D', 'N'), ('D', 'D')], g=lambda X, Y, Z: 0.5 * X + 0 * Y)
out, U, cal = S.certify_convex(P)
record('convex3d', P, U, cal, '<|"Potential" -> Function[u, u^4/4 + u^2/2 - u], "Domain" -> {{0, 1}, {0, 1}, {0, 2}}, "Cells" -> {4, 3, 2}, "BoundaryData" -> Function[{x, y, z}, 0.5 x], "NaturalFaces" -> {{2, "Upper"}}|>')

# 3. calibration mode, 2D wall (constant datum on each Dirichlet face, natural sides)
P = S.Problem(box=[(-3, 3), (0, 1)], N=(12, 2), k=Fr(1, 2), W=W4, bc=[('D', 'D'), ('N', 'N')], g=lambda X, Y: np.where(X < 0, -1.0, 1.0))
X = np.meshgrid(*P.axes, indexing='ij')
U = S.solve_candidate(P, np.tanh(X[0] / np.sqrt(2)))
cal = S.Calibration(P, -1, 1, 4, affine=False)
S.optimise_calibration(P, cal, Qo=4, taus=(1e-2, 1e-3, 1e-4), iters=300)
record('wall2d', P, U, cal, '<|"Potential" -> Function[u, (u^2 - 1)^2/4], "Domain" -> {{-3, 3}, {0, 1}}, "Cells" -> {12, 2}, "BoundaryData" -> Function[{x, y}, If[x < 0, -1, 1]], "NaturalFaces" -> {{2, "Lower"}, {2, "Upper"}}|>')

# 4. calibration mode, 3D, non-constant datum on the faces
P = S.Problem(box=[(-1, 1)] * 3, N=(4, 3, 2), k=Fr(1, 2), W=[1, 0, -2, 0, 1], bc=[('D', 'D')] * 3, g=lambda X, Y, Z: X + 0 * Y)
X = np.meshgrid(*P.axes, indexing='ij')
U = S.solve_candidate(P, X[0])
cal = S.Calibration(P, -1, 1, 4, affine=False)
S.optimise_calibration(P, cal, Qo=4, taus=(1e-2, 1e-3, 1e-4), iters=300)
record('ramp3d', P, U, cal, '<|"Potential" -> Function[u, (u^2 - 1)^2], "Domain" -> {{-1, 1}, {-1, 1}, {-1, 1}}, "Cells" -> {4, 3, 2}, "BoundaryData" -> Function[{x, y, z}, x]|>')

# 5. every code path at once: affine + spline, random coefficients, mixed boundary conditions, datum beyond the t-range
P = S.Problem(box=[(0, 2), (-1, 1)], N=(5, 4), k=Fr(3, 4), W=W4, bc=[('D', 'N'), ('D', 'D')], g=lambda X, Y: 0.7 * X * Y - 0.2)
U = rng.standard_normal(P.fixed.shape) * 0.3
U[P.fixed] = P.gval[P.fixed]
cal = S.Calibration(P, Fr(-1, 2), Fr(3, 4), 3, affine=True)
cal.unpack(0.2 * rng.standard_normal(cal.pack().size))
record('mixed2d', P, U, cal, '<|"Potential" -> Function[u, (u^2 - 1)^2/4], "Stiffness" -> 3/4, "Domain" -> {{0, 2}, {-1, 1}}, "Cells" -> {5, 4}, "BoundaryData" -> Function[{x, y}, 0.7 x y - 0.2], "NaturalFaces" -> {{1, "Upper"}}|>', Qs=4)

# primitives
prim = {}
Pq = S.Problem(box=[(0, 1)] * 2, N=(2, 2), k=Fr(1, 2), W=W4, bc=[('D', 'D')] * 2)
d = [-3.0, -0.75, -0.1, 0.0, 0.3, 1.25, 4.0]
prim['conjUpper'] = dict(W=[str(c) for c in W4], d=d,
                         plus_at_1=[float(x) for x in S.conj_upper(Pq, np.array(d), Fr(1), +1)],
                         minus_at_m1=[float(x) for x in S.conj_upper(Pq, np.array(d), Fr(-1), -1)],
                         plus_at_0=[float(x) for x in S.conj_upper(Pq, np.array(d), Fr(0), +1)],
                         note='rigorous upper bounds for sup{t d - W(t): t >= t0} (plus) and t <= t0 (minus); W = (t^2-1)^2/4')
K = 5
Cb = rng.standard_normal((7, K + 2))
t = np.array([-1.5, -1.0, -0.37, 0.2, 0.6, 1.0, 2.5])
r = S._spline_eval(Cb, K, Fr(-1), Fr(1), I(t))
prim['splineEval'] = dict(K=K, TRange=['-1', '1'], C=Cb.tolist(), t=t.tolist(), lo=r.lo.tolist(), hi=r.hi.tolist(),
                          note='clamped uniform quadratic B-spline with K knot intervals, one coefficient row per t')
index['primitives'] = prim

# the optimiser's objective and gradient (heuristic part; for porting optimise_calibration)
P = S.Problem(box=[(-3, 3), (0, 1)], N=(6, 2), k=Fr(1, 2), W=W4, bc=[('D', 'D'), ('N', 'N')], g=lambda X, Y: np.where(X < 0, -1.0, 1.0))
cal = S.Calibration(P, -1, 1, 3, affine=False)
x0 = 0.3 * rng.standard_normal(cal.pack().size)
cal.unpack(x0)
x0 = cal.pack()                                  # natural faces zeroed
S.optimise_calibration(P, cal, Qo=2, taus=(), iters=1)
fun = S.optimise_calibration.last_fun
opt = dict(wolfram_model='<|"Potential" -> Function[u, (u^2 - 1)^2/4], "Domain" -> {{-3, 3}, {0, 1}}, "Cells" -> {6, 2}, "BoundaryData" -> Function[{x, y}, If[x < 0, -1, 1]], "NaturalFaces" -> {{2, "Lower"}, {2, "Upper"}}|>',
           K=3, TRange=['-1', '1'], Qo=2, packing='C[0].ravel(), C[1].ravel(), S0[0].ravel(), S0[1].ravel() (C order; S0 entries are zero here)',
           x=x0.tolist(), cases=[])
for tau in (1e-1, 1e-2, 1e-3):
    f, g = fun(x0, tau)
    opt['cases'].append(dict(tau=tau, minus_J=float(f), grad_minus_J=g.tolist()))
cal.unpack(x0)
index['optimiser'] = opt
json.dump(index, open(f'{OUT}/index.json', 'w'), indent=1)
print('wrote', OUT)
