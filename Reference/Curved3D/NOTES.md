# Guaranteed stiffness bounds for a 3D elastic body with a curved free surface

Status: first working attempt, Python reference only. The argument is checked by its author only; the stored
fields are also checked by a second program in exact arithmetic. No novelty is claimed (see "What is known").

## The body and the question

- A block (-1,1)^3 of isotropic elastic material, Young's modulus 1, Poisson ratio 3/10.
- A spherical cavity of radius 1/2 at the centre. Its surface is free (no traction).
- The faces y = +-1 and z = +-1 are free.
- The faces x = +-1 are rigid grips moved apart by delta. Two variants:
  - clamped: the material sticks to the grip, u = (+-delta/2, 0, 0);
  - frictionless: only u_x = +-delta/2 is prescribed, the material may slide.
- Stiffness K: the minimal elastic energy is K delta^2 / 2. (Force on a grip = K delta.)

The program returns numbers K_low <= K <= K_up that hold for the exact curved body, with rounding errors enclosed.

## Results

| Grips | Certified bounds | Gap | Grids (stress / displacement) |
|---|---|---|---|
| clamped | 1.817322 <= K <= 1.833487 | 0.89 % | 24 / 24 cells per half-edge |
| frictionless | 1.745105 <= K <= 1.758277 | 0.76 % | 24 / 24 cells per half-edge |

Without the cavity the frictionless value is exactly 2 (uniform tension). So the cavity removes between 12.0 % and
12.8 % of the stiffness, proven.

Earlier, coarser certificates of the same quantities (all valid, each run is independent):

| Grips | Grid | Certified bounds | Gap | Time on 2 cores |
|---|---|---|---|---|
| clamped | 8 | 1.751991 <= K <= 1.839118 | 4.97 % | 37 s |
| clamped | 16 | 1.805106 <= K <= 1.834149 | 1.61 % | 25 min |
| frictionless | 8 | 1.689106 <= K <= 1.762080 | 4.32 % | 40 s |
| frictionless | 16 | 1.734830 <= K <= 1.758619 | 1.37 % | 25 min |

The grid-24 runs took about 25 minutes per body on one core (19 for the stress field, 5 for the displacement).
`verify.py` re-read the stored grid-24 fields and, at a shallower octree depth (3 and 5 instead of 4 and 6), confirmed
1.817059 <= K <= 1.834215 (clamped) and 1.744878 <= K <= 1.758911 (frictionless).

## Why the bounds are true

Notation: Omega is the body, C the elasticity tensor, eps(v) the strain of a displacement v.

**Upper bound.** For every displacement v that satisfies the grip conditions, K delta^2/2 <= energy of v. The program
takes a spline v on a box grid that covers the body, and integrates its energy density over a set Omega+ that contains
Omega. The density is non-negative, so the integral over Omega+ is at least the integral over Omega.

**Lower bound.** Let tau be a symmetric stress field with div tau = 0, zero traction on every free surface (the sphere
included) and, for frictionless grips, zero shear traction on the grips. Then for every admissible v

    energy(v) >= int_Omega tau:eps(v) - (1/2) int_Omega tau:C^-1:tau  =  delta F - (1/2) Q,

where F is the force that tau transmits through a grip and Q = int_Omega tau:C^-1:tau. Scaling tau gives K >= F^2 / Q.
Again Q is replaced by the integral over Omega+, which is larger.

**How tau is made exactly admissible on the curved surface.** tau = tau_p + inc(rho^2 Psi), where
- inc is the Beltrami operator, tau_ij = e_ikl e_jmn d_k d_m Phi_ln. For any twice weakly differentiable symmetric Phi
  the field inc(Phi) is symmetric and divergence-free.
- rho(x) = x^2 + y^2 + z^2 - 1/4 vanishes exactly on the sphere. A potential rho^2 Psi vanishes there together with its
  gradient, so its second derivatives on the sphere are 2 Psi (grad rho)(grad rho), and inc of that has zero traction:
  e_jmn (d_j rho)(d_m rho) = 0. This holds for every Psi, so no curved elements and no approximation of the sphere occur.
- Psi is a tensor field of C^1 quadratic B-splines (six components). On the free flat faces its two outer coefficient
  layers are zero (value and normal derivative vanish, hence zero traction). On the mirror planes each component is even
  or odd as the tensor symmetry demands.
- A field inc(Phi) whose potential vanishes to second order on the free faces carries no net force through a grip
  (integrate the three terms of tau_xx by parts over the grip face). The force is carried by tau_p, an axial stress
  c (1 - indicator of the square column |y|,|z| < 1/2) e_x e_x. It is divergence-free because it does not depend on x,
  it is zero on the column that contains the cavity, and its normal traction on the lateral faces is zero. Its force is
  exactly the cross-section area 4.
- Frictionless grips need zero shear on x = +-1. For that the weight is made flat at the grip: rho uses
  eta(x) = x^2 for |x| <= 3/4 and x^2 - 4 (|x| - 3/4)^2 beyond, which is C^1, has eta'(1) = 0 and does not change the
  zero set. Then even/odd ties of Psi at the grip make it a mirror plane of the potential.

**Geometry.** Each grid cell is split into an octree of boxes. A box is dropped only if its farthest corner lies in the
closed ball, decided in exact integer arithmetic. All other boxes are kept whole. So Omega+ contains Omega.

**Integration.** On every box both integrands are polynomials (degree at most 12 per variable for the stress, 4 for the
strain), so the 7-point and 3-point Gauss rules are exact. The nodes and weights are enclosed in rational intervals of
width below 1e-30; all field values are computed in outward-rounded interval arithmetic; sums are rounded up.

**What is trusted.** Only the evaluation. The solver that finds Psi and v works in ordinary floating point and may be
wrong or sloppy: a worse field gives a weaker bound, never a false one. `verify.py` re-checks a stored pair (Psi, v)
without the solver, including an exact test that the stored coefficients satisfy the ties. `exact_check.py` re-checks
the same stored fields with a second, separate implementation in exact rational arithmetic (see below). The program
refuses set-ups for which the construction is not valid, for example a cavity that sticks out of the column, or
frictionless grips without the flattened weight.

## Tests (`python3 tests.py`, about 5 minutes)

1. The rational Gauss rules agree with numpy.
2. For random Psi: div tau = 0 (finite differences), zero traction on the sphere (1e-15), on the free faces and zero
   shear on mirror planes and frictionless grips.
3. Two independent code paths give the same stress; 7 and 9 Gauss points give the same integral (so the degree count
   is right); the interval bounds sit just above the float values.
4. Green identity with random fields: int_Omega tau:eps(v) equals force times grip displacement. The particular part
   gives 0.5 exactly, the potential part tends to 0 as the octree deepens. This can only hold if equilibrium, all
   boundary conditions and the zero-resultant property are right.
5. Known answers: a solid cube between frictionless grips, exact K = 2: bounds 1.99999999999923 and 2.00000000000049.
   A round bar of radius 3/4 (curved free surface, exact K = pi R^2/2 = 0.883573): the bounds bracket it at every depth.

## Exact rational check (`exact_check.py`)

A second program builds the stress on every grid cell as a polynomial with rational coefficients (the stored floats are
exact binary fractions) and tests the conditions of the lower bound as exact identities, not to a tolerance. It does not
use the evaluation code of `curved3d.py`.

For both final (grid 24) certificates, on all 13 824 cells (9 minutes each):

- (a) div tau = 0 in every cell;
- (b) the normal traction is continuous across all 39 744 faces between cells;
- (c) in all 325 cells cut by the sphere, tau.grad(rho) is divisible by rho, so the traction on the sphere is zero as a
  polynomial identity;
- (d) zero traction on the free faces, zero shear on the mirror planes and on the frictionless grip;
- (e) the potential part transmits exactly zero force through the grip, the whole field exactly 1 per octant;
- (f) on six boxes the exact integral lies below the interval bound, by at most 8.5e-13 (relative);
- (g) the boxes cover the body: full enumeration at depths 3 and 4, about 1.7 million random points next to the sphere
  at depths 5 and 6;
- (h) the displacement is continuous across all faces and meets the grip and mirror conditions exactly.

For both grid-8 certificates, in addition, both bounds were recomputed completely in exact rational arithmetic
(5 minutes each):

| Grips | Exact lower bound | Interval result | Exact upper bound | Interval result |
|---|---|---|---|---|
| clamped | 1.751991202691544 | 1.751991202687701 | 1.839118071762823 | 1.839118071815897 |
| frictionless | 1.689105784678506 | 1.689105784674822 | 1.762080379392570 | 1.762080379443421 |

The interval results are on the safe side of the exact values and agree with them to 11 digits.

Not covered by this check: the mathematical argument itself ("Why the bounds are true"), and for grid 24 the complete
sums of the interval evaluation (six boxes were recomputed exactly, not all).

## Limits

- One family of shapes: a box cut by a surface rho = 0 with rho a quadratic polynomial (sphere, cylinder). Loads only
  through flat grips; the curved surface must be free; no body force; three mirror planes are used.
- A hole that goes right through the body (a plate with a round hole) needs one more ingredient: the share of the force
  and moment carried by each ligament must be free, which the present potential cannot change. Not done.
- The fields are computed by this program on its own grid. It does not read the mesh or the solution of another
  solver. It can still serve as a reference value against which another solver's stiffness is compared;
  `ndsolve_compare.wls` does that for NDSolve's finite element solver (written, not yet run).
- Cost: the lower bound converges like h^1.6. Halving the gap costs roughly 8 times the work.
- The clamped problem has a weak stress singularity along the grip edges, which slows both bounds.
- Not ported to the Wolfram Language.
- The argument has been checked by its author only. `prompt.md` asks for a second check.

## What is known

The ingredients are classical or published: the Prager-Synge / complementary energy bounds; Beltrami stress functions
and their completeness conditions (Gurtin); spline stress functions for the complementary energy; weighted B-splines
with an Airy stress function for plates with free holes (the two-dimensional relative of the weight used here);
verified finite element computations with interval arithmetic on polygonal domains (Liu, Oishi, Nakao, Plum).
Two web searches found no paper that combines them for a curved free surface in 3D with rounding control; that is not
proof that none exists.

## Files

- `curved3d.py` the method; `ia.py` interval arithmetic; `verify.py` re-evaluation of stored fields;
  `exact_check.py` the exact rational check; `tests.py` the tests; `certificates/` the stored fields; `logs/` the run logs.
- `ndsolve_compare.wls` the same body in NDSolve, set against the proven range (not yet run);
  `prompt.md` instructions for a second check.
- Runs (arguments: grid, grips, octree depths):
  - `python3 run_cert.py 8 clamped 2 3 3 5` both bounds on one grid (37 s);
  - `python3 run_dual.py 24 frictionless 2 4` lower bound only; `python3 run_primal.py 24 frictionless 4 6` upper bound only;
  - `python3 verify.py 24 frictionless certificates/psi_frictionless_N24.npz certificates/v_frictionless_N24.npz 4 6`;
  - `python3 exact_check.py frictionless certificates/psi_frictionless_N24.npz --v certificates/v_frictionless_N24.npz`;
  - `python3 exact_check.py clamped certificates/cert_clamped_N8_L2335.npz total 3 5` both bounds in exact arithmetic.
- Needs Python 3 with numpy and scipy. Nothing else.
