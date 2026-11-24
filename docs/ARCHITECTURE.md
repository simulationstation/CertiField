# CertiField architecture

This note is for engineers who want to audit a certificate, extend the paclet, or port a model. The guiding rule:
**every number that leaves CertiField as a bound is produced by arithmetic whose rounding direction is controlled, from
inputs that are exact.** Everything else (grids, multipliers, candidate profiles, calibration constants found by
optimisation) is a heuristic choice that can only make a bound less tight, never wrong.

## Layers

```
Models.wl              public API: CertifiedVariationalBound, VerifyCertificate, CertiFieldModel, presets
  |-- ScalarField1D.wl     1D: calibration lower bound, admissible competitor, exact / Taylor-model energy
  |-- ScalarFieldND.wl     2D/3D scalar fields on boxes: calibration-field lower bound, Q1 competitor, the two modes
  |-- SkyrmeLevelSet.wl    Skyrme lower bound: the level-set dynamic programme
  |-- SkyrmeHedgehog.wl    Skyrme upper bound: hedgehog trial maps and their certified energy
  |-- SkyrmeLevelSetDegree.wl  Skyrme lower bound for degree B (the level-set programme with B times the flux)
  |-- Minorants.wl         certified affine minorants of h_kappa (input to the dynamic programme)
CertifiedBound.wl      the result object: properties, directed decimal display, report
RigorousQuadrature.wl  adaptive interval quadrature (first order and Simpson with 4th-derivative remainder)
RigorousArithmetic.wl  exact enclosures of scalars; outward-rounded packed-array arithmetic; verified roots
```

## Two kinds of rigorous arithmetic

**Exact scalars.** τ-dependent constants of the Skyrme programme (sin, cos at grid levels, π-dependent floors) are exact
rationals: alternating Taylor series stopped with a remainder below 10⁻⁴⁰ (`sinBr`, `cosBr`), and π bracketed by
35-digit rationals. They are rounded to machine numbers in the safe direction by exact comparison (`fLo`, `fHi`).

**Outward-rounded arrays.** The dynamic programme works on packed arrays of machine numbers. Every `+ − × ÷` result is
widened by a relative 2⁻⁴⁹ (eight units in the last place) plus an absolute 10⁻³⁰⁰ (`rUp`, `rDn`). That encloses the exact
result for IEEE-754 operations with at most two roundings. Selections and comparisons use 0/1 masks, which involve no
rounding. Square and cube roots from the math library are not trusted. Each candidate root is checked by squaring or
cubing with directed rounding and nudged until the check passes (`sqrtLo`, `sqrtHi`, `cbrtLo`, `cbrtHi`). Packed arrays
cannot hold infinities, so ±10³⁰⁰ stands in for them and is only ever compared or selected.

**Interval quadrature.** `rigorousIntegralFn` encloses ∫f on each cell by the intersection of the first-order
enclosure h·f([l,r]) and Simpson's rule with the remainder −h⁵/2880·f⁗([l,r]), all in Wolfram Language `Interval`
arithmetic at 30 digits. Simpson is used only on cells where a caller-supplied predicate certifies that f is C⁴ (for
the 1D lower bound: V + c > 0 on the cell, checked by interval evaluation). Refinement is global: cells wider than
tolerance/(2·#cells) are bisected, so cells next to a singularity, where only first-order enclosures are valid, are
refined geometrically without blowing up the cost.

**Overflow.** An overflowing machine operation does not return an infinity in the Wolfram Language: the kernel repeats
it in software arithmetic and returns a number that is not a machine number, in an array that is no longer packed.
The Skyrme programme checks its cell bounds against the sentinel. `ScalarFieldND.wl`, whose inputs are a user's
arrays, tests the result of every interval operation (`ndChk`, a flag test) and stops as `"Inconclusive"`.

## The Skyrme lower bound (SkyrmeLevelSet.wl, Minorants.wl)

The analytic reduction (manuscript, Theorems 1, 3 and 8) shows E(U) ≥ W(τ₁, v) for some v ≥ 0, where W is the value
function of a relaxed level-set problem. The verifier builds, backwards over levels τᵢ, functions Λᵢ ≤ W(τᵢ, ·) that are
affine on the cells of a volume grid. Each transition between levels uses:

* certified affine minorants α + βD ≤ h_κ(D) = 2√((1+2D)(D²+2D+κ)). Validity is equivalent to the cubic
  P(D) = 8D³ + (20−β²)D² + (8+8κ−2αβ)D + 4κ−α² being ≥ 0 on [0, ∞), which `certifyMinorantMask` checks with
  outward-rounded interval evaluation at the cubic's critical point. The 1.8-million-entry table is rebuilt and
  certified in about 10 s, and re-certified whenever it is loaded from the disk cache.
* the second-order step lemma, condition [*]. Where the enclosure cannot confirm it, the step falls back to the
  always-valid first-order bound (84,771 fallbacks at μ = 1, Quick).
* the closed-form floors Ξ, Ξ_μ and the Bogomolny constant, and monotonicity of W in the volume.

The result is min_v Λ₁(v)/(12π²), rounded down. The volume grids in `Data/Skyrme` are heuristic, shipped as exact IEEE
doubles, and checked for consistency (increasing τ in (0, π), exactly representable steps) before use.

**Other degrees (SkyrmeLevelSetDegree.wl).** The degree of the map enters the analytic reduction in one place, the
flux through a level surface, which is B·4π sin²τ for a map of degree ±B (`docs/SKYRME_DEGREE_TWO.md`).
`levelSetLowerBoundDegree` is the verifier above with the three flux integrals of a step multiplied by B; for B = 1
it performs the same operations and returns the same numbers. The grids for B = 2 are in `Data/SkyrmeDegree`
(generator: `Reference/SkyrmeDegree`), together with stored claims that `tools/skyrme_degree_claims.wls` recomputes.
The model `SkyrmeDegreeTwo` pairs this lower bound with twice the certified energy of the hedgehog trial map, which is
the energy of two disjoint copies of it.

## The Skyrme upper bound (SkyrmeHedgehog.wl)

The trial map is a hedgehog with profile f = π − p, where p is a piecewise cubic Hermite function with rational data and
p(0) = 0 exactly. On the last cell f decreases linearly to 0, and it is identically 0 beyond. Any such data give an
admissible degree-one map, so the bound is valid whatever the data. The data (a numerical hedgehog) only affect
tightness. Energy:

* first cell: exact integral of a polynomial majorant (|sin p| ≤ |p|, 1 + cos p ≤ 2);
* middle cells: interval Simpson/first-order enclosures with a global tolerance of 10⁻¹⁰;
* last cell: exact integral of a polynomial majorant, with π bracketed.

## The one-dimensional engine (ScalarField1D.wl)

* **Model normalisation.** Inexact numbers in the model are replaced by the rationals they denote. Potentials must be
  built from smooth functions (`analyticExpressionQ` rejects `Abs`, `Max`, `Piecewise`, interpolation, …).
* **Calibration constant.** `calibrationFloor` proposes c = −inf V with `MinValue`, rounds it up to a rational, and
  *proves* V + c ≥ 0 on ℝ with `Resolve`. Trigonometric polynomials are proved via cos → c, sin → s, c² + s² = 1. On
  finite domains the optimal c is found by golden-section search (LB(c) is concave). Any admissible c is valid.
* **Lower bound.** 2√k·|∫√(V+c)| − c·L by interval quadrature, exact for monotone minimisers.
* **Competitor.** Cubic Hermite interpolation of the candidate at 2001 nodes, values and slopes rounded to multiples of
  10⁻²⁰, boundary values imposed exactly, constant beyond a window chosen where the candidate reaches the vacua.
* **Competitor energy.** The kinetic part is exact: ∫p′² = (6/5Δ² − Δ(m₀+m₁)/5 + 2/15(m₀²+m₁²) − m₀m₁/15)/h per cell.
  The potential part is exact for polynomial V (the integral of V∘p over a cell is precomputed symbolically as a
  polynomial in the cell data). Otherwise it uses a Taylor model of order 8 around the cell midpoint, with the t-integrals of
  the powers computed exactly and the Lagrange remainder enclosed on an exact enclosure of the cell's range.

## Scalar fields in two and three dimensions (ScalarFieldND.wl)

E[u] = ∫ k|∇u|² + W(u) on a box, u prescribed on some faces and free on the others, W a polynomial with rational
coefficients. Statement and proofs: `docs/SCALAR_ND.md`. The module is a port of the Python reference in
`Reference/ScalarND` (function map, conventions and port log: `Reference/ScalarND/PORTING.md`).

* **The certified problem.** The Dirichlet datum is the piecewise multilinear interpolant of the nodal boundary
  values on the grid; the nodal values are machine numbers taken as the exact binary numbers they are, and they are
  stored in the certificate.
* **Lower bound.** For any field S(x, t) = t σ₀(x) + Σ aₘ(x) bₘ(clamp t), with σ₀ and aₘ lowest-order Raviart–Thomas
  fields on the box grid and bₘ quadratic B-splines, every admissible u has
  E[u] ≥ ∫ S(x, g)·ν over the Dirichlet faces − ∫ sup_t (div_x S + |∂_t S|²/(4k) − W(t)).
  `ndLowerBound` evaluates the right-hand side for the given coefficients: the cell integral through the 2ⁿ vertices
  of each cell (the integrand is convex in each coordinate), the supremum over t on subintervals of the knot
  intervals with a curvature allowance, the two tails through certified tables of restricted convex conjugates
  (`conjUpper`, `conjCells`), and the boundary flux on sub-faces with a Lipschitz or Taylor penalty
  (`ndFluxLower`). Sample points, B-spline weights and W at the samples are exact rationals; everything else is
  outward-rounded array arithmetic; the sums over cells are exact rational sums.
* **Choosing the field (heuristic).** Convex mode: the affine part alone, the two-point flux of a cell-centred
  finite-volume solve (`ndFVSolve`, `ndConvexField`); this is the classical complementary-energy bound, tight for
  convex W, and takes seconds. Calibration mode: the spline part alone, found by L-BFGS on a smoothed copy of the
  certified functional, coarse to fine (`ndCalibrate`): the field is embedded exactly into the grid with twice as
  many cells (`ndProlongX`) or the spline space with twice as many knot intervals (`ndRefineT`), and the field with
  the best *certified* value is kept at every stage. A field from elsewhere can be supplied (`"Field"`).
* **Upper bound.** The candidate's nodal values define a piecewise multilinear function with the Dirichlet values
  imposed exactly. Its gradient energy is an exact quadratic form of the nodal values; its potential energy is
  enclosed by the midpoint rule with the remainder |c|/24 Σ hᵢ² sup|W″(u_h)|(∂ᵢu_h)² (`ndUpperBound`). With no
  candidate, or with `"RefineCandidate" -> True`, the package minimises the discrete energy by Newton's method
  (`ndSolveCandidate`).

## Certificates and verification

`CertifiedVariationalBound` returns `CertifiedBound[assoc]`. The reported bounds are machine numbers rounded outward
from the computed enclosures (lower down, upper up); the raw enclosure ends are kept under `"Certificate"`.
`VerifyCertificate` recomputes every rigorous step from the recorded inputs and checks that it reproduces bounds at
least as strong as the stored ones:

* Skyrme: the grid preset, the minorant table and the trial map;
* 1D: the proof of V + c ≥ 0, the calibration constant and the competitor data;
* 2D/3D scalar fields: the stored nodal boundary values (which must agree with the model's boundary data to
  rounding), the stored field and the stored nodal values of the candidate. The stored claims carry a margin of one
  part in 10⁹, because the machine root finder that places the conjugate tables need not return the same last
  bits on another machine; any placement gives a valid bound.

Nothing numerical is copied from the stored result.

## Adding a model

1. Write the analytic lower-bound argument and identify every inequality it needs.
2. Implement the bound with `RigorousArithmetic.wl` / `RigorousQuadrature.wl` primitives, keeping heuristics (grids,
   multipliers) separate from the certified evaluation.
3. Convert candidates into explicitly admissible competitors and compute their energy rigorously.
4. Register the model in `Models.wl`: `$CertiFieldModels`, `CertiFieldModelInfo`, a description with `"Theorem"`,
   `"AdmissibleClass"`, `"TrustedBase"`, `"NotClaimed"`, and a branch in `verifyCertificate`.
5. Add tests that enclose a known exact value and that check failure modes.
