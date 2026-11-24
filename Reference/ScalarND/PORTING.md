# Porting notes: `scalarnd.py` → `Kernel/ScalarFieldND.wl`

## 1. State of the port

`Kernel/ScalarFieldND.wl` was first written without a Wolfram kernel and then executed, corrected and completed with
Wolfram Engine 15.0 (Windows, x86-64) on 23 November 2026. Both modes now run in the Wolfram Language without the
Python:

* the model specification, the Q1 competitor (a Newton solve for the candidate) and its certified energy;
* the certified lower bound for any field (affine part, spline part, or both), including the boundary flux for
  non-constant multilinear data and the conjugate tables;
* convex mode: finite-volume solve and the affine field;
* calibration mode: the optimiser, the two embeddings and the coarse-to-fine driver (Section 6);
* import of fields written by the Python reference (`export_field_dir`) and of its stored certificates
  (`NAME.json` + `NAME.npz`);
* registration in `Models.wl`, a `VerifyCertificate` branch, tests.

Tests: `wolframscript -file tools/run_tests.wls --nd` (63 tests, about a minute) and `--nd-slow` (10 tests: the
acceptance runs of Section 6 and every stored reference certificate, about nine minutes). The existing fast suite
passes unchanged.

The Python reference remains the specification of the certified formula: the two implementations evaluate the same
right-hand side of Theorem 1 and agree to rounding wherever they are meant to (Section 4). They differ in the points
listed in Section 5, which never weaken a bound. What the first execution found is logged in Section 7, and
observations about the reference and the note are in Section 8.

## 2. Function map

| Python (`scalarnd.py`) | Wolfram Language | Rigorous? |
|---|---|---|
| `Problem` | `scalarFieldNDSpec` | exact data |
| `solve_candidate` (L-BFGS) | `ndSolveCandidate` (Newton with sparse matrices), called by `ndCompetitor` | no |
| `upper_bound` | `ndUpperBound` | yes |
| `Calibration` | field association `<|"TRange", "K", "Affine", "S0", "C"|>`, `ndPackField`, `ndRespectsNaturalFaces` | data |
| `conj_upper`, `conj_cells`, `tail_inf` | `conjUpper`, `conjCells`, `tailInf` | yes |
| `conj_float` | `conjArgmax` (placement of the fine partition only) | no |
| `_spline_eval` | `splineEval` | yes |
| `flux_lower` | `ndFluxLower` | yes |
| `lower_bound` | `ndLowerBound` | yes |
| `fv_energy`, `fv_solve` (L-BFGS) | `ndFVSolve` (damped Newton with sparse matrices) | no |
| `affine_from_candidate` | `ndConvexField` | no |
| `optimise_calibration` | `ndCalSetup`, `ndCalValueGradient`, `ndLBFGS`, `ndWolfe`, `ndCalStage` | no |
| `refine_t`, `prolong_x` | `ndRefineT`, `ndProlongX` | no |
| `certify_calibration` | `ndCalibrate` (with `ndAutoLevels`, `ndAutoTRange` for the automatic choices) | no |
| `claim_strings` | the margin at the end of `certifyScalarFieldND` | exact |
| `export_field_dir` | `ndImportField`, `ndImportArray` | data |
| `load_certificate`, `verify_certificate` | `ndImportCertificate`, `ndVerifyReferenceCertificate` (reference certificates); `verifyScalarFieldND` (the paclet's own) | yes |

## 3. Conventions

* Arrays are indexed (x₁, x₂, x₃) with the last index fastest (NumPy C order = `Flatten`/`ArrayReshape` order).
  Raw files are little-endian float64 in that order.
* Node arrays have dimensions `Cells + 1`; cell arrays `Cells`; the face array of direction i has dimensions `Cells`
  with entry i increased by one, and holds the component of the field along +eᵢ on each face orthogonal to eᵢ.
  The cell with multi-index c has lower face c and upper face c + eᵢ in that array.
* Spline part: on the faces of direction i, `C[i]` has a last dimension of length K + 2 holding B-spline
  coefficients c₀ … c_{K+1} (Python, 0-based). On knot interval κ = 0 … K−1 with local coordinate θ ∈ [0, 1],

      s(θ) = c_κ (1−θ)²/2 + c_{κ+1} (−2θ² + 2θ + 1)/2 + c_{κ+2} θ²/2,
      ds/dt = [ (c_{κ+1} − c_κ)(1−θ) + (c_{κ+2} − c_{κ+1}) θ ] / Δ,      d²s/dt² = (c_κ − 2c_{κ+1} + c_{κ+2}) / Δ².

  The value at knot j is (c_j + c_{j+1})/2. Outside [t_lo, t_hi] the spline is constant.
  Greville abscissae ξ_m = t_lo + (m − ½)Δ satisfy Σ ξ_m b_m(t) = t, so inside the range the affine part t σ₀ is
  the same as adding σ₀ ξ_m to the coefficients (this is how both evaluators treat it).
* Outward normal sign: −1 on the lower face, +1 on the upper face.
* Exact data: the box, k, the coefficients of W, t_lo, t_hi are rationals. Nodal values and field coefficients are
  machine numbers taken as the exact binary numbers they are.
* Rounding. The Python widens each operation by one ulp; the paclet's `rUp`/`rDn` widen by 2⁻⁴⁹ relative (eight
  ulps) plus 10⁻³⁰⁰. Results therefore agree to about 10⁻¹³ relative, not bit for bit. Sums over cells are exact
  rational sums in the port (`exactTotal`) and `math.fsum` widened by two ulps in the reference.
* Overflow and underflow. NumPy returns an infinity on overflow, which the reference catches with `isfinite`. The
  Wolfram kernel instead repeats an overflowing machine operation in software arithmetic, without a message; the
  result is not a machine number and an array holding it is no longer packed. The enclosure argument does not cover
  that, so every interval helper of the port tests its result (`ndChk`, a flag test) and stops as `"Inconclusive"`.
  For this to be meaningful, inputs are packed before use and `Developer`ToPackedArray` is never applied to the
  result of an arithmetic operation. Underflow to subnormal numbers happens routinely (exact zeros widened to
  10⁻³⁰⁰, then multiplied by 2⁻⁴⁹); IEEE gradual underflow is exact to 2⁻¹⁰⁷⁵, inside the absolute widening, and the
  kernel's warning about it is suppressed (`ndQuiet`).
* Stored claims. Both implementations weaken the computed bounds by one part in 10⁹ before storing them
  (`claim_strings` in the Python; the end of `certifyScalarFieldND` in the port). Verification recomputes the
  bounds, and the floating-point root finder that places the conjugate tables (`conj_float`, `conjArgmax`) is not
  bitwise reproducible across machines. Any placement is valid; the margin keeps the comparison from failing on
  the last digits.
* The datum. The paclet's certificate stores the nodal boundary values (`"BoundaryNodalValues"`).
  `VerifyCertificate` uses them as the datum and requires them to agree with the model's boundary data to 10⁻¹²,
  because a function of the coordinates need not evaluate to the same last bits on another machine.

## 4. Test vectors (`testvectors/`)

`index.json` has one entry per case with the model as a Wolfram Language expression (`wolfram_model`) and the
reference's outputs (`expected`). Each case directory holds `field.json` and raw arrays `gval.f64`, `U.f64`,
`S0_i.f64`, `Ci.f64`.

| Case | What it exercises |
|---|---|
| `convex2d` | affine field, rectangular 6×5 grid, non-constant multilinear datum |
| `convex3d` | affine field, 4×3×2 grid, one natural face |
| `wall2d` | spline field, constant datum on each Dirichlet face, natural sides |
| `ramp3d` | spline field in 3D, datum varying on the faces (sub-faces, Taylor bound) |
| `mixed2d` | affine + spline with random coefficients, mixed boundary conditions, datum beyond the t-range, k = 3/4 |
| `primitives` | `conjUpper` at seven slopes for three (t₀, side); `splineEval` at seven points including both clamps |
| `optimiser` | value and gradient of the smoothed objective at a random point for three τ |

`expected` contains `lower` (with `Subintervals` = Q_s used), `flux_lower`, `upper_enclosure`.

How the port compares (`Tests/ScalarFieldND.wlt`):

* `upper_enclosure`: 10⁻⁹ relative in all five cases (measured: 10⁻¹³).
* `flux_lower` and `lower` for `convex2d`, `convex3d`, `wall2d`: 10⁻⁹ relative; 10⁻⁶ for lower bounds of fields with
  an affine part (the conjugate table is interpolated and the two implementations place its nodes slightly
  differently). Measured: 10⁻¹³ to 10⁻¹⁴.
* `flux_lower` for `ramp3d` and `mixed2d` (datum varying on the faces): the port is tighter than the reference by
  design (Section 5a), so the test is two-sided against the *exact* flux, which the test computes in rational
  arithmetic by integrating the clamped spline in closed form: reference − 10⁻⁹ ≤ port ≤ exact. The cell integral
  (lower bound minus flux) still agrees with the reference to 10⁻⁹ (10⁻⁶ with an affine part).
* `primitives`: 10⁻⁹; `optimiser`: 10⁻⁹ (measured: 10⁻¹⁵ for the value and the gradient).

## 5. Where the port differs from the reference

None of these weakens a bound; (a) makes the port's bound tighter.

a. **Knot intervals in the boundary flux.** Lemma 5 takes the Lipschitz and curvature constants over the knot
   intervals met by the datum on a sub-face. The reference finds them as floor((g − t_lo)/Δ) − 1 … + 1 in floating
   point, a superset by one interval on each side. The port decides them by exact comparisons of machine numbers
   with enclosures of the knots: interval κ is included when (lower enclosure of t_κ) < g_max and (upper enclosure
   of t_{κ+1}) > g_min, still a superset of the intervals whose interior meets [g_min, g_max], and empty when the
   datum lies outside the t-range (where the clamped spline is constant). Reason: with the paclet's wider rounding
   the reference's rule picks up one more interval whenever a datum value sits exactly on a knot, which made the
   port *weaker* than the reference by 6×10⁻⁸ (2D ramp) and 3×10⁻⁷ (3D ramp), enough to fail re-verification of
   those two stored certificates. With the exact rule the port's lower bounds are higher than the reference's on
   the same fields: +5×10⁻⁶ relative (`ramp2d`), +3×10⁻⁵ (`ramp3d_32`), +4×10⁻⁷ (`saddle2d`), and all eleven stored
   certificates verify.
b. **Overflow is detected** (Section 3).
c. **Sums are exact** (Section 3).
d. **Candidate.** `ndSolveCandidate` minimises the same discrete energy (Q1, two-point Gauss) by Newton's method
   with Armijo steps. For convex W the Newton system is solved as it is (conjugate gradients with an incomplete
   factorisation above 5000 unknowns). For other W the direction comes from conjugate gradients truncated at the
   first direction of negative curvature, which is always a descent direction: plain Newton steps converge to
   saddle points as readily as to minima, and the first version of the port did exactly that on the saddle example
   (energy 13.02 instead of 10.95). The package's own start is the finite-volume solution for convex W and
   otherwise the boundary function evaluated at all nodes plus a fixed perturbation of relative size 10⁻³, so that
   a symmetric start does not stay on a symmetric saddle point. The certified upper bounds agree with the
   reference's to nine or ten digits on all eleven examples.
e. **Optimiser.** `ndLBFGS` is a plain L-BFGS (20 correction pairs, two-loop recursion) with a strong-Wolfe line
   search (bracketing and zoom), not SciPy's L-BFGS-B. The objective, the schedule of (τ, iterations) and the
   limit of three evaluations per iteration are the reference's. Fields differ in the last digits, the certified
   bounds by parts in 10⁶.
f. **Automatic choices** (new in the port, heuristic): convex mode if W″ ≥ 0 on the line (decided exactly),
   calibration mode otherwise; the t-range is the hull of the Dirichlet data and of the critical points of W,
   rounded outward to a binary fraction; the scale of τ is max(1, range of W on the t-range); the levels halve the
   grid up to three times and K down to 8 and then double them alternately (grid first). For the 64×64, K = 32 ramp
   this is exactly the reference's list of levels.
g. **Certificates** store the nodal boundary values (Section 3).

## 6. The optimiser (calibration mode)

Specification: `optimise_calibration`, `refine_t`, `prolong_x`, `certify_calibration` in `scalarnd.py`.

* Unknowns: the coefficients `C[i]` (spline part alone; `affine=False`). Natural faces are masked to zero.
* Objective (to maximise; concave): flux − Σ_cells |c| 2⁻ⁿ Σ_vertices τ·log-mean-exp over
  { Φ̃(κ, j) : κ < K, j = 0 … Q_o } ∪ { two tails }, with Φ̃ = Φ + (w²/8)(sup W″ − D″)⁺, w = Δ/Q_o. The curvature
  term is essential: without it the optimiser hides peaks between samples (see `check_theorem.py`, where this is
  used as an attack).
* Flux: linear in the coefficients, 2-point Gauss rule on each Dirichlet face.
* Solver: L-BFGS, 20 correction pairs. Schedule of (τ, iterations) in `certify_calibration` and `$ndSchedule`:
  first level (10⁻², 300), (10⁻³, 400), (10⁻⁴, 600); intermediate levels (10⁻⁴, 300), (3·10⁻⁵, 400); last level
  (10⁻⁴, 300), (3·10⁻⁵, 600), (10⁻⁵, 600). Every τ is multiplied by `scale`, of the order of the barrier height of
  W (1 for (u² − 1)²/4 and (u² − 1)², 4 for 4(u² − 1)² in the examples).
* Multilevel: `prolong_x` (cells doubled in every direction: copy on old faces, average on new faces in the normal
  direction, repeat in the tangential directions) and `refine_t` (K doubled: c′_{2j} = (3c_j + c_{j+1})/4,
  c′_{2j+1} = (c_j + 3c_{j+1})/4). Both are exact embeddings of the field, so the right-hand side of Theorem 1 is
  unchanged. The certified value can still move a little, because the allowances of Lemma 3 and Lemma 5 are
  evaluated on a different partition; the driver therefore re-evaluates it after each embedding.
* The field with the best *certified* value is kept after every stage within a level.
* In the port the objective is evaluated on matrices (cells × samples): the divergence coefficients of every cell
  are gathered once (cells × (K + 2)), multiplied by the three B-spline tables, and the 2ⁿ soft-max passes share
  them. One evaluation takes 0.13 s for the 64×64, K = 32 ramp (283 000 unknowns).

Acceptance (run by `--nd-slow`, numbers of 23 November 2026):

| | Reference (stored claim) | Wolfram Language |
|---|---|---|
| test vector `optimiser`: −J and its gradient at the stored point | | equal to 10⁻¹⁵ (required: 10⁻⁹) |
| `wall2d`, 96×4, K = 32 | 0.9442106791 ≤ inf E ≤ 0.9451112063 | 0.9442107497 ≤ inf E ≤ 0.9451112063 |
| exact wall value 0.9449674266 | enclosed | enclosed |
| `ramp2d`, 64×64, K = 32 | 8.275401140 ≤ inf E ≤ 8.291701687 | 8.275406985 ≤ inf E ≤ 8.291701687 |

Both lower bounds are within 0.1 % of the reference's (the criterion), in fact within 10⁻⁶, and the upper bounds
agree to ten digits. The fields come from two different optimisers, so the lower bounds are not expected to agree
further; an earlier run of the port, before its objective was rearranged, gave 8.275416 for the ramp.

Things tried in the reference that did not work: optimising the affine and spline parts together (slow; the parts
overlap inside the range), and freezing the affine part at the finite-volume solution for non-convex W (the tails
bind and the bound collapses).

## 7. Port log: what the first execution found

The draft passed 35 of its 37 tests on the first run. The items of the former list "where the draft is most likely
wrong" that turned out to be right: the slicing helpers, threading of arrays of different depth, `splineEval`
(gather, both clamps, the last piece), the bracket repair in `conjCells`, the `KroneckerProduct` order and the
boundary source in `ndFVSolve`, associations with list keys, `Return` inside `If`, and the JSON import.

What was changed, and why:

1. **Display of numbers of modulus 10 or more** (`Kernel/CertifiedBound.wl`, `decimalString`; in the paclet since
   0.2.0, not in the draft). `IntegerString[n, 10, len]` keeps only the last `len` digits, so 10.6792 was printed as
   "0.6792", 1234570 as "0", and a relative gap of 124 % as "5%", in `"Summary"`, `"Theorem"` and `"Report"`. The
   stored bounds were never affected. No existing example reached 10; `saddle2d` would have. Fixed (the digit
   string is padded, never shortened); test `decimal-strings-keep-leading-digits`.
2. **Underflow warnings.** With unpacked arguments the kernel printed `General::munfl` for subnormal products. The
   values were right. Arguments are now packed and the warning is suppressed (Section 3).
3. **Overflow was not handled** (Section 3). Test `nd-machine-overflow-is-refused`.
4. **Boundary flux weaker than the reference** on knot-aligned data (Section 5a). The five small test vectors did
   not show it; re-verifying the stored certificates did.
5. **`ndFVSolve` did not stop** on fine grids with zero data: its test, max|gradient| < 10⁻¹²(1 + max|b|), is below
   rounding level there (128×128: 60 iterations of failing line searches, 6 s). It now stops when the Newton step
   is at rounding level or the energy no longer decreases.
6. **Candidate.** `ndCompetitor[spec, Automatic]` took node means of the finite-volume solution. It now runs
   `ndSolveCandidate` (Section 5d), also on a user's candidate when `"RefineCandidate" -> True`. Only the boundary
   values of the boundary function are required to be real numbers (it may be singular inside the box).
7. **Speed.** `tailInf` is remembered (it was recomputed, one second each, at every call of `ndLowerBound`);
   `ndFieldQ` no longer unpacks the field (`ndPackField`); the two Newton solvers use conjugate gradients above
   5000 unknowns (the 48³ example went from 56 s to 6 s).
8. **Verification** recomputed the datum from the model's function (Section 3, last item). A certificate whose
   stored boundary values or field were changed is rejected (two new tests).
9. **Tests.** `nd-rough-candidate-gives-honest-wide-gap` used 0.26 sin πx sin πy, guessed without a kernel: its
   energy is positive, a gap of 124 %, outside the expected 1 % to 50 %. The code was right; the amplitude is now
   0.08 (gap 6 %), and the 0.26 case is kept as a test of `"RefineCandidate"`. The flux tests for varying data are
   two-sided against the exact flux (Section 4), also for forty random problems with data on the knots, between
   them and beyond the range. New tests: exactness of the two embeddings, the optimiser's test vector, calibration
   mode through the public interface in 2D and 3D, four stored reference certificates, an `NDSolve` candidate, the
   saddle-point case of Section 5d, single-cell grids, overflow and tampering.

## 8. Observations about the reference and the note

Nothing was found that makes a bound wrong. The three check scripts pass on Windows (Python 3.11.9, NumPy 2.4.6,
SciPy 1.17.1; logs in `logs/`), and the port, a second implementation with different arithmetic, recomputes all
eleven stored certificates. The proofs of the note were read during the port and no gap was noticed; that is not
an independent check of them.

1. **The boundary flux of the reference is looser than Lemma 5 allows** (Section 5a): `flux_lower` takes one knot
   interval more on each side than the datum meets. Valid, but it costs up to 3×10⁻⁵ of the bound in the examples
   and more on coarse grids (0.07 % of the flux in test vector `ramp3d`, 10 % in `mixed2d`, whose field is random). The exact rule of
   the port could be taken over (`ka`, `kb` from comparisons with `frac_enclosure` of the knots). Not changed in the
   reference, so that its stored test vectors and claims stay as they were produced.
2. **The conjugate tables have a slack of order 10⁻⁷** in both implementations (the port follows the reference
   here, so that the two agree to rounding). `conj_upper` puts 6000 points on a fine part around the maximisers and
   1500 on each coarse part. When the maximisers coincide, the fine part extends 0.25·span = 2.5×10⁻⁴(R − t₀)
   beyond them, which is less than half a coarse step, (R − t₀)/3000. On the first coarse subinterval the curvature
   allowance w²/8·sup W″ then exceeds what the function has lost. Measured: T(0) ≤ 1.9×10⁻⁷ for W = u² (exact
   value 0); the tail values −inf W are bounded by 7.8×10⁻⁷ for 4(u² − 1)² and 4.9×10⁻⁸ for (u² − 1)²/4 (exact: 0).
   A floor of 2×10⁻³(R − t₀) for the span, or a graded partition, removes it. The effect on the examples is below
   10⁻⁶ relative.
3. **Data exactly at an end of the t-range** (g = t_lo on a face, as in the ramp at x = −1) are treated as outside
   the range by both implementations, because the outward-rounded g_min falls below t_lo; the Lipschitz penalty is
   then used where the Taylor penalty would be valid. This could be decided exactly. Small.
4. `scalarnd.py` and `check_exact.py` still refer to the note as `NOTE.md`; it is `docs/SCALAR_ND.md`.
5. On another platform the optimiser of the reference takes a different path: the calibration rows of
   `check_exact.py` differ from the Linux log in the sixth digit (for instance 0.896002 against 0.896009). Both
   runs pass; the rigorous evaluators agree to rounding.

## 9. Later

* Potentials V(x, u): cell-wise constant in x works unchanged; a source f(x)u can be reduced to that with the
  Poincaré inequality on each cell (cost O(h²)); see `docs/SCALAR_ND.md`, Section 9.
* Non-polynomial potentials: the method needs W at exact sample points, bounds for W″ on intervals (symbolic
  derivative + `Interval`), and the tail bounds, which `Resolve` can prove as in `ScalarField1D.wl`.
* Exactly equilibrated fluxes for W(u) = −f u (the Poisson problem).
* Higher order in t (cubic splines) to reduce the number of knot intervals; the convexity argument of Lemma 3 then
  has to be replaced by a bound on Φ″.
* Speed of calibration mode: the objective is bound by memory bandwidth (2ⁿ passes over cells × samples). The
  vertex passes are independent and could run in parallel; LibraryLink would help as for the Skyrme programme.

## 10. ScalarFieldND 2 (branch `scalar-nd-2`): log

Scope: position-dependent linear terms, a different potential in each cell, potentials on an interval that are
not polynomials, a prescribed mean, periodic boundary conditions, potentials without superlinear growth. The
statements and proofs are in `docs/SCALAR_ND_2.md`; the statements of `docs/SCALAR_ND.md` keep their numbers.

### 23 November 2026

**Done.**
* `docs/SCALAR_ND_2.md`: Theorem 8 (Theorem 1 with a constant matrix K, W(x,u), a value range and periodic
  faces; the earlier proof checked line by line), Corollary 9 (prescribed mean), Proposition 10 (cells to vertices
  with the term a(x)u), Lemma 11 (the subinterval that ends at a singular point of an entropy term: the closed
  form max(Phi(beta), Phi(alpha) + c w)), Lemma 12 (logarithm from the four operations), Proposition 13 (the
  competitor: linear term, entropy terms, tables, exact mean by a shift), Lemmas 14 and 15 and Proposition 16
  (linear potentials through the Poincare inequality on a box), the six known answers for the tests.
* Decisions. (1) A companion note, not new sections of the first one, so that no number changes. (2) For
  potentials without growth, the residual route (Proposition 16), not exact equilibrium: reasons in Section 5 of
  the note. (3) Potentials on an interval and tables are evaluated in calibration mode only, with the spline
  range equal to the value range, so they have no tails. (4) A matrix K is covered by the theorem and not
  implemented. (5) "A different potential in each cell" is two things: a finite list of materials with an integer
  per cell, and a table per cell that is piecewise linear in u.

* `ia.py`: `ln` (Lemma 12), `xlnx` (s ln s on intervals that may start at 0), `sqrt_up` (verified square root),
  `pi_lower` (Machin). New functions only; nothing existing changed. `check_ia.py` compares each with exact
  rational arithmetic (607 + 303 + 3000 points): all enclosures hold, relative width of ln 4e-15.

* `scalarnd.py`, first part of step 2. Potentials are objects (`Poly`, `PolyLog` with the entropy terms and the
  rule of Lemma 11 in `sub_bound`); `Problem` takes `lin` (nodal or cell values of a), `mean`, `mat` (materials)
  and `('P', 'P')` in `bc`, with `tie`/`tied` for periodic arrays and methods `Wenc`, `sub_bound`, `Wenc_cells`,
  `W2rng_cells` that return scalars for one polynomial and arrays over the cells otherwise. `upper_bound` works on
  nodal intervals (`nodal_enclosure`: the shift that makes the mean exact), adds the term a(x)u exactly, and for
  potentials that are not plain polynomials takes the better of the midpoint rule and the vertex rule
  (Proposition 13). `check_unchanged.py` (new) recomputes the five test vectors of version 0.3.0 and compares
  bit for bit: identical. `verify_all.py`: all stored certificates verified.
  **Not yet exercised by any test: every new branch of `upper_bound`** (linear term, mean, `PolyLog`, materials).
  They are written from the propositions and have not been run.

* `scalarnd.py`, second part of step 2 (rigorous lower bound). `Calibration` has `lam` (the multiplier of the
  mean), `opt_lam`, `tie` and `respects_periodic_faces`. `lower_bound`: the per-vertex term (lam − a_v) t, the
  shifted tails (`_tails_shifted`), the end rule of Lemma 11 through `Problem.sub_bound`, no tails when the
  potential has a value range, + lam m |Ω|. `poisson_lower`: Proposition 16, both cases.
* `check_theorem2.py` (new, seconds): known answers with fields written down by hand, so only the evaluators
  run. Output of the last run: source f = 4k, lower bound −0.16666666666666824 against −1/6 in 2D, in 2D with
  periodic sides and a on the nodes, and in 3D (exact to rounding, as required); source f = c x enclosed at three
  grids; convex potential with prescribed mean, |Ω|W(m) enclosed with the lower bound 2×10⁻⁷ below it (the slack
  of the conjugate tables known from version 0.3.0); source with prescribed mean and no Dirichlet face, bounds do
  not cross; competitor energy with entropy terms, enclosure 10⁻¹⁴ wide around the floating-point value.
  `check_unchanged.py` after these edits: identical bit for bit.

### 23 November 2026, state at the stop

**Tested:** what the two scripts above cover, and `check_ia.py` as run before the last edit of `sqrt_up`.

**Written but never executed** (treat as untested code):
* `lower_bound` with K > 0 and any new feature: the term (lam − a_v) t at the samples, `_tails_shifted`, the end
  rule (the branch with `ca`, `cb`), `PolyLog` in `sub_bound`/`enc`, potentials on a value range (no tails).
* `upper_bound`: materials, the vertex rule, a candidate whose mean has to be shifted (`nodal_enclosure` was
  only run with candidates that already had the right mean).
* `Problem.spec` / `from_spec` with the new keys; `Calibration.pack`/`unpack` with `opt_lam`.

**Known gaps and caveats.**
* `sqrt_up` returned after a very long loop for subnormal arguments (it hung the second case of test 1); it now
  returns 1e-150 for arguments below 1e-300. `check_ia.py` has not been re-run since.
* In test 2 the lower bound converges at first order only. The cells there are h × 1/2 and Proposition 16 uses
  h_max. Use cells of equal sides in tests, or sharpen Lemma 14 for a coefficient that varies in one direction.
* `refine_t`, `prolong_x`, `save_certificate`, `load_certificate`, `export_field_dir` do not carry `lam`, `lin`,
  `mat`, `mean` yet; `certify_convex` and `certify_calibration` know nothing of the new features.
* The heuristics are untouched: `discrete_energy`, `solve_candidate`, `fv_energy`, `fv_solve`,
  `affine_from_candidate`, `optimise_calibration` ignore the linear term, the mean, periodic faces and materials
  (for a problem that is not plain they would compute with the wrong energy; nothing rigorous depends on them).
* The table kind of potential (a table per cell) is not written at all.
* `docs/SCALAR_ND_2.md`: no lemma is half-proved. References [1] and [2] there are from memory and unverified;
  the literature search has not been done.
* Not started: `check_exact2.py`, falsification, stored certificates, test vectors, the Wolfram Language port
  (`Kernel/ScalarFieldND.wl` and `Tests/` are as in version 0.3.0), the three examples, README,
  `docs/INDUSTRY.md`.

**In progress.** Nothing is half-edited: every file is in a state that imports and passes the scripts above.

**Next, in this order.**
1. (done, see above)
2. `scalarnd.py`, the rigorous evaluators first. Plan that keeps every existing result bit for bit: potentials
   become objects (`Poly`, `PolyLog`, `Table`, a list with a material map) behind methods of `Problem` that return
   a scalar for one polynomial, exactly as now, and arrays over the cells otherwise; `Problem` gains `lin`
   (nodal or cell values of a), `mean`, and `'P'` in `bc`; periodic arrays keep their full shape with the last
   slice tied to the first; `Calibration` gains `lam`. `upper_bound` takes nodal intervals (Proposition 13);
   `lower_bound` gets the per-vertex term, the shifted tails (stack the 2^n shifted arguments into one call of
   `conj_cells`), the end rule of Lemma 11, and the branch for linear potentials (Proposition 16).
3. `check_exact2.py`: the same formulas in exact rational arithmetic with no shared code (logarithms by rational
   series with remainder, square roots by integer square roots), on small grids, for random fields.
4. The heuristics: finite-volume solve with a(x), the mean constraint and periodic faces (direct solve for linear
   potentials); the optimiser with the per-vertex term, the multiplier and tied faces; the candidate.
5. `check_theorem2.py`: the known answers of Section 7 of the note, then falsification as in `check_theorem.py`.
6. Stored certificates (format 2 for problems that use a new feature; format 1 untouched), test vectors.
7. The Wolfram Language port and `Tests/ScalarFieldND2.wlt` (flag `--nd2`); the three examples; README;
   `docs/INDUSTRY.md`; the literature search that the note lists as open.

**Checks that must keep passing after every step:** `verify_all.py`, `check_exact.py`, `check_theorem.py`, and
`tools/run_tests.wls --nd` with `Tests/ScalarFieldND.wlt` byte-identical.
