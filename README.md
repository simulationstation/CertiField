# CertiField

**Certified optimality gaps for variational problems, in Wolfram Language.**

Every simulation that minimises an energy tells you what it found. None of them tells you how far the true optimum
could be below it. CertiField closes that gap with proof: give it a model and a candidate solution from any solver,
and it returns

* a **certified upper bound**: the energy of an explicitly admissible competitor built from your candidate, computed in
  exact rational or outward-rounded interval arithmetic, and
* a **certified lower bound that holds for every admissible field**: not just for the fields on your mesh or in your
  ansatz, but over the whole infinite-dimensional admissible class,

packaged as a `CertifiedBound` object that anyone can re-check in their own kernel with `VerifyCertificate`.

```wolfram
In[1]:= PacletInstall["https://github.com/simulationstation/CertiField/releases/download/v0.3.0/CertiField-0.3.0.paclet"];
        Needs["CertiField`"]

In[2]:= cert = CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> 1];     (* about 30 s *)

In[3]:= cert["Summary"]
Out[3]= "SkyrmeDegreeOne: Certified  1.41402 <= inf E <= 1.41623  (relative gap <= 0.156%)"

In[4]:= cert["Theorem"]
Out[4]= Every finite-energy degree-one map U: R^3 -> S^3 with pion mass mu = 1 has E(U)/(12 Pi^2) >= 1.4140234.
        The hedgehog trial map has E/(12 Pi^2) <= 1.4162244.

In[5]:= VerifyCertificate[cert]
Out[5]= True
```

When a bound cannot be established, CertiField says so (`"Inconclusive"`, or a `Failure` explaining why) instead of
reporting a number it cannot prove.

Three kinds of problem are supported: the Skyrme model (a computer-assisted theorem about all degree-one maps),
scalar field theories in one dimension (kinks, domain walls), and **your own energy ∫ k|∇u|² + W(u) on a box in two
or three dimensions**, convex or not, with a candidate from any solver
([below](#bring-your-own-energy-scalar-fields-in-two-and-three-dimensions)).

---

## Why this matters

Numerical minimisation produces **upper bounds** for free: any trial state's energy is at least the minimum. What it
never produces is the other half. Mesh refinement, convergence plots and ensembles of random starts are evidence, not
proof, and they say nothing about configurations your discretisation or ansatz cannot represent.

A certified **lower bound** over the full admissible class turns "we found a good solution" into a theorem:

> *no admissible configuration whatsoever has energy more than X% below this one.*

That is what you need when a simulated optimum is used as a reference value, a benchmark, a stability argument, a
design margin, or the input to a proof.

## Headline result: the Skyrme model

The Skyrme model describes baryons as topological solitons. Whether the spherically symmetric "hedgehog" is the
global energy minimiser in the degree-one (single-baryon) sector is a long-standing open question: extensive numerics
support it, but it is unproved. The classical Faddeev–Bogomolny bound E ≥ 12π²|B| is far from tight; in the massless
model the hedgehog's energy is 23% above it.

CertiField implements the **level-set lower bound** of [[1]](#references) (coarea formula on the level sets of the
distance to the vacuum, isoperimetry, and a convexified level energy, reduced to a dynamic programme and certified
with outward-rounded arithmetic), extended to the model with a pion mass [[2]](#references). The result is a
computer-assisted theorem about **all** finite-energy degree-one maps, with **no symmetry assumption**:

| Pion mass μ | Certified lower bound, E/12π² (all degree-one maps) | Hedgehog trial map (certified upper bound) | Certified gap |
|---|---|---|---|
| 1 | **1.414324** | 1.416225 | **≤ 0.135%** |
| 0.65 | 1.330923 | 1.335151 | ≤ 0.317% |
| 0.526 | 1.301718 | 1.307666 | ≤ 0.455% |
| 0 (massless) | 1.192040 | 1.231447 | ≤ 3.20% |

These lower bounds use the Publication-resolution grids. CertiField reproduces the μ = 1 value in Wolfram Language in
about an hour: 1.4143249628948 against the reference verifier's 1.4143249629108, over the same 517,125,854 transitions.
Upper bounds are computed by CertiField. All entries are rounded in the safe direction. Battye–Sutcliffe units:
E = ∫ |dU|² + |Λ²dU|² + 2μ²(1 − cos f).

In plain terms: **at μ = 1, no degree-one Skyrme field of any shape can have less than 99.865% of the hedgehog's
energy.** The hedgehog is either the minimiser or within 0.135% of it.

Running the certificate end to end in CertiField (one core, Wolfram Engine 15.0):

| μ | Resolution | Certified bounds, E/12π² | Gap | Time |
|---|---|---|---|---|
| 1 | Quick | 1.41402 ≤ inf E ≤ 1.41623 | 0.156% | 36 s (includes certifying the minorant table once) |
| 1 | Standard | 1.41420 ≤ inf E ≤ 1.41623 | 0.143% | 77 s |
| 1 | Publication | 1.41432 ≤ inf E ≤ 1.41623 | 0.135% | 59 min |
| 0.65 | Quick | 1.33062 ≤ inf E ≤ 1.33516 | 0.340% | 27 s |
| 0.526 | Quick | 1.30142 ≤ inf E ≤ 1.30767 | 0.478% | 28 s |
| 0 | Quick | 1.19174 ≤ inf E ≤ 1.23145 | 3.23% | 33 s |

**Other pion masses.** `"PionMass"` accepts any rational μ with 0 ≤ μ < 2. For a mass without shipped data the grids of
the lower bound and the hedgehog trial map are computed on request (about two minutes at Quick resolution), and the
certificate stores them. `Method -> "Monotonicity"` instead uses the lower bound of the largest shipped mass below μ,
which holds because the energy increases with the pion mass. Details: `docs/SKYRME_ANY_MASS.md`.

**What is and is not claimed.** The lower bound is a theorem (analytic reduction plus certified computation). It does
*not* prove that the hedgehog is the minimiser. It bounds how much better any competitor could possibly be. Every
`CertifiedBound` carries a `"NotClaimed"` field saying exactly this.

### Application: degree-two Skyrmions exist

Esteban's criterion says the degree-two infimum I₂ is attained if I₂ < 2I₁. The preprint [[1]](#references) proves
this with computer assistance in the massless model. An explicit axially symmetric map gives I₂ ≤ 2.358835·12π²,
and the level-set bound gives I₁ ≥ 1.189162·12π², leaving a margin of 0.0194·12π².

CertiField re-derives the lower-bound half independently, in a different language and with different arithmetic.
The μ = 0 Quick certificate takes 33 seconds and gives **I₁ ≥ 1.191748·12π²**, sharper than the bound used in the
preprint. With the preprint's upper bound, the margin for I₂ < 2I₁ grows to more than 0.0246·12π².

```wolfram
In[6]:= c0 = CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> 0];
        {c0["LowerBound"], 2 c0["LowerBound"] > 2.358835}
Out[6]= {1.1917480092669916, True}
```

### Reproduced independently, in your kernel

The Skyrme certificate was first produced by a Python/NumPy verifier. CertiField re-implements it from scratch in
Wolfram Language: exact rational enclosures for every τ-dependent constant, outward-rounded packed-array arithmetic,
and a minorant table that is rebuilt and **certified entry by entry inside the kernel** (1.8 million entries, about 10
s). On identical inputs the two implementations agree structurally and numerically:

| μ = 1, Quick resolution | Python reference | CertiField (Wolfram Language) |
|---|---|---|
| Certified lower bound | 1.414023462272984 | 1.4140234622695356 |
| Transitions checked | 5,025,904 | 5,025,904 |
| Second-order fallbacks | 84,771 | 84,771 |
| Far-target bindings | 0 | 0 |

The 3×10⁻¹² difference comes from CertiField's slightly more conservative rounding. On every preset checked (all four Quick grids, and
Standard and Publication at μ = 1) the transition counts are identical. The fallback counts are identical for every
Quick and Standard preset with μ > 0. At Publication they differ by 6 out of 576,134, because CertiField rebuilds the
minorant table itself and borderline steps can fall either way. At μ = 0 there is no mass term to give slack, so
many step conditions sit exactly at the edge. There, CertiField's wider rounding margin triggers 649 extra
first-order fallbacks (0.013% of transitions), and its bound is 1.2×10⁻⁷ lower. Being more conservative is always
safe.

## Bring your own candidate: one-dimensional field theories

For scalar field theories E[u] = ∫ k u′² + V(u) dx (kinks, sine-Gordon solitons, magnetic domain walls, phase-field
interfaces), the candidate can be anything: a formula, an `NDSolve` solution, or a list of data points. The lower bound
comes from a calibration (Bogomolny-type) argument that is exact whenever the minimiser is monotone, so a good
candidate gives a near-zero certified gap.

```wolfram
In[7]:= CertifiedVariationalBound["SineGordonKink"]["Summary"]
Out[7]= "SineGordonKink: Certified  7.999999999999956 <= inf E <= 8.000000000000212  (relative gap <= 3.19x10^-12 %)"

In[8]:= sol = NDSolveValue[{u''[x] == u[x]^3 - u[x], u[0] == -1, u[3] == 1/2}, u, {x, 0, 3}];
        CertifiedVariationalBound[
          CertiFieldModel["ScalarField1D", <|"Potential" -> Function[s, (s^2 - 1)^2/4],
             "Domain" -> {0, 3}, "BoundaryValues" -> {-1, 1/2}|>], sol]["Summary"]
Out[8]= "ScalarField1D: Certified  0.80350006723246 <= inf E <= 0.80350006723404  (relative gap <= 1.97x10^-10 %)"
```

| Model | Exact minimum | Certified enclosure | Relative gap | Time |
|---|---|---|---|---|
| φ⁴ kink | 2√2/3 | 0.9428090415820432 ≤ E ≤ 0.9428090415820903 | 5×10⁻¹⁴ | 0.5 s |
| Sine-Gordon kink | 8 | 7.999999999999956 ≤ E ≤ 8.000000000000212 | 3×10⁻¹⁴ | 10 s |
| Bloch domain wall (A = 2, K = ½) | 4√(AK) = 4 | 3.999999999999922 ≤ E ≤ 4.000000000000106 | 5×10⁻¹⁴ | 6 s |
| Double well, finite domain, `NDSolve` candidate | (unknown) | 0.80350006723246 ≤ E ≤ 0.80350006723404 | 2×10⁻¹² | 1.3 s |
| φ⁴ kink, deliberately poor candidate tanh x | 2√2/3 | 0.942809 ≤ E ≤ 1.000000 | 5.7% | 0.6 s |

The last row is the point: a poor candidate gets an honest, wide gap. With `"GapTolerance" -> 10^-3` the same call
returns status `"Inconclusive"`.

## Bring your own energy: scalar fields in two and three dimensions

`ScalarFieldND` certifies

    E[u] = ∫_Ω k |∇u|² + W(u) dx

on a box Ω in 2D or 3D, with u prescribed on some faces and free on the others. W is any polynomial with rational
coefficients; it does not have to be convex. The candidate can be a function of the coordinates, an `NDSolve`
solution, or an array of nodal values. Without one, the package computes its own.

```wolfram
In[9]:= model = CertiFieldModel["ScalarFieldND", <|
          "Potential" -> Function[u, u^4/4 + u^2/2 - u],
          "Domain" -> {{0, 1}, {0, 1}, {0, 1}}, "Cells" -> 32, "BoundaryData" -> 0|>];
        cert = CertifiedVariationalBound[model];          (* 3 s *)
        cert["Summary"]
Out[9]= "ScalarFieldND: Certified  -0.00980025 <= inf E <= -0.00976023  (relative gap <= 0.410%)"

In[10]:= cert["Theorem"]
Out[10]= Every u in H^1 of the box {{0, 1}, {0, 1}, {0, 1}} whose trace on the Dirichlet faces is the piecewise
         multilinear interpolant of the nodal boundary values on the 32x32x32 grid has E[u] >= -0.0098002440;
         the competitor has E <= -0.0097602333.

In[11]:= sol = NDSolveValue[{-Laplacian[u[x, y], {x, y}] + u[x, y]^3 + u[x, y] == 1,
            DirichletCondition[u[x, y] == 0, True]}, u, {x, 0, 1}, {y, 0, 1}];
         square = CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, u^4/4 + u^2/2 - u],
            "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 64, "BoundaryData" -> 0|>];
         CertifiedVariationalBound[square, sol]["Summary"]          (* your own candidate *)
Out[11]= "ScalarFieldND: Certified  -0.0167643 <= inf E <= -0.0167543  (relative gap <= 0.0591%)"

In[12]:= ramp = CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, 4 (u^2 - 1)^2],
            "Domain" -> {{-1, 1}, {-1, 1}}, "Cells" -> 32, "BoundaryData" -> Function[{x, y}, x]|>];
         CertifiedVariationalBound[ramp]["Summary"]          (* a double well, not convex: under a minute *)
Out[12]= "ScalarFieldND: Certified  8.24547 <= inf E <= 8.30526  (relative gap <= 0.720%)"
```

The lower bound holds for **every** function in H¹ with the given boundary values, not only for functions on the grid.
It is the classical calibration inequality: for any vector field S(x, t) with matching normal components,

    E[u] ≥ ∫ S(x, g)·ν over the Dirichlet faces − ∫_Ω sup_t ( div_x S + |∂_t S|²/(4k) − W(t) ) dx,

evaluated with outward-rounded arithmetic. The field S is found numerically and never trusted. When W is convex it
is the flux of a finite-volume solve (the complementary-energy bound) and the certificate takes seconds, also in
3D. When W is not convex the field is optimised coarse to fine, which takes from under a minute to about twenty
minutes; `"Verbose" -> True` shows the certified bound rising. The mode is chosen automatically. Statement and
proofs: [`docs/SCALAR_ND.md`](docs/SCALAR_ND.md).

| Potential W | Domain and data | Grid | Certified bounds | Gap | Time |
|---|---|---|---|---|---|
| u² − u | u = 0 on the boundary of the unit square | 128×128 | −0.01602601 ≤ inf E ≤ −0.01602361 | 0.015 % | 1.7 s |
| u² − u | u = 0 on the boundary of the unit cube | 48×48×48 | −0.00950668 ≤ inf E ≤ −0.00948922 | 0.184 % | 6.1 s |
| u⁴/4 + u²/2 − u | u = 0 on the boundary of the unit cube | 32×32×32 | −0.00980025 ≤ inf E ≤ −0.00976023 | 0.410 % | 2.7 s |
| (u² − 1)²/4 (double well) | wall in (−3,3)×(0,1)², free sides | 48×2×2 | 0.944198 ≤ inf E ≤ 0.945551 | 0.143 % | 54 s |
| 4(u² − 1)² (double well) | u = x on the boundary of (−1,1)² | 32×32 | 8.245472 ≤ inf E ≤ 8.305251 | 0.720 % | 43 s |
| 4(u² − 1)² (double well) | the same | 64×64 | 8.275406 ≤ inf E ≤ 8.291702 | 0.197 % | 6 min |
| (u² − 1)² (double well) | u = x on the boundary of (−1,1)³ | 16×16×16 | 8.005335 ≤ inf E ≤ 8.040144 | 0.433 % | 7 min |
| (u² − 1)² (double well) | the same | 32×32×32 | 8.015122 ≤ inf E ≤ 8.027319 | 0.152 % | 18 min |

Times are for the whole certificate (Wolfram Engine 15.0, 16-core processor): candidate, field, and both certified
evaluations. `VerifyCertificate` then recomputes both bounds in 1 to 15 seconds. For the first two rows
the exact minimum is known from a sine series (−0.01602506 and −0.00949862), and for the wall it is the
one-dimensional value 0.94496743. All three lie inside the certified intervals.

**What is and is not claimed.** The bounds are theorems about the stated problem, whose boundary datum is the
piecewise multilinear interpolant of the given boundary values on the grid. They do not say that the candidate is the
minimiser. The inequality is classical and no novelty is claimed for it; the proofs in `docs/SCALAR_ND.md` have been
checked by one author only. Scalar fields on boxes with polynomial potentials only; no source terms f(x)u yet, and
therefore not the plain Poisson problem; the method does not extend to vector-valued fields.

## Linear elasticity in two dimensions

`LinearElasticity` certifies the minimal energy of a linearly elastic body in plane strain,

    E[u] = ∫_Ω ½ ε(u) : C : ε(u) − b·u dx,

on a rectangle Ω, with the displacement prescribed on some faces and no traction on the others. The material is
isotropic (Lamé constants) or given by its 3×3 stiffness matrix, and the body force b is constant. The candidate
can be a function of the coordinates, an `NDSolve` solution, or an array of nodal values. Without one, the package
computes the finite-element solution on the grid.

```wolfram
In[13]:= block = CertiFieldModel["LinearElasticity", <|"Lame" -> {3/2, 1},
            "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 32, "Degree" -> 2,
            "TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}},
            "BoundaryData" -> Function[{x, y}, {x/10, 0}]|>];
         cert = CertifiedVariationalBound[block];          (* under a second *)
         cert["Summary"]
Out[13]= "LinearElasticity: Certified  0.0150234 <= inf E <= 0.0150298  (relative gap <= 0.0418%)"

In[14]:= VerifyCertificate[cert]
Out[14]= True
```

This block is clamped on its left face, pulled by 0.1 on its right face, and free on the other two. The lower
bound holds for **every** displacement in H¹ with these boundary values. It is the complementary-energy inequality
of Prager and Synge: for any stress field σ that is in equilibrium with the body force and has no traction on the
free faces,

    E[u] ≥ ∫ (σν)·g over the Dirichlet faces − ½ ∫_Ω σ : C⁻¹ : σ dx.

Equilibrium must hold exactly, not to rounding. The stress field is therefore derived from an Airy stress function
in a C¹ finite-element space, which satisfies equilibrium identically for any nodal values, and its nodal values
on the free faces are fixed in exact rational arithmetic. The right-hand side is then a polynomial in machine
numbers with rational coefficients, and it is evaluated with outward-rounded arithmetic. The stress field comes
from one sparse linear solve and is never trusted. The gap also bounds the energy-norm error of the candidate and
of the stress field. Statement and proofs: [`docs/ELASTICITY.md`](docs/ELASTICITY.md).

| Problem | Grid, degree | Certified bounds | Gap | Time |
|---|---|---|---|---|
| block: clamped on the left, pulled on the right, free sides | 16×16, 1 | 0.01501720 ≤ inf E ≤ 0.01506422 | 0.313 % | 0.3 s |
| the same | 32×32, 2 | 0.01502346 ≤ inf E ≤ 0.01502974 | 0.042 % | 0.6 s |
| the same | 64×64, 2 | 0.01502578 ≤ inf E ≤ 0.01502813 | 0.016 % | 1.8 s |
| column under its own weight, clamped at the bottom, three free faces | 16×32, 2 | −0.45145662 ≤ inf E ≤ −0.45112300 | 0.074 % | 0.5 s |
| bar in pure bending, free sides | 4×2, 2 | 0.0312499999 ≤ inf E ≤ 0.0312500001 | 2×10⁻⁷ % | 0.2 s |

`"Degree"` is the polynomial degree of the candidate in each variable (1 or 2). The gap closes slowly on the block
because the corners where a clamped face meets a free face carry stress singularities. For the bent bar the exact
minimum 1/32 is known, and it lies inside the certified interval.

**What is and is not claimed.** The bounds are theorems about the stated problem, whose boundary datum is the
piecewise polynomial of degree `"Degree"` through the given boundary values on the grid. They do not say that the
candidate is the minimiser. The inequality and the use of stress functions are classical ("dual analysis" with
equilibrium finite elements), and no novelty is claimed for them; what the package adds is rigorous rounding and a
certificate that can be re-checked. The proofs in `docs/ELASTICITY.md` have been checked by their author only. Rectangles in two dimensions, constant stiffness and
body force only; no prescribed tractions, no roller conditions, not yet three dimensions.

## How it works

```
 your candidate ──► admissible competitor ──► certified energy ──────────► upper bound
 (formula, NDSolve,  (exact rational cubic     (exact rational arithmetic
  data, trial map)    Hermite, boundary         or interval quadrature)
                      values imposed exactly)
                                                                              CertifiedBound
 the model ──────► analytic lower-bound ─────► certified computation ─────► lower bound    + theorem text
                    argument (calibration,      (outward rounding, exact       (all fields)  + VerifyCertificate
                    coarea + isoperimetry)       constants, proved positivity)
```

* **Upper bounds** are valid however poor the candidate; only their tightness depends on it. CertiField converts the
  candidate into an admissible competitor and computes that competitor's energy rigorously. Polynomial potentials are
  integrated exactly in rational arithmetic; other potentials use a Taylor model of order 8 with a Lagrange remainder.
  The hedgehog energy uses Simpson enclosures with a fourth-derivative remainder in interval arithmetic. In two and
  three dimensions the competitor is the piecewise multilinear function through the candidate's nodal values: exact
  gradient energy, midpoint rule with remainder for the potential. In linear elasticity it is the piecewise
  polynomial displacement through the candidate's nodal values, whose energy is an exact quadratic form.
* **Lower bounds** come from a model-specific mathematical argument, evaluated with certified arithmetic:
  * one dimension: the calibration k u′² + V + c ≥ 2√(k(V+c))|u′|, with V + c ≥ 0 on ℝ *proved* by quantifier
    elimination, and the integral of √(V+c) enclosed by adaptive interval quadrature;
  * scalar fields in two and three dimensions: a calibration field S(x, u) (Raviart–Thomas in x, B-splines in u)
    whose coefficients come from a finite-volume solve or an optimiser and are never trusted; the bound for those
    coefficients is evaluated with exact sample points and outward-rounded arrays, cell by cell and face by face;
  * linear elasticity: the complementary-energy bound with a stress field from an Airy stress function, which is in
    equilibrium identically and traction-free exactly; the bound is a polynomial with rational coefficients in the
    stored machine numbers, evaluated cell by cell with outward rounding;
  * Skyrme: the level-set dynamic programme with τ-scalars as exact rationals (Taylor series with remainder, 35-digit
    π), every array operation outward-rounded, square and cube roots verified by squaring/cubing, and every entry of the
    affine-minorant table certified by a cubic positivity test.
* **Trusted base** (also stored in each certificate): IEEE-754 arithmetic for + − × ÷, Wolfram Language exact and
  `Interval` arithmetic, `Resolve` for polynomial and trigonometric positivity, and the analytic reduction
  (the lower-bound theorem) for each model.

### Why Wolfram Language

The certificate needs exact rationals, interval arithmetic, symbolic differentiation, quantifier elimination and fast
compiled array code, all in one system that the person checking the result already has. That combination is unusual,
and it is what lets a third party re-verify a computer-assisted theorem with one function call.

## Install and use

Requires Wolfram Language 13 or later (developed and tested on 15.0, including the free Wolfram Engine).

```wolfram
PacletInstall["https://github.com/simulationstation/CertiField/releases/download/v0.3.0/CertiField-0.3.0.paclet"]
(* or, from a clone of this repository: PacletDirectoryLoad["path/to/CertiField"] *)
Needs["CertiField`"]

cert = CertifiedVariationalBound["Phi4Kink"]          (* built-in reference candidate *)
cert = CertifiedVariationalBound["Phi4Kink", Function[x, Tanh[x/Sqrt[2]]]]    (* your own candidate *)
cert["Report"]                                         (* theorem, admissible class, methods, trusted base *)
cert[{"LowerBound", "UpperBound", "RelativeGap"}]
VerifyCertificate[cert]
```

| Function | Purpose |
|---|---|
| `CertifiedVariationalBound[model, candidate]` | certified two-sided bounds; returns a `CertifiedBound` |
| `CertiFieldModel[name, params]` | a model specification (`"ScalarField1D"` takes any smooth potential) |
| `CertifiedBound[...]["prop"]` | `"LowerBound"`, `"UpperBound"`, `"RelativeGap"`, `"Status"`, `"Theorem"`, `"AdmissibleClass"`, `"TrustedBase"`, `"NotClaimed"`, `"Report"`, ... |
| `VerifyCertificate[cert]` | recomputes every rigorous step from the recorded inputs |
| `CertiFieldModelInfo[name]` | energy, parameters, admissible class and methods of a model |
| `$CertiFieldModels` | `ScalarField1D`, `Phi4Kink`, `SineGordonKink`, `MagneticDomainWall`, `ScalarFieldND`, `LinearElasticity`, `SkyrmeDegreeOne` |

Options: `"Resolution"` (`"Quick"`, `"Standard"`, `"Publication"` for Skyrme), `"GapTolerance"`, `"Nodes"`,
`"Verbose"`, `WorkingPrecision`; `Method` (`"Monotonicity"`, `"GeneratedGrid"`) for `SkyrmeDegreeOne`; for `ScalarFieldND` also `Method` (`Automatic`, `"Convex"`, `"Calibration"`) and
`"RefineCandidate"` (use the given candidate as the starting point of the package's own minimisation).
`LinearElasticity` takes `"RefineCandidate"` (replace the given candidate by the finite-element solution).
`CertiFieldModelInfo[name]` lists the parameters of each model.

**Timings** (one core, Wolfram Engine 15.0): 1D models take 0.5–10 s. Skyrme at Quick resolution takes 25–35 s per
pion mass, including the minorant table, which is certified once (about 10 s), cached on disk and re-certified on
load. Standard takes about 80 s, and Publication about an hour (59 min at μ = 1). `ScalarFieldND` (several threads)
takes 1 to 7 s for convex energies up to 48³ cells and 39 s to 18 min for the non-convex examples above.
`LinearElasticity` takes 0.2 to 2 s for grids up to 64×64.

## Validation

* The Wolfram Language dynamic programme reproduces the independent Python verifier on identical inputs. At μ = 1 it
  checks the same 5,025,904 transitions with the same 84,771 fallbacks, and the values agree to 3×10⁻¹². Transition
  counts match on every preset checked, including Publication at μ = 1 (517,125,854 transitions, bound within
  1.6×10⁻¹¹).
* `VerifyCertificate` re-derives the μ = 1 Skyrme certificate from its recorded inputs in 26 s.
* The exact minima of the φ⁴, sine-Gordon and domain-wall models lie inside the certified enclosures, which are
  10⁻¹³ wide.
* In two and three dimensions the minima known from a sine series (unit square and unit cube) and from the
  one-dimensional first integral (domain wall) lie inside the certified intervals.
* `ScalarFieldND` reproduces its Python reference (`Reference/ScalarND`) on identical inputs, to 10⁻¹³ for the
  competitor energies and the lower bounds, and recomputes all eleven certificates stored by it. Its boundary-flux
  bound is tested from both sides against the exact flux in rational arithmetic. The reference passes an exact
  rational cross-check and falsification tests against known minima.
* `LinearElasticity` reproduces its Python reference (`Reference/Elasticity`) on identical inputs to 10⁻¹¹ and
  recomputes the six certificates stored by it. An independent evaluation with symbolic polynomials in rational
  arithmetic confirms that the stress fields are in equilibrium exactly, that their tractions are continuous and
  vanish on the free faces exactly, and that the exact values of both bounds lie inside the certified enclosures.
  Minima known as rational numbers (homogeneous strain, uniaxial tension, a quadratic displacement with body force,
  pure bending) lie inside the certified intervals, and random admissible stress fields and displacements never
  give a lower bound above an energy.
* Tampering with a stored bound, a stored field or stored boundary values makes `VerifyCertificate` return `False`.
* Non-smooth potentials, boundary values that are not vacua, unprovable positivity and unsupported parameters return
  `Failure` objects, never numbers.
* The test suite is in `Tests/` (`wolframscript -file tools/run_tests.wls`; add `--skyrme` for the Skyrme suite,
  `--nd` for the two- and three-dimensional suite, `--nd-slow` for its long runs, `--elasticity` for linear
  elasticity).

## Roadmap

* **Scalar fields in 2D and 3D**: source terms f(x)u and the Poisson problem (exactly equilibrated fluxes),
  potentials that are not polynomials, domains other than boxes and `NDSolve`'s own finite-element meshes.
* **Linear elasticity**: three dimensions (Beltrami stress functions), prescribed tractions and roller conditions,
  stiffness and body forces that vary in space.
* **Chiral magnets**: two-dimensional magnetic skyrmions and domain walls with topological lower bounds.
* Higher-charge Skyrmions, and LibraryLink acceleration of the dynamic programme.

## Repository layout

```
Kernel/          the paclet (CertiField.wl loads the modules below)
  RigorousArithmetic.wl   exact enclosures, outward-rounded array arithmetic, verified roots
  RigorousQuadrature.wl   adaptive interval quadrature with smoothness guards
  Minorants.wl            certified affine minorants of h_kappa (the Skyrme table)
  SkyrmeLevelSet.wl       the level-set dynamic programme (lower bound)
  SkyrmeHedgehog.wl       hedgehog trial maps and their certified energy (upper bound)
  ScalarField1D.wl        calibration lower bound, admissible competitors, exact/Taylor-model energies
  ScalarFieldND.wl        2D/3D scalar fields: calibration-field lower bound, Q1 competitor, convex and calibration modes
  Elasticity.wl           2D linear elasticity: complementary-energy lower bound from an Airy stress function, Q_p competitor
  CertifiedBound.wl       the result object, safe decimal display, reports
  Models.wl               model registry, public entry points, VerifyCertificate
Data/Skyrme/     grids for 4 pion masses x 3 resolutions (heuristic inputs: they affect tightness, never validity)
Reference/ScalarND/  the Python reference of ScalarFieldND: implementation, checks, test vectors, stored certificates
Reference/Elasticity/  the Python reference of LinearElasticity: implementation, checks, test vectors, stored certificates
Tests/           VerificationTest suites
tools/           test runner and syntax checks
docs/            architecture (ARCHITECTURE.md), the statements and proofs behind ScalarFieldND (SCALAR_ND.md) and
                 LinearElasticity (ELASTICITY.md)
```

## References

1. *The inequality I₂ < 2I₁ for the Skyrme problem and the existence of a minimiser of degree two*, preprint (2026),
   Zenodo. [doi:10.5281/zenodo.24294210](https://doi.org/10.5281/zenodo.24294210). This preprint introduces the
   level-set inequality and the dynamic-programming certificate that CertiField implements.
2. *A level-set lower bound for the Skyrme energy*, manuscript in preparation. It extends [1] to the model with a pion
   mass. The Publication grids in `Data/Skyrme` reproduce its certified values; Quick and Standard are coarser presets.
3. M. J. Esteban, *A direct variational approach to Skyrme's model for meson fields*, Comm. Math. Phys. 105 (1986),
   571–591. Erratum, *Existence of 3D Skyrmions. Complete version*, arXiv:math-ph/0401042 (2004).
4. R. A. Battye and P. M. Sutcliffe, Skyrmions with massive pions, Phys. Rev. C 73 (2006). Source of the units and
   pion-mass conventions.

References for the two- and three-dimensional bounds (calibrations, duality, complementary energy) are listed in
`docs/SCALAR_ND.md`, those for linear elasticity (Prager and Synge, dual analysis, equilibrium finite elements) in
`docs/ELASTICITY.md`.

## License

CertiField is source-available under the [PolyForm Noncommercial License 1.0.0](LICENSE). It is free for research,
teaching, personal projects, and use by educational, nonprofit and government institutions. **Commercial use requires
a separate commercial license** from SimulationStation; to arrange one, open an issue titled "Commercial license".

## About

CertiField is developed by SimulationStation. For collaboration, commercial licensing or certified-bound projects on
your own models, open an issue.
