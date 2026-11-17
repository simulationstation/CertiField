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
In[1]:= PacletInstall["https://github.com/simulationstation/CertiField/releases/download/v0.2.0/CertiField-0.2.0.paclet"];
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
  The hedgehog energy uses Simpson enclosures with a fourth-derivative remainder in interval arithmetic.
* **Lower bounds** come from a model-specific mathematical argument, evaluated with certified arithmetic:
  * one dimension: the calibration k u′² + V + c ≥ 2√(k(V+c))|u′|, with V + c ≥ 0 on ℝ *proved* by quantifier
    elimination, and the integral of √(V+c) enclosed by adaptive interval quadrature;
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
PacletInstall["https://github.com/simulationstation/CertiField/releases/download/v0.2.0/CertiField-0.2.0.paclet"]
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
| `$CertiFieldModels` | `ScalarField1D`, `Phi4Kink`, `SineGordonKink`, `MagneticDomainWall`, `SkyrmeDegreeOne` |

Options: `"Resolution"` (`"Quick"`, `"Standard"`, `"Publication"` for Skyrme), `"GapTolerance"`, `"Nodes"`,
`"Verbose"`, `WorkingPrecision`.

**Timings** (one core, Wolfram Engine 15.0): 1D models take 0.5–10 s. Skyrme at Quick resolution takes 25–35 s per
pion mass, including the minorant table, which is certified once (about 10 s), cached on disk and re-certified on
load. Standard takes about 80 s, and Publication about an hour (59 min at μ = 1).

## Validation

* The Wolfram Language dynamic programme reproduces the independent Python verifier on identical inputs. At μ = 1 it
  checks the same 5,025,904 transitions with the same 84,771 fallbacks, and the values agree to 3×10⁻¹². Transition
  counts match on every preset checked, including Publication at μ = 1 (517,125,854 transitions, bound within
  1.6×10⁻¹¹).
* `VerifyCertificate` re-derives the μ = 1 Skyrme certificate from its recorded inputs in 26 s.
* The exact minima of the φ⁴, sine-Gordon and domain-wall models lie inside the certified enclosures, which are
  10⁻¹³ wide.
* Tampering with a stored bound makes `VerifyCertificate` return `False`.
* Non-smooth potentials, boundary values that are not vacua, unprovable positivity and unsupported parameters return
  `Failure` objects, never numbers.
* The test suite is in `Tests/` (`wolframscript -file tools/run_tests.wls`; add `--skyrme` for the Skyrme suite).

## Roadmap

* **Certified error bounds for finite-element solutions** of convex elliptic problems (elasticity, diffusion) via
  complementary energy and equilibrated fluxes, applied to `NDSolve`'s own FEM output.
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
  CertifiedBound.wl       the result object, safe decimal display, reports
  Models.wl               model registry, public entry points, VerifyCertificate
Data/Skyrme/     grids for 4 pion masses x 3 resolutions (heuristic inputs: they affect tightness, never validity)
Tests/           VerificationTest suites
tools/           test runner and syntax checks
docs/            architecture and method notes
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

## License

CertiField is source-available under the [PolyForm Noncommercial License 1.0.0](LICENSE). It is free for research,
teaching, personal projects, and use by educational, nonprofit and government institutions. **Commercial use requires
a separate commercial license** from SimulationStation; to arrange one, open an issue titled "Commercial license".

## About

CertiField is developed by SimulationStation. For collaboration, commercial licensing or certified-bound projects on
your own models, open an issue.
