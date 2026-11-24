# Porting notes: `elasticity.py` → `Kernel/Elasticity.wl`

## 1. State

**Done.** `Kernel/Elasticity.wl` (model `LinearElasticity`, two dimensions) with the specification, the exact local
matrices, the exact nodal values on the traction-free faces, both certified bounds, the two heuristic solves, import
of fields and certificates written by the reference, the certificate and its re-verification; registration in
`Kernel/Models.wl`; `Tests/Elasticity.wlt` (90 tests, `tools/run_tests.wls --elasticity`).

The reference evaluates on differences as well (Section 6); its certificates and test vectors were regenerated.

**In progress.** Documentation: the results table of the Wolfram Language model in `docs/ELASTICITY.md`,
`docs/ARCHITECTURE.md`, `README.md`, `Reference/Elasticity/README.md`, the example, the build smoke test and the
CI step.

**Next.** Final runs of the existing suites (fast, `--nd`, `--skyrme`).

**Not part of this port.** Three dimensions: `linearElasticitySpec` refuses a three-dimensional domain with a clear
message.

## 2. Function map

| Python (`elasticity.py`) | Wolfram Language (`Kernel/Elasticity.wl`) | Rigorous? |
|---|---|---|
| `Problem` | `linearElasticitySpec` | exact data |
| `imat`, `IH`, `JH`, `lagrange`, `HERM` | `elIMat`, `elLagrange`, `$elHermite` | exact |
| `Problem.stiffness_local` | `elStiffnessLocal` | exact |
| `Problem.airy_local` | `elAiryLocal` | exact |
| `Problem.face_local` | `elFaceLocal` | exact |
| `Problem.traces`, `nullspace_and_particular` | `elTraces` (`RowReduce`, `MatrixRank`) | exact |
| `nodal_enclosure` | `elNodalEnclosure` | yes |
| `_quadratic` | `elQuadratic` | yes |
| `upper_bound` | `elUpperBound` | yes |
| `lower_bound` | `elLowerBound` | yes |
| `solve_candidate`, `_stiffness_global` | `elSolveCandidate`, `elGlobal1D` | no |
| `solve_stress` | `elSolveStress`, `elGlobal1D`, `elToKron`, `elFromKron` | no |
| `Stress` | `<|"X" -> ..., "Parameters" -> ...|>`, `elPackField` | data |
| `certify`, `claim_strings` | `certifyLinearElasticity` | — |
| `export_dir` | `elImportField`, `elImportArray` | data |
| `save_certificate`, `load_certificate` | `elImportCertificate` | data |
| `verify_certificate` | `elVerifyReferenceCertificate`, `verifyLinearElasticity` | yes |

## 3. Conventions

* Displacement nodal values: array of dimensions {2, M₁, M₂}, Mᵢ = p Nᵢ + 1, component first, then the nodes of the
  Q_p grid with the last index fastest (NumPy C order). Raw files are little-endian float64 in that order.
* Airy nodal values: array of dimensions {4, N₁+1, N₂+1}; the first index is 1 + d₁ + 2d₂ for the derivative
  ∂₁^{d₁}∂₂^{d₂}φ (φ, ∂₁φ, ∂₂φ, ∂₁₂φ), physical derivatives. The entries at the nodes of the traction-free faces are
  ignored: there the nodal values are ψ_b + Σ P_k ψ_k.
* Local unknowns of a cell: {k₁, k₂}, kᵢ = 2 (end of the cell in direction i) + (1 for a derivative in direction i),
  as in the reference.
* The polynomials ψ_b, ψ_k. Both implementations bring the same linear system (16 unknowns, the coefficient of
  x₁^m x₂^n in column 4m + n) to reduced row echelon form, take the particular solution with the free unknowns zero
  and one null vector per free column in increasing order, and keep those whose nodal values are independent of each
  other and of those of 1, x₁, x₂ (in that order). The reduced row echelon form is unique, so the parameters P_k
  of a field mean the same in both.
* Voigt order of strains and stresses: (11, 22, 12) with the engineering shear strain ∂₁u₂ + ∂₂u₁.
* Stored claims: both implementations weaken the computed bounds by one part in 10⁹ of the larger modulus.

## 4. Test vectors (`testvectors/`)

`index.json` has one entry per case with the model as a Wolfram Language expression (`wolfram_model`) and the
reference's outputs (`expected`: `lower`, `boundary`, `intQ`, `upper_enclosure`, the number of parameters and of nodes
on the traction-free faces). Each case directory holds `field.json` and the raw arrays `gval.f64`, `U.f64`, `X.f64`,
`par.f64`.

| Case | What it exercises |
|---|---|
| `dirichlet_p1` | all faces Dirichlet, Q₁, 3×2 grid, datum varying on the faces |
| `bodyforce_p2` | all faces Dirichlet, Q₂, body force |
| `block_p1` | two opposite traction-free faces (three parameters) |
| `threefree` | three traction-free faces, anisotropic stiffness, body force, Q₂, box not at the origin |
| `random` | random nodal values and parameters, body force, two opposite traction-free faces |
| `corner` | two adjacent traction-free faces, body force |

The port reproduces `upper_enclosure`, `boundary`, `intQ` and `lower` to better than 10⁻¹¹ relative on all six
(tolerance in the tests: 10⁻⁹).

## 5. Differences between the two implementations

1. **Slope differences on a cell face.** The reference scales the two tangential slopes by the cell length and then
   subtracts; the port subtracts and then scales. Both enclose the same number.
2. **Sums.** The reference sums the cell enclosures with `math.fsum` and two outward steps; the port sums them
   exactly in rational arithmetic and combines the boundary term and the cell term exactly.
3. **Nodal values on the traction-free faces.** The reference encloses ψ_b + Σ P_k ψ_k by interval arithmetic; the
   port evaluates it exactly in rational arithmetic and encloses the result.
4. **Heuristic solves.** The reference assembles the two sparse matrices cell by cell; the port builds them from
   one-dimensional matrices with Kronecker products. Both solve the same linear systems.
5. **Rounding.** The reference widens each operation by one unit in the last place, the port by 2⁻⁴⁹ relative plus
   10⁻³⁰⁰ (`RigorousArithmetic.wl`); every result of an interval operation is tested for machine overflow
   (`ndChk` of `ScalarFieldND.wl`).

## 6. Log

* 23 November 2026. First execution: the rigorous evaluators agreed with the reference on all six test vectors
  at once (differences below 5×10⁻¹²), and the public interface reproduced the block and gravity examples of the
  reference.
* 23 November 2026. **Cancellation found and removed in both implementations.** Re-verification of the stored
  16×16 and 32×32 certificates in the Wolfram Language failed by 10⁻¹⁰. Cause: both evaluators formed the quadratic
  forms X·A·X and u·K·u from the nodal values directly. The products are then larger than the result by the factor
  (|φ|/(h²|σ|))² (stress) and (|u|/(h|∇u|))² (displacement), and outward rounding of the cancelling terms costs
  accuracy that grows like h⁻⁴ under refinement. The port widens every operation by 2⁻⁴⁹, sixteen times more than
  the reference, so it showed first (10⁻¹⁰ in the lower bound of the 16×16 block, which is 0.0150172); the
  reference lost 2×10⁻¹⁰ on the 32×32 block and 10⁻⁹ on the 16×32 column (lower bound −0.4514566). Never a wrong
  bound, only a weaker one. Fix, in both: work with
  differences. The Airy nodal values of a cell minus the affine function through the three values at its lower
  corner, the displacement nodal values minus the value at the lower corner, and on a cell face the corresponding
  one-dimensional differences. This is exact, because the element matrices annihilate affine Airy functions and
  constant displacements; both facts are checked in rational arithmetic when the local matrices are built (an
  assertion in the reference, a `Failure` in the port). `check_exact.py`, `check_theorem.py` and `verify_all.py`
  pass; the stored certificates and test vectors were regenerated (the arrays are bit-identical, the claims and
  expected values changed in the last digits).
* 23 November 2026. All stored certificates of the reference are recomputed in the Wolfram Language with bounds at
  least as strong as their claims. Test suite written: 90 tests pass.

## 7. Observations

* The first row of the results table in `docs/ELASTICITY.md` had its lower bound rounded to nearest
  (0.01500031) instead of down (0.01500030). The table is now produced with outward rounding.
* A reference in three dimensions that forms X·A·X from the nodal values directly will lose accuracy in the same
  way, and sooner (more terms, the same factor h⁻⁴).
