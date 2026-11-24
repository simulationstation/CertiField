# Linear elasticity on boxes: the certified bounds

This note states and proves what the model `LinearElasticity` will certify, and describes how the two bounds are
computed. It is the analytic reduction behind `Reference/Elasticity/elasticity.py` (the Python reference).

**State of the work (23 November 2026).** The Python reference for two dimensions is complete and tested
(Sections 7 and 8). The Wolfram Language model (`Kernel/Elasticity.wl`, registration, `VerifyCertificate`, tests
against `Reference/Elasticity/testvectors`) is not written yet, and three dimensions are not started.

**Status.** Two dimensions (plane strain). The inequality is the complementary-energy principle of Prager and Synge
(1947) and the construction of the stress field by an Airy function in a C¹ finite-element space is classical; no
novelty is claimed for either. The short proofs below have been checked by their author only. What the package adds
is exact static admissibility of the stress field in floating-point data, outward-rounded evaluation, and a
certificate that can be re-checked. Numbers quoted as *certified* come from outward-rounded interval arithmetic on
exact inputs; numbers quoted as *exact* are rational numbers. Three dimensions are not implemented (Section 9).

## 1. The problem

* Ω = (a₁,b₁) × (a₂,b₂), divided into N₁ × N₂ congruent rectangles (*cells*) with edge lengths hᵢ.
* Each of the four faces of Ω is a *Dirichlet face* or a *traction-free face*. Γ_D is the union of the Dirichlet
  faces and is not empty; Γ_N is the union of the traction-free faces.
* Displacements u : Ω → ℝ², strains ε(u) = ½(∇u + ∇uᵀ), written as e(u) = (∂₁u₁, ∂₂u₂, ∂₁u₂ + ∂₂u₁).
* The stiffness is a constant symmetric positive definite 3×3 matrix D with rational entries: the stress in the
  same notation is s = (σ₁₁, σ₂₂, σ₁₂) = D e, and the energy density is ½ eᵀDe. For an isotropic material in plane
  strain, D = [[λ+2μ, λ, 0], [λ, λ+2μ, 0], [0, 0, μ]] with μ > 0 and λ + μ > 0. S = D⁻¹ is the compliance.
* b ∈ ℚ² is a constant body force.
* g is the Dirichlet datum on Γ_D. In the implementation it is **the trace of the piecewise polynomial function of
  degree p (p = 1 or 2) in each variable that interpolates the given nodal values on the boundary grid**; it
  coincides with the user's function exactly when that function is such a polynomial on each boundary cell face.

Admissible class and energy:

    𝒜 = { u ∈ H¹(Ω; ℝ²) : u = g on Γ_D (trace) },        E[u] = ∫_Ω ½ e(u)ᵀ D e(u) − b·u dx.

Nothing is prescribed on Γ_N (no traction acts there).

## 2. The lower bound

Call a field of symmetric 2×2 matrices σ, with s = (σ₁₁, σ₂₂, σ₁₂), *statically admissible* if

1. σ is a polynomial on every cell,
2. the normal traction σν is single-valued on every interior cell face,
3. div σ + b = 0 in every cell,
4. σν = 0 on Γ_N.

**Theorem 1 (Prager–Synge).** For every statically admissible σ and every u ∈ 𝒜,

    E[u]  ≥  ∫_{Γ_D} (σν)·g dℋ¹  −  ½ ∫_Ω sᵀ S s dx.

The right-hand side does not depend on u. Any statically admissible σ gives a valid bound; a poor one gives a poor
bound.

*Proof.* Pointwise, ½ eᵀDe − s·e + ½ sᵀSs = ½ (e − Ss)ᵀ D (e − Ss) ≥ 0, and s·e(u) = σ : ∇u because σ is symmetric.
Hence ½ e(u)ᵀDe(u) ≥ σ : ∇u − ½ sᵀSs almost everywhere. On a cell c, σ is smooth and u ∈ H¹(c), so
∫_c σ : ∇u = ∫_{∂c} (σν_c)·u − ∫_c (div σ)·u = ∫_{∂c} (σν_c)·u + ∫_c b·u by 3. Sum over the cells. Every interior
face occurs twice with opposite normals; by 2 and because the trace of u is single-valued these contributions
cancel. On Γ_N the traction vanishes by 4, and on Γ_D the trace of u is g. Therefore
∫_Ω σ : ∇u = ∫_{Γ_D} (σν)·g + ∫_Ω b·u, and

    E[u] ≥ ∫_{Γ_D} (σν)·g + ∫_Ω b·u − ½ ∫_Ω sᵀSs − ∫_Ω b·u.  ∎

The two terms with b·u cancel only because equilibrium holds exactly. A residual div σ + b ≠ 0, however small, would
leave a term that is unbounded on 𝒜. This is why conditions 3 and 4 are enforced exactly and not to rounding.

**Remark (what the gap measures).** If u_h ∈ 𝒜 and σ is statically admissible, the same computation gives
E[u_h] − (right-hand side) = ½ ∫ (e(u_h) − Ss)ᵀ D (e(u_h) − Ss), and with the minimiser u* and σ* = D e(u*),

    E[u_h] − LB(σ) = ½ ‖e(u_h) − e(u*)‖²_D + ½ ‖s − s*‖²_S

(the hypercircle identity). The certified gap is therefore also a certified bound for the energy-norm errors of the
candidate and of the stress field.

## 3. Statically admissible stresses from an Airy function

**The space.** V_h is the space of tensor products of C¹ piecewise cubic functions on the two grids (the
Bogner–Fox–Schmit element): on each cell a polynomial of degree at most three in each variable, determined by the
values of φ, ∂₁φ, ∂₂φ and ∂₁₂φ at the nodes. Every φ ∈ V_h is a finite sum Σ p(x₁) q(x₂) with p, q ∈ C¹ piecewise
cubic, so φ, ∂₁φ, ∂₂φ and ∂₁₂φ are continuous on Ω̄.

For φ ∈ V_h let σ^φ be the field with σ₁₁ = ∂₂₂φ, σ₂₂ = ∂₁₁φ, σ₁₂ = −∂₁₂φ, and let σ^b be the field with
σ₁₁ = −b₁x₁, σ₂₂ = −b₂x₂, σ₁₂ = 0.

**Lemma 2.** σ = σ^φ + σ^b satisfies conditions 1, 2 and 3.

*Proof.* 1 is clear. In a cell φ is a polynomial, so div σ^φ = (∂₁∂₂₂φ − ∂₂∂₁₂φ, −∂₁∂₁₂φ + ∂₂∂₁₁φ) = 0, and
div σ^b = −b: this is 3. For 2, σ^b is continuous. On a cell face in {x₁ = const} the traction of σ^φ is
±(∂₂₂φ, −∂₁₂φ) = ±(∂₂(∂₂φ), −∂₂(∂₁φ)): derivatives along the face of the functions ∂₂φ and ∂₁φ, which are
continuous across the face and whose restrictions to it are therefore the same from both sides. The faces in
{x₂ = const} are treated in the same way. ∎

**Lemma 3 (traction-free faces).** Let F ⊂ {x₁ = a} be a face of Ω. The traction of σ^φ + σ^b vanishes on F if and
only if there are constants α, β, γ with

    φ(a, x₂) = b₁a x₂²/2 + α + βx₂      and      ∂₁φ(a, x₂) = γ      on F,

and this holds if and only if the nodal values on F are those of these two polynomials: φ = b₁a x₂²/2 + α + βx₂,
∂₂φ = b₁a x₂ + β, ∂₁φ = γ, ∂₁₂φ = 0 at every node of F. The same holds on a face in {x₂ = c} with the indices
exchanged.

*Proof.* The traction is ±(∂₂₂φ − b₁a, −∂₂(∂₁φ)). It vanishes on F exactly when ∂₂₂φ(a,·) = b₁a and ∂₁φ(a,·) is
constant, which is the first statement. On F, φ(a,·) is the C¹ piecewise cubic function with the nodal values and
slopes (φ, ∂₂φ), and ∂₁φ(a,·) is the one with (∂₁φ, ∂₁₂φ). A C¹ piecewise cubic function equals a given polynomial
of degree at most three on all of F if and only if its nodal values and slopes are those of the polynomial
(uniqueness of cubic Hermite interpolation on each cell). ∎

**How condition 4 is met exactly.** The nodal values of Lemma 3 are in general not machine numbers (they contain
products such as βx₂), so they cannot simply be stored. The implementation computes, in rational arithmetic,

* one polynomial ψ_b of degree at most three in each variable with ∂₂₂ψ_b = b₁a, ∂₁₂ψ_b = 0 on every
  traction-free face in {x₁ = a}, and correspondingly on the others (ψ_b = 0 if b = 0 or Γ_N = ∅), and
* polynomials ψ₁, …, ψ_r of the same kind with zero right-hand sides, chosen so that their nodal values on Γ_N,
  together with those of the affine functions, are linearly independent,

and defines the field by: the four nodal values at a node on Γ_N are those of ψ_b + P₁ψ₁ + … + P_rψ_r, with
parameters P_k that are machine numbers; the nodal values at the other nodes are machine numbers. The nodal values
on Γ_N are thus exact rational expressions in machine numbers, and by Lemma 3 the traction vanishes on Γ_N
exactly. (In two dimensions nothing is lost: by Lemma 3 the admissible face data are polynomials of degree at most
two with compatible values at shared corners, and these are the traces of polynomials of degree at most three in
each variable. With one connected traction-free part of the boundary r = 0; with two, r = 3: the resultant force
and moment carried through the body.) Adding an affine function to φ does not change σ^φ; if Γ_N = ∅ the nodal
values φ, ∂₁φ, ∂₂φ at one node are set to zero.

**The certified lower bound.** For such a field,

    LB(σ) = ∫_{Γ_D} (σν)·g − ½ ∫_Ω sᵀSs

is a polynomial of degree two in the nodal values and the parameters, with rational coefficients: on each cell
sᵀSs is a product of polynomials, and on each Dirichlet cell face (σν)·g is a product of polynomials, so the
integrals are exact combinations of the integrals of products of the one-dimensional basis functions and their
derivatives. `lower_bound` evaluates this polynomial with outward-rounded interval arithmetic; the nodal values on
Γ_N enter as enclosures of their exact values, all others as the exact binary numbers they are.

## 4. The upper bound

Let u_h be the continuous function that is a polynomial of degree p in each variable on every cell (Q_p, p = 1 or
2), with the given nodal values of g at the nodes of Γ_D. Then u_h ∈ 𝒜, and E[u_h] is a polynomial of degree two
in the nodal values with rational coefficients (products of polynomials integrated over cells). `upper_bound`
returns an enclosure of it; its upper end is the certified upper bound for inf E. Any nodal values give a valid
bound.

## 5. Choosing the fields (heuristic)

Both bounds are quadratic in their unknowns, so each field is the solution of one sparse linear system in ordinary
floating point: the finite-element solution in Q_p for the candidate, and the maximiser of LB over the free nodal
values and parameters for the stress field. Neither is trusted. The nodal values are then taken as the exact
binary numbers they are.

## 6. Arithmetic and trusted base

Rigorous code uses +, −, × of IEEE-754 doubles with outward rounding, comparisons, and exact rational arithmetic
(the coefficients of the two quadratic polynomials, the grid, D, S, b and the polynomials ψ are rational). No
linear solver and no library function enters a certified number. Trusted: IEEE-754 arithmetic; exact rational
arithmetic; this note.

## 7. Checks

All with the Python reference (`Reference/Elasticity`).

* **Exact cross-check** (`check_exact.py`). For seven problems (all faces Dirichlet; one, two adjacent, two
  opposite and three traction-free faces; with and without body force; isotropic and anisotropic stiffness; p = 1
  and 2), with the computed fields and with random ones, both fields are rebuilt cell by cell as polynomials in
  rational arithmetic, with no code shared with the evaluators. Verified exactly: div σ + b = 0 in every cell; the
  traction is the same from both sides of every interior cell face; it vanishes identically on the traction-free
  faces; the candidate is continuous. The exact values of both bounds lie inside the computed enclosures, which
  are 10⁻¹² wide. The identity E[u_h] − LB(σ) = ½∫(e − Ss)ᵀD(e − Ss) of Section 2 holds as an equality of rational
  numbers, which it would not if the integration by parts of the proof failed for these fields.
* **Known minima** (`check_theorem.py`, part A). Homogeneous strain with all faces Dirichlet; uniaxial tension with
  free sides; a quadratic displacement with constant body force; pure bending with free sides. Data, grids and
  compliance are dyadic, so the stored nodal values are exactly the data. The exact rational minimum lies inside
  the certified interval in every case, and both bounds equal it to 10⁻¹² (the exact fields are in the spaces).
* **Falsification** (`check_theorem.py`, part B). Random and perturbed stress fields and displacements for four
  problems: no pair has LB(σ) > E[u]; the smallest difference found is 1.7×10⁻⁵.
* **Stored certificates** (`verify_all.py`) recompute both bounds from the stored fields.

## 8. Results (Python reference)

Unit square, λ = 3/2, μ = 1 (Poisson ratio 0.3), clamped on the left face, displaced by (0.1, 0) on the right face,
no traction on the other two (the corners where the two kinds of face meet carry stress singularities):

| Grid | p | Certified: lower ≤ inf E ≤ upper | Gap | Time |
|---|---|---|---|---|
| 8×8 | 1 | 0.01500031 ≤ inf E ≤ 0.01512044 | 0.794 % | 0.1 s |
| 16×16 | 1 | 0.01501720 ≤ inf E ≤ 0.01506422 | 0.313 % | 0.1 s |
| 32×32 | 1 | 0.01502346 ≤ inf E ≤ 0.01504168 | 0.122 % | 0.2 s |
| 64×64 | 1 | 0.01502578 ≤ inf E ≤ 0.01503277 | 0.047 % | 1.2 s |
| 16×16 | 2 | 0.01501720 ≤ inf E ≤ 0.01503404 | 0.112 % | 0.1 s |
| 32×32 | 2 | 0.01502346 ≤ inf E ≤ 0.01502974 | 0.042 % | 0.4 s |
| 64×64 | 2 | 0.01502578 ≤ inf E ≤ 0.01502813 | 0.016 % | 1.9 s |

The gap falls by a factor of about 2.6 per refinement, the rate h^1.4 that the corner singularities allow. A
column under its own weight, clamped at the bottom with three free faces (16×32, p = 2):
−0.45145662 ≤ inf E ≤ −0.45112300, gap 0.074 %.

## 9. What is not covered

* **Not claimed:** that the candidate is the minimiser. The gap bounds E[u_h] − inf E, and by the remark in
  Section 2 also the energy-norm error of the candidate.
* Three dimensions. The same theorem holds with Beltrami stress functions (σ = curl (curl Φ)ᵀ with Φ a symmetric
  matrix of C¹ piecewise polynomials: the traction on a cell face then contains only tangential derivatives of Φ
  and of its first derivatives, so it is continuous). What is missing is the exact treatment of traction-free
  faces, where the admissible traces are no longer a finite-dimensional family as in Lemma 3.
* Conditions that prescribe only some components of u on a face (rollers), prescribed non-zero tractions, body
  forces that are not constant, and stiffness that varies in space.
* Domains other than rectangles.
* Corners where a Dirichlet face meets a traction-free face carry stress singularities. Both bounds remain valid,
  but the gap closes slowly under refinement there.
* Nearly incompressible materials (λ/μ large): the bounds remain valid; the candidate in Q_p locks and the gap
  grows.

## References

1. W. Prager, J. L. Synge, *Approximations in elasticity based on the concept of function space*, Quart. Appl. Math.
   5 (1947), 241–269. (The inequality of Theorem 1 and the hypercircle identity.)
2. B. Fraeijs de Veubeke, *Displacement and equilibrium models in the finite element method*, in Stress Analysis
   (O. C. Zienkiewicz, G. S. Holister, eds.), Wiley, 1965, 145–197; *Upper and lower bounds in matrix structural
   analysis*, AGARDograph 72, Pergamon, 1964, 165–201. (Dual analysis: bounding the energy from both sides with
   a displacement model and an equilibrium model.)
3. J. P. Moitinho de Almeida, E. A. W. Maunder, *Equilibrium Finite Element Formulations*, Wiley, 2017.
4. F. K. Bogner, R. L. Fox, L. A. Schmit, *The generation of interelement compatible stiffness and mass matrices by
   the use of interpolation formulas*, Proc. Conf. Matrix Methods in Structural Mechanics, Wright-Patterson AFB,
   1965. (The C¹ rectangle used for the Airy function; cited from memory.)
5. F. Key, L. Freinberger, *A spline-based stress function approach for the principle of minimum complementary
   energy*, Comput. Methods Appl. Mech. Engrg. 448 (2025) 118492, arXiv:2506.19534. (Stress functions in smooth
   spline spaces for the same principle, without certification.)

A brief web search (November 2026) confirmed the bibliographic data of items 1, 2, 3 and 5 and found the method
described as standard under the names "dual analysis" and "equilibrium finite elements". It did not find an
implementation that enforces static admissibility exactly in floating-point data and evaluates the bound with
rigorous arithmetic. That search is not a literature review.
