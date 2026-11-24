# Linear elasticity on boxes in three dimensions: the certified bounds

This note extends `docs/ELASTICITY.md` (two dimensions) to boxes in ℝ³. It states and proves what the reference
implementation `Reference/Elasticity/elasticity3d.py` certifies and describes how the two bounds are computed.

**Status.** Python reference only; there is no Wolfram Language model for three dimensions yet. The inequality is
the complementary-energy principle of Prager and Synge (1947); stress functions for the equilibrium equations in
three dimensions go back to Maxwell, Morera and Beltrami, and bounding the energy from both sides with a
displacement model and an equilibrium model is standard in structural mechanics ("dual analysis"). No novelty is
claimed for any of this. The short proofs below have been checked by their author only. What the reference adds is
exact static admissibility of the stress field in floating-point data, outward-rounded evaluation, and stored
certificates that can be re-checked. Numbers quoted as *certified* come from outward-rounded interval arithmetic on
exact inputs; numbers quoted as *exact* are rational numbers. Statements marked [SKETCH] concern how good the bound
can be, never whether it is valid.

## 1. The problem

* Ω = (a₁,b₁) × (a₂,b₂) × (a₃,b₃), divided into N₁ × N₂ × N₃ congruent boxes (*cells*) with edge lengths hᵢ.
* Each of the six faces of Ω is a *Dirichlet face* or a *traction-free face*. Γ_D is the union of the Dirichlet
  faces and is not empty; Γ_N is the union of the traction-free faces.
* Displacements u : Ω → ℝ³. Strains and stresses are written as vectors of six components in the order
  (11, 22, 33, 23, 13, 12):

      e(u) = (∂₁u₁, ∂₂u₂, ∂₃u₃, ∂₂u₃ + ∂₃u₂, ∂₁u₃ + ∂₃u₁, ∂₁u₂ + ∂₂u₁),      s = (σ₁₁, σ₂₂, σ₃₃, σ₂₃, σ₁₃, σ₁₂).

* The stiffness is a constant symmetric positive definite 6×6 matrix D with rational entries, s = D e, and the
  energy density is ½ eᵀDe = ½ σ : ε. For an isotropic material the upper left 3×3 block of D has λ + 2μ on the
  diagonal and λ elsewhere, the lower right block is μ times the identity, and the other entries vanish (μ > 0,
  3λ + 2μ > 0). S = D⁻¹ is the compliance.
* b ∈ ℚ³ is a constant body force.
* g is the Dirichlet datum on Γ_D: **the trace of the continuous function that is a polynomial of degree p (p = 1 or
  2) in each variable on every cell and takes the given nodal values at the nodes of Γ_D.** It coincides with the
  user's function exactly when that function is such a polynomial on each boundary cell face.

Admissible class and energy:

    𝒜 = { u ∈ H¹(Ω; ℝ³) : u = g on Γ_D (trace) },        E[u] = ∫_Ω ½ e(u)ᵀ D e(u) − b·u dx.

## 2. The lower bound

A field σ of symmetric 3×3 matrices is *statically admissible* if

1. σ is a polynomial on every cell,
2. the traction σν is single-valued on every interior cell face,
3. div σ + b = 0 in every cell,
4. σν = 0 on Γ_N.

**Theorem 1 (Prager–Synge).** For every statically admissible σ and every u ∈ 𝒜,

    E[u]  ≥  LB(σ) := ∫_{Γ_D} (σν)·g dℋ²  −  ½ ∫_Ω sᵀ S s dx.

*Proof.* The proof of Theorem 1 of `docs/ELASTICITY.md` applies word for word: ½eᵀDe − s·e + ½sᵀSs =
½(e − Ss)ᵀD(e − Ss) ≥ 0 pointwise, s·e(u) = σ : ∇u because σ is symmetric, and cell by cell
∫_c σ : ∇u = ∫_{∂c} (σν_c)·u + ∫_c b·u by 3. In the sum over the cells the interior faces cancel by 2 and because
the trace of u is single-valued, Γ_N drops out by 4, and on Γ_D the trace of u is g. ∎

As in two dimensions, for u_h ∈ 𝒜 the difference is E[u_h] − LB(σ) = ½∫(e(u_h) − Ss)ᵀ D (e(u_h) − Ss), and it
equals ½‖e(u_h) − e(u*)‖²_D + ½‖s − s*‖²_S with the minimiser u* and its stress s*.

## 3. Statically admissible stresses from stress functions

### 3.1 The space

Let 𝒮ᵢ be the space of C¹ functions on [aᵢ,bᵢ] that are cubic polynomials on each cell interval of direction i, and
V_h = 𝒮₁ ⊗ 𝒮₂ ⊗ 𝒮₃ the span of the products p(x₁)q(x₂)r(x₃) (tricubic Hermite elements). A *type* is a
d ∈ {0,1}³, and ∂^d = ∂₁^{d₁}∂₂^{d₂}∂₃^{d₃}. A function φ ∈ V_h is determined by its eight *nodal values* ∂^dφ(n)
at every node n: φ = Σ_n Σ_d ∂^dφ(n) β_{n₁}^{d₁}(x₁) β_{n₂}^{d₂}(x₂) β_{n₃}^{d₃}(x₃), where β_j⁰ and β_j¹ are the
one-dimensional functions in 𝒮ᵢ with (β_j^d)^{(m)}(x_{j′}) = δ_{jj′} δ_{dm} for d, m ∈ {0,1}.

* **Fact A.** For φ ∈ V_h and α ∈ {0,1}³ the derivative ∂^αφ is continuous on Ω̄ (each product
  p^{(α₁)}q^{(α₂)}r^{(α₃)} is).
* **Fact B.** Let x_i = c be a grid plane and m ∈ {0,1}. The restriction of ∂ᵢ^mφ to the plane is the element of
  the tensor product of the two other spaces 𝒮 with the nodal values { ∂^dφ(n) : d_i = m } at the nodes n of the
  plane. (Insert x_i = c in the sum above and use the property of the β.)
* **Fact C.** A polynomial of degree at most three in each variable belongs to V_h, and the corresponding statement
  holds on a grid plane. Hence two such polynomials with the same nodal values at the nodes of a grid plane (the
  four values with d_i = 0) agree on that plane.

### 3.2 The stress field

Unknowns are six functions Φ = (Φ₁₁, Φ₂₂, Φ₃₃, Φ₂₃, Φ₁₃, Φ₁₂) in V_h, the entries of a symmetric matrix. Let σ^Φ be
the field with components σ_ij = Σ ε_ikl ε_jmn ∂_k∂_m Φ_ln (the representation of Beltrami), that is

    σ₁₁ = ∂₂₂Φ₃₃ + ∂₃₃Φ₂₂ − 2∂₂₃Φ₂₃            σ₂₃ = ∂₁₃Φ₁₂ + ∂₁₂Φ₁₃ − ∂₂₃Φ₁₁ − ∂₁₁Φ₂₃
    σ₂₂ = ∂₃₃Φ₁₁ + ∂₁₁Φ₃₃ − 2∂₁₃Φ₁₃            σ₁₃ = ∂₁₂Φ₂₃ + ∂₂₃Φ₁₂ − ∂₁₃Φ₂₂ − ∂₂₂Φ₁₃
    σ₃₃ = ∂₁₁Φ₂₂ + ∂₂₂Φ₁₁ − 2∂₁₂Φ₁₂            σ₁₂ = ∂₂₃Φ₁₃ + ∂₁₃Φ₂₃ − ∂₁₂Φ₃₃ − ∂₃₃Φ₁₂

(each line is obtained from the one above by the cyclic permutation 1 → 2 → 3 → 1). With Φ₂₃ = Φ₁₃ = Φ₁₂ = 0 these
are the stress functions of Maxwell, with Φ₁₁ = Φ₂₂ = Φ₃₃ = 0 those of Morera. Let σ^b be the field
−diag(b₁x₁, b₂x₂, b₃x₃).

**Lemma 2.** σ = σ^Φ + σ^b satisfies conditions 1, 2 and 3.

*Proof.* 1 is clear. For 3: in a cell the Φ are polynomials and

    ∂₁σ₁₁ + ∂₂σ₁₂ + ∂₃σ₁₃ = (∂₁₂₂Φ₃₃ + ∂₁₃₃Φ₂₂ − 2∂₁₂₃Φ₂₃) + (∂₂₂₃Φ₁₃ + ∂₁₂₃Φ₂₃ − ∂₁₂₂Φ₃₃ − ∂₂₃₃Φ₁₂)
                              + (∂₁₂₃Φ₂₃ + ∂₂₃₃Φ₁₂ − ∂₁₃₃Φ₂₂ − ∂₂₂₃Φ₁₃) = 0;

the other two rows follow by the cyclic permutation, and div σ^b = −b. For 2: σ^b is continuous. On a cell face in
a plane {x₁ = c} the traction of σ^Φ is ±(σ₁₁, σ₁₂, σ₁₃). Reading the formulas, every term of these three components
is a derivative in the directions x₂, x₃ (of total order two or one) of one of the functions

    Φ₂₂, Φ₃₃, Φ₂₃, Φ₁₃, Φ₁₂,   ∂₁Φ₂₂, ∂₁Φ₃₃, ∂₁Φ₂₃                                                       (R₁)

restricted to the plane: for example ∂₁₃Φ₂₃ = ∂₃(∂₁Φ₂₃) and ∂₃₃Φ₁₂. By Fact A these eight functions are
continuous across the face, so the polynomials on the two sides have the same restriction to the face and
therefore the same derivatives along it. The planes {x₂ = c} and {x₃ = c} follow by the cyclic permutation, with
the lists

    Φ₃₃, Φ₁₁, Φ₁₃, Φ₁₂, Φ₂₃,   ∂₂Φ₃₃, ∂₂Φ₁₁, ∂₂Φ₁₃                                                       (R₂)
    Φ₁₁, Φ₂₂, Φ₁₂, Φ₂₃, Φ₁₃,   ∂₃Φ₁₁, ∂₃Φ₂₂, ∂₃Φ₁₂.                                                      (R₃)  ∎

(In words: on a plane with normal eᵢ the traction involves the components of Φ without the index i together with
their normal derivatives, the components with exactly one index i without normal derivative, and Φ_ii not at all.)

### 3.3 Traction-free faces

**Lemma 3.** Let F be a face of Ω in the plane {x_i = c} and let Ψ = (Ψ₁₁, …, Ψ₁₂) be six polynomials of degree at
most three in each variable. Suppose that at every node n of F

    ∂^dΦ_f(n) = ∂^dΨ_f(n)      for all components f and types d such that ∂ᵢ^{d_i}Φ_f occurs in the list (Rᵢ).

Then σ^Φν = σ^Ψν on F.

*Proof.* By Fact C, Ψ_f ∈ V_h. By Fact B the restrictions to the plane of the functions in the list (Rᵢ), formed
with Φ and with Ψ, are elements of the tensor product space of the plane, and they have the same nodal values by
assumption; hence they coincide on F. By the proof of Lemma 2 the traction on F is a combination of derivatives
along F of these restrictions. ∎

**Corollary 4 (how condition 4 is met exactly).** Let Ψ_b, Ψ₁, …, Ψ_r be sextuples of polynomials of degree at most
three in each variable with rational coefficients such that, as identities between polynomials on every
traction-free face,

    σ^{Ψ_b}ν = −σ^bν        and        σ^{Ψ_k}ν = 0   (k = 1, …, r),

and let P₁, …, P_r be real numbers. If the nodal values of Φ named in Lemma 3 are, at every node of every
traction-free face, those of Ψ_b + P₁Ψ₁ + … + P_rΨ_r, then σ^Φ + σ^b satisfies condition 4. All other nodal values
of Φ are arbitrary.

*Proof.* Lemma 3 with Ψ = Ψ_b + ΣP_kΨ_k, and linearity. ∎

In the implementation the P_k and the free nodal values are machine numbers, and the prescribed nodal values are
the exact rational expressions Ψ_b + ΣP_kΨ_k in them; they are never rounded. The lower bound of Theorem 1 is
then valid for every choice of the machine numbers.

### 3.4 The polynomials Ψ (*states*)

They are found in exact rational arithmetic and depend on the box, on which faces are traction-free and on b, not
on the grid. A sextuple of polynomials of degree at most three in each variable has 6·64 coefficients.

* The conditions of Corollary 4 are linear equations for the coefficients (equate the coefficients of the
  monomials in the two variables of each traction-free face). Ψ_b is the solution with all free unknowns equal to
  zero in the reduced row echelon form of this system. If the system has no solution the problem is refused.
* 𝒫 is the space of solutions of the homogeneous system, 𝒩 ⊂ 𝒫 the space of sextuples with σ^Ψ = 0 identically
  (*null sextuples*). Adding an element of 𝒩 to Φ does not change σ^Φ. Ψ₁, …, Ψ_r are elements of the echelon basis
  of 𝒫 whose face data (the polynomials in the list (Rᵢ) on each traction-free face) are linearly independent of
  each other and of the face data of 𝒩, as many as there are.

Any choice of the Ψ_k satisfying the identities of Corollary 4 gives valid bounds; the choice only affects how
good the best bound in the family can be. Section 3.5 says what is expected.

### 3.5 What the states are for (tightness only) [SKETCH]

Nothing in this subsection is used for the validity of a bound.

* *Gauge.* σ^Φ = 0 for Φ = ½(∇w + ∇wᵀ) with any vector field w; these fields change Φ and not the stress. This is
  why the Hessian of the dual problem is singular, and why the nodal values at a traction-free face can often be
  taken to be zero: if the stress vanishes on the far side of a face (extend it by zero across the traction-free
  face), Φ is of this form there, and subtracting a field of this form defined everywhere leaves a Φ that vanishes
  with its first derivatives on the face.
* *Rings.* That argument needs w to be single-valued in the region outside Γ_N. If Γ_N consists of the four faces
  around one axis and the two end faces are Dirichlet faces, that region is a ring, and w can have a jump which is
  a rigid motion (six parameters). Correspondingly r = 6 is expected: six states carrying the resultant force and
  moment transmitted through the body. Examples for a bar along x₁ with the cross-section (−w,w) × (−h,h): tension,
  Φ₃₃ = x₂²/2; bending, Φ₃₃ = x₂³/6; torsion, Φ₃₃ = −2x₁x₃(x₂³/3 − w²x₂), Φ₂₂ = 2x₁x₂(x₃³/3 − h²x₃).
  For every other arrangement of traction-free faces r = 0 is expected, and then all prescribed nodal values are
  those of Ψ_b alone.
* *Why six functions.* With the three functions of Maxwell alone the fields that change Φ and not the stress are
  far fewer (w₁ may depend on x₂ and x₃ only in a restricted way), and between two opposite traction-free faces
  the prescribed nodal values then allow shear resultants ∫σ₁₂dx₂, ∫σ₂₃dx₂ that are functions of one variable only.
  A bound computed with Maxwell's functions alone is valid, but for bending-type loads it does not converge to
  the infimum under refinement. Section 8 shows this numerically.

### 3.6 The certified lower bound

For a field as in Corollary 4, LB(σ) is a polynomial of degree two in the nodal values with rational coefficients:
on each cell sᵀSs is a product of polynomials, and on each Dirichlet cell face (σν)·g is a product of polynomials.
The coefficients are products of one-dimensional integrals of products of the cubic Hermite and Lagrange basis
functions and their derivatives, computed exactly. `lower_bound` evaluates the polynomial with outward-rounded
interval arithmetic; the prescribed nodal values on Γ_N enter as enclosures of their exact values, all others as
the exact binary numbers they are.

## 4. The upper bound

Let u_h be the continuous function that is a polynomial of degree p in each variable on every cell (p = 1 or 2),
with the given nodal values of g at the nodes of Γ_D. Then u_h ∈ 𝒜, and E[u_h] is a polynomial of degree two in the
nodal values with rational coefficients. `upper_bound` returns an enclosure of it; its upper end is the certified
upper bound for inf E. Any nodal values give a valid bound.

## 5. Choosing the fields (heuristic)

Both bounds are quadratic in their unknowns. The candidate is the finite-element solution in Q_p, the stress field
the maximiser of LB over the free nodal values and the parameters P_k. Both are computed by preconditioned
conjugate gradients in ordinary floating point, the second on a singular system (Section 3.5) with a small
regularisation. Neither is trusted; the nodal values are then taken as the exact binary numbers they are.

## 6. Arithmetic and trusted base

As in two dimensions: rigorous code uses +, −, × of IEEE-754 doubles with outward rounding, comparisons, and
exact rational arithmetic (the coefficients of the two quadratic polynomials, the grid, D, S, b, the polynomials
Ψ). No linear solver and no library function enters a certified number. Trusted: IEEE-754 arithmetic; exact
rational arithmetic; this note.
