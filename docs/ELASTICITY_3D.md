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
on the grid. A sextuple of polynomials of degree at most three in each variable has 6·64 coefficients; the
coefficient of x₁^γ₁ x₂^γ₂ x₃^γ₃ in the f-th function has the number 64f + 16γ₁ + 4γ₂ + γ₃.

* The conditions of Corollary 4 are linear equations for the coefficients (equate the coefficients of the
  monomials in the two variables of each traction-free face). The reduced row echelon form of this system is
  computed; it is unique. Ψ_b is the solution with all free unknowns equal to zero. If the system has no solution
  the problem is refused; this does not happen for any of the 63 arrangements of the faces (Section 7). The space
  𝒫 of solutions of the homogeneous system gets the echelon basis, one vector per free unknown in increasing
  order.
* **Default rule.** Call Γ_N a *ring* if it consists of the four faces around one axis x_a and the two faces normal
  to e_a are Dirichlet faces. For a ring, Ψ₁, …, Ψ₆ are the first six elements of the echelon basis whose
  *resultants* are linearly independent: the force ∫σe_a and the moment ∫x × σe_a over a cross-section, six
  rational numbers for each element, computed exactly. In every other case r = 0, and the prescribed nodal values
  are those of Ψ_b alone (zero if b = 0).
* **Alternative rule** (`states='all'`, used for comparison only): all elements of the echelon basis whose face
  data (the polynomials in the list (Rᵢ) on each traction-free face) are linearly independent of each other and of
  the face data of the sextuples with σ^Ψ = 0 identically. For this rule the ranks are computed modulo two primes.

Any choice of the Ψ_k satisfying the identities of Corollary 4 gives valid bounds, whatever the rule; a wrong
rank could only change which valid states are used. The choice affects how good the best bound in the family can
be, and Section 3.5 says what is expected.

### 3.5 What the states are for (tightness only) [SKETCH]

Nothing in this subsection is used for the validity of a bound.

* *Gauge.* σ^Φ = 0 for Φ = ½(∇w + ∇wᵀ) with any vector field w; these fields change Φ and not the stress. This is
  why the system for the stress field is singular, and why the nodal values at a traction-free face can be taken
  to be zero: if the stress is extended by zero across the face, Φ is of this form on the far side, and subtracting
  a field of this form defined everywhere leaves a Φ that vanishes with its first derivatives on the face.
* *Rings.* That argument needs w to be single-valued in the region outside Γ_N. Outside a ring it need not be:
  going once around, w can change by a rigid motion, which has six parameters (the classical statement for
  multiply connected bodies is due to Volterra). Correspondingly six states are needed, carrying the resultant
  force and moment through the cross-section. Examples for a bar along x₁ with the cross-section (−w,w) × (−h,h):
  tension, Φ₃₃ = x₂²/2; bending, Φ₃₃ = x₂³/6; torsion, Φ₃₃ = −2x₁x₃(x₂³/3 − w²x₂), Φ₂₂ = 2x₁x₂(x₃³/3 − h²x₃).
* *The alternative rule* finds 18 states for a ring, 12 for two opposite traction-free faces (with or without one
  further face between them), and none otherwise. The additional ones are face data of fields ½(∇w + ∇wᵀ) with
  different polynomials w on different faces. A w that interpolates between them is not a polynomial, so they are
  not face data of polynomial null sextuples; but the free nodal values can nearly reproduce them, they make the
  system for the stress field badly conditioned, and in the experiments of Section 8 they do not change the bound
  for two of three loads and improve it on coarse grids for the third. This is why they are not used by default.
* *Why six functions.* With the three functions of Maxwell alone far fewer fields change Φ without changing the
  stress, and zero nodal values on two opposite traction-free faces then exclude stress fields that the problem
  needs. Section 8 shows the effect: for a twisted end face the bound with Maxwell's functions and the default
  rule is useless.

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

Both bounds are quadratic in their unknowns. The candidate is the finite-element solution in Q_p (a sparse direct
solve, or conjugate gradients on large grids). The stress field is a maximiser of LB over the free nodal values and
the parameters P_k. That system is singular (Section 3.5). It is solved by an eigen-decomposition when it has fewer
than 3000 unknowns and otherwise by conjugate gradients with the diagonal as preconditioner, after adding 10⁻¹⁰
times the diagonal. (A multiple of the identity does not work: the diagonal entries differ by ten orders of
magnitude between the types of nodal values.) Neither field is trusted; the nodal values are then taken as the
exact binary numbers they are.

## 6. Arithmetic and trusted base

As in two dimensions: rigorous code uses +, −, × of IEEE-754 doubles with outward rounding, comparisons, and
exact rational arithmetic (the coefficients of the two quadratic polynomials, the grid, D, S, b, the polynomials
Ψ). No linear solver and no library function enters a certified number. Trusted: IEEE-754 arithmetic; exact
rational arithmetic; this note.
