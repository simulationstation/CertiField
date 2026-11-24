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
  the problem is refused; this did not happen in any arrangement of the faces tried (Section 7). The space
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
* *The alternative rule* finds 18 states for a ring, 12 for two opposite traction-free faces, 6 when a third face
  joins them, and none in the other arrangements tried. The additional ones are face data of fields ½(∇w + ∇wᵀ) with
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

## 7. Checks

All with the Python reference (`Reference/Elasticity`, files with `3d` in their names).

* **Exact cross-check** (`check_exact3d.py`, log in `logs3d/`). For twelve problems on grids of two to eight cells
  (all faces Dirichlet; one, two adjacent, two opposite, three, four and five traction-free faces; with and
  without body force; isotropic and fully anisotropic stiffness; p = 1 and 2; six, Maxwell's three and Morera's
  three stress functions; both rules for the states), with the computed fields and with random ones, both fields
  are rebuilt cell by cell as polynomials in three variables in rational arithmetic. The stress formulas are typed
  in from Section 3.2, not taken from the table that the reference generates from the ε-formula. Verified exactly
  in all 24 cases: div σ + b = 0 in every cell; the traction is the same from both sides of every interior cell
  face; it vanishes identically on the traction-free faces; the candidate is continuous. The exact values of both
  bounds lie inside the computed enclosures, which are 10⁻¹¹ wide or less. The identity
  E[u_h] − LB(σ) = ½∫(e − Ss)ᵀD(e − Ss) holds as an equality of rational numbers.
* **Known minima** (`check_theorem3d.py`, part A). Homogeneous strain with all faces Dirichlet (six, three and three
  stress functions); uniaxial tension of a bar with four traction-free faces along each of the three axes; a
  quadratic displacement with constant body force; pure bending of a bar with four traction-free faces. Data,
  grids and compliance are dyadic, so the stored nodal values are exactly the data. The exact rational minimum
  lies inside the certified interval in all ten cases, and both bounds equal it to 10⁻⁹ (the exact fields are in
  the spaces).
* **Falsification** (part B). The block of Section 8 has the same datum on every grid, so every lower bound is
  compared with every upper bound: optimal and perturbed fields on seven grids, with six and with three stress
  functions and both rules for the states (63 stress fields, 84 candidates). Four further problems (body force, anisotropic stiffness, five
  and two traction-free faces, a twisted and sheared end face) are tested within each grid. No pair has
  LB(σ) > E[u].
* **States** (part C; this part of the script had not finished when the work was interrupted, see the log). For
  the eight arrangements tried by hand, with and without body force, the particular sextuple Ψ_b exists, a ring
  has six states with independent resultants, and the alternative rule finds 12 states for two opposite
  traction-free faces, 6 when a third face joins them, 18 for a ring and none otherwise.
* **Stored certificates** (`verify_all3d.py`) recompute both bounds from the stored fields.

## 8. Results (Python reference)

Certified numbers; times on one desktop machine, several of them measured while other jobs were running.

**E-T4.** Unit cube, λ = 3/2, μ = 1, clamped on the face x₁ = 0, displaced by (0.1, 0, 0) on the face x₁ = 1, no
traction on the four other faces; p = 2, six stress functions, default rule for the states (r = 6):

| Grid | Certified: lower ≤ inf E ≤ upper | Gap | Time |
|---|---|---|---|
| 2×2×2 | 0.01356377 ≤ inf E ≤ 0.01400734 | 3.167 % | 9 s |
| 3×3×3 | 0.01368516 ≤ inf E ≤ 0.01392955 | 1.754 % | 10 s |
| 4×4×4 | 0.01373380 ≤ inf E ≤ 0.01389667 | 1.172 % | 10 s |
| 5×5×5 | 0.01375984 ≤ inf E ≤ 0.01387910 | 0.859 % | 10 s |
| 6×6×6 | 0.01377585 ≤ inf E ≤ 0.01386834 | 0.667 % | 14 s |

The gap is below 1 % from 5×5×5 on. It falls roughly like h^1.4, as in two dimensions: the edges where the clamped
faces meet the traction-free faces carry stress singularities. Of the time, about 4 s is the exact local matrix of
the stress functions (384 × 384 rationals) and 2 s the rigorous evaluation; the rest is the two heuristic solves
(conjugate gradients for the stress field: 954, 1422 and 2288 iterations on the last three grids).

**Other loads on the same cube** (4×4×4, p = 2): end face displaced sideways by (0, 0.1, 0),
0.00331920 ≤ inf E ≤ 0.00338967 (2.08 %); end face rotated about the axis,
0.00072267 ≤ inf E ≤ 0.00073014 (1.02 %).

**Other arrangements.** Pure bending of a bar with four traction-free faces (2×2×2, exact minimum 1/32):
0.03124999996 ≤ inf E ≤ 0.03125000004. A column 1 × 1 × 2 clamped at the bottom under its own weight, five
traction-free faces (3×3×6): −0.49544773 ≤ inf E ≤ −0.48957355 (1.20 %). A plate 1 × ½ × 1 clamped on its four
edges under its own weight, two traction-free faces: gap 21.0 % (2×1×2), 8.85 % (4×2×4), 3.93 % (6×3×6).

**Six functions against three, and the two rules for the states** (lower bounds; the upper bound is the same):

| Problem | Grid | Six functions, default rule | Six functions, all states | Maxwell's three, all states | Maxwell's three, default rule |
|---|---|---|---|---|---|
| cube, pulled | 3×3×3 | 0.01368516 | 0.01368516 | 0.01368510 | 0.01368510 |
| cube, end displaced sideways | 3×3×3 | 0.00329906 | 0.00329906 | 0.00329680 | 0.00329677 |
| cube, end rotated | 2×2×2 | 0.00070966 | 0.00071860 | 0.00070660 | 0.00015452 |
| cube, end rotated | 3×3×3 | 0.00071902 | 0.00072146 | 0.00071746 | 0.00015615 |
| plate under its weight | 4×2×4 | −0.01212271 | −0.01206074 | −0.01451176 | |
| plate under its weight | 6×3×6 | −0.01178063 | −0.01178063 | −0.01422838 | |

(Upper bounds: 0.01392955, 0.00341047, 0.00074705, 0.00073401, −0.01113753, −0.01133525.) With six functions the
additional states of the alternative rule matter only on coarse grids. With Maxwell's three functions the bound
for the plate stays 25 % away under refinement, and with the default rule the bound for the rotated end face is
useless: these are the cases that Section 3.5 predicts. All four columns are valid lower bounds.

## 9. What is not covered

* **Not claimed:** that the candidate is the minimiser. The gap bounds E[u_h] − inf E and, by Section 2, the
  energy-norm errors of the candidate and of the stress field.
* **Not claimed:** that the bound converges to the infimum under refinement. Section 3.5 gives the reasons to
  expect it and is a sketch; the tables of Section 8 are the evidence.
* A Wolfram Language model. `testvectors3d/` holds inputs, expected outputs and the polynomial states for a port.
* Conditions that prescribe only some components of u on a face (rollers), prescribed non-zero tractions, body
  forces that are not constant, stiffness that varies in space, domains other than boxes.
* Edges where a Dirichlet face meets a traction-free face carry stress singularities. Both bounds remain valid,
  but the gap closes slowly under refinement, and the grids are uniform.
* Nearly incompressible materials: the bounds remain valid; the candidate locks and the gap grows.
* Cost. The stress field has 48 unknowns per node. The dual solve by conjugate gradients needs a number of
  iterations that grows with the grid; beyond about 8³ cells a better preconditioner or nested iteration would be
  needed.

## References

1. W. Prager, J. L. Synge, *Approximations in elasticity based on the concept of function space*, Quart. Appl. Math.
   5 (1947), 241–269.
2. M. E. Gurtin, *The linear theory of elasticity*, in Handbuch der Physik VIa/2, Springer, 1972. (The stress
   functions of Maxwell, Morera and Beltrami; the standard reference. Not consulted directly: cited through
   the secondary sources found.)
3. M. E. Gurtin, *A generalization of the Beltrami stress functions in continuum mechanics*, Arch. Rational Mech.
   Anal. 13 (1963), 321–329. (Completeness for self-equilibrated stress fields.)
4. R. Rostamian, *The completeness of Maxwell's stress function representation*, J. Elasticity 9 (1979), 349–356.
5. V. Volterra, *Sur l'équilibre des corps élastiques multiplement connexes*, Ann. Sci. École Norm. Sup. (3) 24
   (1907), 401–517. (Many-valued displacements in multiply connected bodies: the six parameters of Section 3.5.)
6. M. Křížek, *An equilibrium finite element method in three-dimensional elasticity*, Apl. Mat. 27 (1982), 46–75.
   (An equilibrium element on tetrahedra for the dual problem.)
7. J. P. Moitinho de Almeida, E. A. W. Maunder, *Equilibrium Finite Element Formulations*, Wiley, 2017.
8. F. Key, L. Freinberger, *A spline-based stress function approach for the principle of minimum complementary
   energy*, Comput. Methods Appl. Mech. Engrg. 448 (2025) 118492, arXiv:2506.19534. (Smooth stress functions in two
   dimensions.)

A brief web search (November 2026) confirmed the bibliographic data of items 3 to 6 and 8 (items 1 and 7 were
confirmed for the two-dimensional note). Completeness of the representations (items 2 to 4) concerns how good the
bound can be and is not used here. The search found equilibrium finite elements for three dimensions built
directly in stress spaces (items 6, 7), and stress functions in smooth spline spaces in two dimensions (item 8).
It did not find an implementation with Beltrami stress functions in a C¹ tensor-product space, nor one that
enforces static admissibility exactly in floating-point data and evaluates the bound with rigorous arithmetic.
That search is not a literature review, and the absence of a result is not a claim of novelty.
