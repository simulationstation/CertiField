# Scalar fields in two and three dimensions: the certified bounds

This note states and proves what `ScalarFieldND` certifies, and describes how the two bounds are computed. It is the
analytic reduction behind `Reference/ScalarND/scalarnd.py` (the Python reference implementation) and
`Kernel/ScalarFieldND.wl` (the Wolfram Language port).

**Status.** The statements below marked *Theorem*, *Proposition* and *Lemma* are proved here; the proofs have been
checked by their author only. The underlying inequality is classical (calibrations, the complementary-energy
principle; see [References](#references)) and no novelty is claimed for it. Numbers quoted as *certified* come from
outward-rounded interval arithmetic on exact inputs. Numbers quoted as *reference* are ordinary floating-point values
(series, quadrature) used only for comparison. The numbers of Sections 7 and 8.1 come from the Python reference,
those of Section 8.2 from the Wolfram Language port. The port has been executed and tested
(`Tests/ScalarFieldND.wlt`). It covers two and three dimensions (one dimension has its own model, `ScalarField1D`),
agrees with the reference on identical inputs, and recomputes every certificate stored by the reference
(Section 7).

## 1. The problem

* Ω = (a₁,b₁) × … × (aₙ,bₙ) ⊂ ℝⁿ, n = 1, 2, 3, divided into N₁ × … × Nₙ congruent boxes (*cells*) with edge lengths
  hᵢ = (bᵢ − aᵢ)/Nᵢ.
* Each of the 2n faces of Ω is either a *Dirichlet face* or a *natural face*. Γ_D is the union of the Dirichlet faces.
* k > 0 is the stiffness and W : ℝ → ℝ the potential (continuous; C² wherever second derivatives are used).
* g is the Dirichlet datum on Γ_D. In the implementation it is **the piecewise multilinear interpolant of the given
  nodal values on the boundary grid**, so it coincides with the user's function exactly when that function is
  multilinear on each boundary cell face (constants and x, y, z, xy, … are; a general function is replaced by its
  interpolant, and the certificate says so).

Admissible class and energy:

    𝒜 = { u ∈ H¹(Ω) : u = g on Γ_D (trace) },        E[u] = ∫_Ω k |∇u|² + W(u) dx.

Nothing is prescribed on the natural faces. No symmetry, regularity or bound on u is assumed.

## 2. The lower bound

**Fields.** Let S(x,t) = Σₘ aₘ(x) βₘ(t) (finite sum), where

* each aₘ : Ω̄ → ℝⁿ is a polynomial vector field on every cell, its normal component aₘ·ν is single-valued on every
  interior cell face, and aₘ·ν = 0 on the natural faces;
* each βₘ : ℝ → ℝ is globally Lipschitz and C¹ outside a finite set Z.

For x in the interior of a cell and t ∉ Z put

    Φ(x,t) = div_x S(x,t) + |∂_t S(x,t)|² / (4k) − W(t),          M(x) = sup_{t ∉ Z} Φ(x,t).

**Theorem 1.** If ∫_Ω M⁺ dx < ∞, then for every u ∈ 𝒜 the negative part of W(u) is integrable and

    E[u]  ≥  ∫_{Γ_D} S(x, g(x)) · ν dℋⁿ⁻¹  −  ∫_Ω M(x) dx.

The right-hand side does not depend on u. The field S is arbitrary: a poor S gives a poor bound, never a wrong one.

*Proof.* Fix u ∈ 𝒜 and put F(x) = S(x,u(x)) = Σₘ aₘ(x) βₘ(u(x)).

1. *Chain rule.* For a Lipschitz function β that is C¹ off a finite set Z and u ∈ H¹, β∘u ∈ H¹ and ∇(β∘u) = β′(u)∇u
   almost everywhere, both sides being 0 almost everywhere on {u ∈ Z} (there ∇u = 0 a.e.). Hence on each cell c,
   F ∈ H¹(c; ℝⁿ) and, almost everywhere in c,

       div F = (div_x S)(x,u) + ∂_t S(x,u) · ∇u   on {u ∉ Z},        div F = (div_x S)(x,u)   on {u ∈ Z}.

2. *Pointwise inequality.* On {u ∉ Z}, Young's inequality q·p ≤ |q|²/(4k) + k|p|² gives

       div F ≤ (div_x S)(x,u) + |∂_t S(x,u)|²/(4k) + k|∇u|² = Φ(x,u) + W(u) + k|∇u|² ≤ M(x) + W(u) + k|∇u|².

   On {u ∈ Z}: for t₀ ∈ Z and t → t₀ with t ∉ Z, Φ(x,t) ≥ (div_x S)(x,t) − W(t) → (div_x S)(x,t₀) − W(t₀) because the
   βₘ and W are continuous. So M(x) ≥ (div_x S)(x,t₀) − W(t₀), and again div F ≤ M + W(u) + k|∇u|². Thus

       div F ≤ M + W(u) + k|∇u|²        almost everywhere in Ω.                                        (∗)

3. *Divergence theorem.* On a cell c, ∫_c div F = ∫_{∂c} F·ν_c with the trace F = Σₘ aₘ βₘ(u|_{∂c}); the trace commutes
   with βₘ because βₘ is Lipschitz. Sum over the cells. Each interior face occurs twice with opposite normals, and on
   it both aₘ·ν and the trace of u are single-valued (u ∈ H¹(Ω)), so these contributions cancel. On natural faces
   aₘ·ν = 0. On Γ_D the trace of u is g. Hence Σ_c ∫_c div F = ∫_{Γ_D} S(x,g)·ν.

4. *Integration.* M is bounded below on each cell by the continuous function Φ(·,t), any fixed t ∉ Z, so ∫M is
   well defined in (−∞, ∞). From (∗), W(u) ≥ div F − M⁺ − k|∇u|², which is integrable, so W(u)⁻ ∈ L¹ and E[u] is well
   defined in (−∞, ∞]. Integrating (∗) and using step 3 gives the claim. ∎

**Remarks.**

* The same proof works for a potential V(x,t) that also depends on x (continuous in t): replace W(t) by V(x,t) in Φ.
* *Special cases.* S = 0 gives E ≥ |Ω| inf W. S = t σ(x) gives M = |σ|²/(4k) + W*(div σ) with the convex conjugate
  W*(d) = sup_t (t d − W(t)), i.e. the classical complementary-energy bound

      E[u] ≥ ∫_{Γ_D} g σ·ν − ∫_Ω |σ|²/(4k) − ∫_Ω W*(div σ),

  with equality for convex W at σ = 2k∇u*, u* the minimiser. In one dimension S = s(t) − c x with
  s′ = 2√(k(W + c)) gives M = 0 and the bound of `ScalarField1D`.
* *Tightness.* For scalar problems that are convex in the gradient, the supremum of bounds of this type over all
  admissible fields equals the infimum of the energy (no duality gap) under the hypotheses of Bouchitté and Fragalà
  [2]. Only the elementary direction, Theorem 1, is used here; the size of the gap actually achieved is always
  reported.
* Theorem 1 is special to scalar u. For vector-valued fields (elasticity, magnetisation, Skyrme fields) fields of this
  kind do not in general close the gap.

## 3. The discrete fields and their rigorous evaluation

The implementation uses

    S(x,t) = t σ₀(x) + Σₘ aₘ(x) bₘ(clamp(t; t_lo, t_hi)),

* σ₀ and the aₘ lowest-order Raviart–Thomas fields on the box grid: the i-th component is affine in xᵢ and constant in
  the other coordinates on each cell, so a field is one number per cell face (its normal component there), the normal
  component is single-valued, and the divergence is constant on each cell;
* bₘ, m = 0,…,K+1, the uniform C¹ quadratic B-splines on [t_lo, t_hi] with K knot intervals of width Δ;
  clamp(t) = min(max(t, t_lo), t_hi). So Z = {t_lo, t_hi}.

On a cell c write sᵢ^±(t) for the normal component of S on its two faces orthogonal to eᵢ,
D_c(t) = Σᵢ (sᵢ⁺(t) − sᵢ⁻(t))/hᵢ for the divergence, and a dot for d/dt.

**Proposition 2 (cells to vertices).** For v ∈ {−,+}ⁿ let

    Φ_c(v,t) = D_c(t) − W(t) + (1/4k) Σᵢ ( ṡᵢ^{vᵢ}(t) )²,          M_c(v) = sup_{t ∉ Z} Φ_c(v,t).

Then ∫_c M dx ≤ |c| 2⁻ⁿ Σ_v M_c(v).

*Proof.* In local coordinates ξ ∈ [0,1]ⁿ, ∂_t Sᵢ(x,t) = (1−ξᵢ) ṡᵢ⁻(t) + ξᵢ ṡᵢ⁺(t). So for fixed t, Φ(x,t) is
D_c(t) − W(t) plus a sum of squares of functions that are affine in one coordinate each: it is convex in each ξᵢ
separately. A function on [0,1]ⁿ that is convex in each variable separately lies below the multilinear interpolant of
its vertex values (apply convexity one coordinate at a time). Hence Φ(x,t) ≤ Σ_v λ_v(ξ) Φ_c(v,t) ≤ Σ_v λ_v(ξ) M_c(v),
where λ_v ≥ 0 are the multilinear weights. Take the supremum over t and integrate, using ∫_c λ_v = |c|/2ⁿ. ∎

**Lemma 3 (supremum over t inside the range).** On one knot interval the ṡᵢ^± are affine in t and D_c is a quadratic
polynomial. For a subinterval J = [α,β] of it, of width w,

    sup_J Φ_c(v,·) ≤ max( Φ_c(v,α), Φ_c(v,β) ) + (w²/8) sup_J ( W″ − D_c″ )⁺.

*Proof.* Φ_c(v,·) = C + P with C = (1/4k)Σ(ṡ)² convex and P = D_c − W of class C². C lies below its chord. Let
r = P − (chord of P), so r(α) = r(β) = 0 and r″ = P″. If r has a positive maximum at τ ∈ (α,β), then r′(τ) = 0, and
expanding r about τ and evaluating at the nearer endpoint e (|e − τ| ≤ w/2) gives 0 = r(τ) + r″(ζ)(e − τ)²/2, hence
r(τ) ≤ (w²/8) sup_J (−P″)⁺. The chord of Φ = C + P is at most its larger end value. ∎

**Lemma 4 (supremum over t outside the range).** For t > t_hi the spline part is constant in t, so

    Φ_c(v,t) = D_c^spl(t_hi) + (1/4k) Σᵢ (σ₀ᵢ^{vᵢ})² + ( t d_c − W(t) ),

where D_c^spl is the divergence of the spline part and d_c = div σ₀ on c. Its supremum over t > t_hi is the first two
terms plus T⁺(d_c), T⁺(d) = sup_{t ≥ t_hi} (t d − W(t)); similarly for t < t_lo with T⁻(d) = sup_{t ≤ t_lo} (t d − W(t)).
T^± are convex functions of d (suprema of affine functions).

*Rigorous evaluation of T⁺ for a polynomial W of even degree m ≥ 2 with leading coefficient a_m > 0.* For |d| ≤ d_max
and t ≥ R = max(1, t_hi, (Σ_{j<m} j|a_j| + d_max)/(m a_m)) one has W′(t) ≥ t^{m−2}(m a_m t − Σ_{j<m} j|a_j|) ≥ |d|, so
t d − W(t) is non-increasing on [R, ∞) and the supremum is over [t_hi, R]. There, on any subinterval of width w,
t d − W(t) ≤ (larger end value) + (w²/8) sup (W″)⁺ (Lemma 3 with C = 0). This gives upper bounds at table nodes d_j.
By convexity, for d between two nodes T⁺(d) is at most the chord through the node bounds, and over an interval
[d₁, d₂] it is at most the larger of the two end bounds.

If the affine part is absent (σ₀ = 0) this reduces to T⁺(0) = − inf_{t ≥ t_hi} W. If W is not of this type the affine
part is not used, and the bound is finite only when W is bounded below outside [t_lo, t_hi].

**Lemma 5 (boundary flux).** On a boundary cell face f ⊂ Γ_D orthogonal to eᵢ, with normal component
t σ_f + p_f(clamp t) and outward sign ε = ±1,

    ∫_f S(x,g)·ν = ε |f| σ_f ḡ_f + ε ∫_f p_f(clamp g(x)) dℋⁿ⁻¹,

where ḡ_f, the mean of g over f, is the mean of its 2ⁿ⁻¹ vertex values because g is multilinear on f. Let f′ ⊂ f be a
box on which g has mean ḡ and max_{f′} |g − ḡ| ≤ δ. The function p̂ = p_f∘clamp is Lipschitz with constant
L₁ = sup |ṗ_f| over the knot intervals met by g(f′), and if g(f′) ⊂ [t_lo, t_hi] it is C¹ there with |p̈_f| ≤ L₂. Then

    (1/|f′|) ∫_{f′} ε p̂(g) ≥ ε p̂(ḡ) − min( L₁ δ , L₂ δ²/2 )

(the second alternative only when g(f′) ⊂ [t_lo, t_hi]), by the Lipschitz bound and by Taylor's formula at ḡ, whose
first-order term integrates to zero. If g is constant on f the integral is ε |f| p̂(g) exactly. The implementation
cuts each face into rⁿ⁻¹ congruent boxes (r = 8), on each of which g is again multilinear.

Only the knot intervals whose interior meets [min g, max g] on f′ matter for L₁ and L₂ (outside [t_lo, t_hi] the
function p̂ is constant). Both implementations take a superset, from enclosures of min g and max g: the Wolfram
Language port includes interval κ when (lower enclosure of t_κ) < max g and (upper enclosure of t_{κ+1}) > min g,
by exact comparisons; the Python reference takes ⌊(g − t_lo)/Δ⌋ − 1 … + 1 in floating point, one interval more on
each side, and so a slightly larger penalty.

**Proposition 6 (the affine part alone, K = 0).** For S = t σ₀, M(x) = |σ₀(x)|²/(4k) + W*(d_c) on c, with
W* = max(T⁻, T⁺) for t_lo = t_hi, and ∫_c |σ₀|² = |c| Σᵢ ( (σᵢ⁻)² + σᵢ⁻σᵢ⁺ + (σᵢ⁺)² ) / 3 exactly.

**The certified lower bound.** Theorem 1, Proposition 2 and Lemmas 3–5 give, for any coefficients,

    inf_𝒜 E ≥ LB(S) := Σ_{f ⊂ Γ_D} (lower bound of Lemma 5) − Σ_c |c| 2⁻ⁿ Σ_v M̄_c(v),

M̄_c(v) = the largest of the subinterval bounds of Lemma 3 (each knot interval cut into Q_s equal parts) and the two
tail values of Lemma 4; for K = 0, Proposition 6 replaces the vertex sum. `lower_bound` evaluates this with
outward-rounded interval arithmetic: sample points, B-spline weights and W at the samples are exact rationals; W″ on
a subinterval is enclosed by interval Horner evaluation.

## 4. The upper bound

**Proposition 7.** Let u_h be the continuous function that is multilinear on every cell (Q1) with nodal values U,
equal to the given nodal values of g at the nodes of Γ_D. Then u_h ∈ 𝒜, and on each cell c

* ∫_c |∂ᵢu_h|² = |c| Σ_{e,e′} d_e d_{e′} Π_{j≠i} m(e_j,e′_j), where d_e are the difference quotients of U along the
  2ⁿ⁻¹ edges of c parallel to eᵢ and m(0,0) = m(1,1) = 1/3, m(0,1) = m(1,0) = 1/6 (∂ᵢu_h is multilinear in the other
  coordinates);
* | ∫_c W(u_h) − |c| W(ū_c) | ≤ (|c|/24) Σᵢ hᵢ² sup_c |W″(u_h)| (∂ᵢu_h)², where ū_c, the value of u_h at the centre of
  c, is the mean of its 2ⁿ vertex values. On c, u_h lies between the smallest and the largest vertex value and
  |∂ᵢu_h| ≤ max_e |d_e|.

*Proof.* The first item is the integral of the square of a multilinear function. The second is the midpoint rule
applied one coordinate at a time, | ∫_{−h/2}^{h/2} φ − h φ(0) | ≤ (h³/24) sup |φ″|, with F = W∘u_h and
∂ᵢᵢF = W″(u_h)(∂ᵢu_h)² because ∂ᵢᵢu_h = 0. ∎

`upper_bound` returns an enclosure [lo, hi] of E[u_h]; hi is the certified upper bound for inf E. Any nodal values
give a valid bound. The candidate only affects how tight it is. When no candidate is given, both implementations
minimise the Q1 energy with two-point Gauss quadrature for the potential (L-BFGS in the reference; in the port
Newton's method, with conjugate gradients truncated at negative curvature when W is not convex, so that the
iteration descends to a minimum and not to a saddle point); the result is a candidate like any other.

## 5. Choosing the field (heuristic)

Any field gives a valid bound, so this part is ordinary floating point and is never trusted.

**Convex mode (K = 0).** With the vertex rule for |σ₀|² the functional σ₀ ↦ flux − Σ_c |c| ( Σᵢ((σᵢ⁻)² + (σᵢ⁺)²)/(8k) +
W*(d_c) ) is stationary exactly when σ₀·eᵢ = 2k (p_R − p_L)/hᵢ on interior faces, σ₀·ν = 4k (g_f − p_c)/hᵢ on
Dirichlet faces and σ₀·ν = 0 on natural faces, with W′(p_c) = d_c: the cell-centred finite-volume scheme for
2kΔp = W′(p). (This is the known equivalence of lowest-order mixed elements with trapezoidal quadrature and
cell-centred differences [6].) The implementation solves this finite-volume problem (by minimising its energy,
starting from the cell means of the candidate), forms these fluxes, and evaluates the bound with the exact
Proposition 6. The field is not optimised further. In the examples the bound is tight to second order in the mesh
width when W is convex; for non-convex W it is still valid but is the bound for the convexified potential, since
W* is also the conjugate of the convex envelope of W.

**Calibration mode (K > 0, spline part alone).** The right-hand side of Theorem 1 is concave in the coefficients:
the flux is linear in them and M is a supremum of convex functions of them. A smoothed copy (the supremum over t
replaced by τ·log-mean-exp over the samples, with the same curvature allowance as Lemma 3 so that nothing can hide
between samples, and the flux by a Gauss rule) is maximised by L-BFGS for a decreasing sequence of τ, coarse to fine:
the field is embedded exactly into the space with twice as many cells per direction (Raviart–Thomas spaces are
nested) or twice as many knot intervals (B-spline subdivision), and the field with the best certified value is kept
at every stage. For a tight bound the range [t_lo, t_hi] should contain the values of the datum and the wells of W
that competitors can reach; validity does not depend on it.

**Automatic choices in the Wolfram Language port.** Convex mode is used when W is convex and of degree at least
two (W″ ≥ 0 on ℝ, decided exactly), calibration mode otherwise. Unless given, the range is the hull of the Dirichlet data and of the
critical points of W, rounded outward to a binary fraction of its length: when W has even degree and a positive
leading coefficient it decreases to the left of its first critical point and increases to the right of its last
one, so truncating a competitor to that range does not increase its energy. The smoothing parameters τ are
multiplied by max(1, range of W on [t_lo, t_hi]), and the levels halve the grid up to three times and K down to 8
before doubling them alternately. All of this is heuristic.

Optimising the affine and spline parts together is allowed by the theory and by the verifier but converges slowly
(the two parts overlap inside the range), so the two modes are kept separate.

## 6. Arithmetic and trusted base

Rigorous code uses only +, −, ×, ÷ of IEEE-754 doubles, comparisons, and exact rational arithmetic. In the Python
reference every operation is widened by one unit in the last place in each direction (`numpy.nextafter`); sums of
many terms use `math.fsum` widened by two units; rationals are converted to enclosing doubles by exact comparison;
stored doubles (nodal values, field coefficients) are taken as the exact binary numbers they are. No elementary
function from a math library enters a certified number. A floating-point root finder is used in one place, to
choose where the tables of Lemma 4 are refined; the bound holds for any choice.

The Wolfram Language port uses the paclet's outward-rounded array arithmetic (every operation widened by 2⁻⁴⁹
relative plus 10⁻³⁰⁰ in each direction) and exact rationals; sums over cells are exact rational sums. Two points
are specific to the Wolfram kernel. An overflowing machine operation does not return an infinity: the kernel
repeats it in software arithmetic and returns a number that is not a machine number. The enclosure argument does
not cover that, so the result of every interval operation is tested and the evaluation stops as inconclusive if
it happens. Underflow to subnormal numbers is harmless: gradual underflow is exact to 2⁻¹⁰⁷⁵, well inside the
absolute widening.

Trusted: IEEE-754 arithmetic for the four operations; the proofs in this note; for the Wolfram Language port also
exact rational arithmetic of the kernel, and that the kernel returns either the IEEE result of a machine operation
or a number that is not a machine number.

## 7. Checks

The first four were performed with the Python reference (`Reference/ScalarND`), on Linux and again on Windows; the
last one with the Wolfram Language port.

* **Known answers enclosed.** In every case with an independent reference value the certified interval contains it
  (table below).
* **Falsification** (`check_theorem.py`). For problems with known minima, the certified value LB(S) was evaluated for
  optimised fields, for several hundred random perturbations of them, for fields produced by an optimiser that is
  allowed to cheat between samples (no curvature allowance, one sample per knot interval), and for fields mixing the
  affine and spline parts, in one, two and three dimensions. None exceeded the reference; the closest came within
  1.1×10⁻⁴ below it.
* **Exact cross-check** (`check_exact.py`). An independent evaluator in exact rational arithmetic (no floating point,
  no shared code) recomputes the formula of Section 3 on small grids. The interval value never exceeds the exact one
  and agrees with it to 10⁻¹⁴ in calibration mode and to 10⁻⁶ or better in convex mode (the difference there is the
  chord of the conjugate table). The exact energy of the Q1 function lies inside the computed enclosure.
* **Stored certificates.** Each example stores the problem, the candidate and the field as exact doubles;
  `verify_certificate` recomputes both bounds from them and compares them with the stored claims. A claim is the
  computed bound weakened by one part in 10⁹ and rounded outward to ten digits, so the comparison does not depend on
  the last bits of the root finder mentioned in Section 6, which differ between platforms (perturbing its roots by
  a relative 10⁻⁶ moved the bounds of the examples by less than 10⁻¹² relative).
* **Second implementation** (`Tests/ScalarFieldND.wlt`, `Tests/ScalarFieldNDSlow.wlt`). The Wolfram Language
  port follows the reference function by function, but it is a second implementation in another language, with
  different rounding, exact summation and its own solvers and optimiser; a mistake in the evaluation of the
  formula of Section 3 would have to be made twice to go unnoticed. On the reference's test vectors it reproduces the competitor energies and, where the datum is constant
  on each face, the lower bounds to 10⁻¹³ (to 10⁻⁶ for fields with an affine part, whose conjugate tables the two
  implementations place differently). Where the datum varies on the faces its flux bound is tighter than the
  reference's (Lemma 5, last paragraph); there it is tested from both sides: it is at least the reference's value
  and at most the exact flux, which the test computes in rational arithmetic by integrating the clamped spline in
  closed form. The port recomputes all eleven stored certificates, and both of its bounds are at least as strong as
  the stored claims. The known answers (the sine-series values and the wall value) lie inside the intervals it
  certifies with its own candidates and fields.

## 8. Results

E = ∫ ½|∇u|² + W(u) throughout (k = ½).

### 8.1 Python reference (one core)

| Potential W | Domain and data | Mode | Grid | K | Certified: lower ≤ inf E ≤ upper | Gap | Time | Reference value |
|---|---|---|---|---|---|---|---|---|
| u² − u | u = 0 on ∂(0,1)² | convex | 64×64 | – | -0.01602884 ≤ inf E ≤ -0.01601929 | 0.060 % | 1 s | -0.01602506 |
| u² − u | the same | convex | 128×128 | – | -0.01602601 ≤ inf E ≤ -0.01602361 | 0.015 % | 2 s | -0.01602506 |
| u² − u | u = 0 on ∂(0,1)³ | convex | 32×32×32 | – | -0.00951680 ≤ inf E ≤ -0.00947745 | 0.416 % | 2 s | -0.00949862 |
| u² − u | the same | convex | 48×48×48 | – | -0.00950668 ≤ inf E ≤ -0.00948922 | 0.184 % | 11 s | -0.00949862 |
| u⁴/4 + u²/2 − u | u = 0 on ∂(0,1)² | convex | 64×64 | – | -0.01676427 ≤ inf E ≤ -0.01675442 | 0.059 % | 1 s | none known |
| u⁴/4 + u²/2 − u | u = 0 on ∂(0,1)³ | convex | 32×32×32 | – | -0.00980025 ≤ inf E ≤ -0.00976023 | 0.410 % | 2 s | none known |
| (u² − 1)²/4 | (−3,3)×(0,1): u = −1 and +1 on the ends, sides free | calibration | 96×4 | 32 | 0.944210 ≤ inf E ≤ 0.945112 | 0.096 % | 37 s | 0.94496743 |
| (u² − 1)²/4 | (−3,3)×(0,1)²: the same | calibration | 48×2×2 | 32 | 0.944198 ≤ inf E ≤ 0.945551 | 0.143 % | 36 s | 0.94496743 |
| 4(u² − 1)² | u = x on ∂(−1,1)² | calibration | 64×64 | 32 | 8.275401 ≤ inf E ≤ 8.291702 | 0.197 % | 7 min | none known |
| 4(u² − 1)² | u = xy on ∂(−1,1)² | calibration | 64×64 | 32 | 10.679217 ≤ inf E ≤ 10.952917 | 2.499 % | 7 min | none known |
| (u² − 1)² | u = x on ∂(−1,1)³ | calibration | 32×32×32 | 16 | 8.014875 ≤ inf E ≤ 8.027319 | 0.156 % | 21 min | none known |

Times are for the whole run on one core (Python/NumPy): candidate, field, and both certified evaluations. Recomputing
both bounds from a stored certificate takes 1 to 10 seconds. Reference values are floating-point evaluations of a sine
series (first four rows) and of the one-dimensional first integral (the two walls, whose minimum is the cross-section
times the 1D value). Every reference value lies inside its certified interval.

Coarser levels of the calibration runs, for orientation: the 2D ramp at 32×32, K = 16 has a certified gap of 0.74 %
after 27 s; the 3D ramp at 16×16×16, K = 16 has 0.45 % after under 4 minutes.

The saddle example has two symmetric minimisers. Its gap closes slowly as K grows (5.1 % → 3.4 % → 2.5 % along the
last three levels), which is the behaviour to expect when the minimiser is not unique.

### 8.2 Wolfram Language port

The same problems through `CertifiedVariationalBound`, with the candidate and the field computed in the Wolfram
Language (no input from the Python). The grids, the t-range [−1, 1], the numbers of knot intervals and the levels
are those of Section 8.1; they are the automatic choices of the port except for the two walls (their levels are
given as in the reference) and the saddle (whose candidate is started from xy + 0.3(1 − x²)(1 − y²), as in the
reference, to leave the symmetric saddle point).

| Potential W | Domain and data | Mode | Grid | K | Certified: lower ≤ inf E ≤ upper | Gap | Time | Reference value |
|---|---|---|---|---|---|---|---|---|
| u² − u | u = 0 on ∂(0,1)² | convex | 32×32 | – | -0.01604019 ≤ inf E ≤ -0.01600195 | 0.239 % | 1.1 s | -0.01602506 |
| u² − u | the same | convex | 64×64 | – | -0.01602884 ≤ inf E ≤ -0.01601929 | 0.060 % | 1.2 s | -0.01602506 |
| u² − u | the same | convex | 128×128 | – | -0.01602601 ≤ inf E ≤ -0.01602361 | 0.015 % | 1.7 s | -0.01602506 |
| u² − u | u = 0 on ∂(0,1)³ | convex | 16×16×16 | – | -0.00957197 ≤ inf E ≤ -0.00941375 | 1.681 % | 1.4 s | -0.00949862 |
| u² − u | the same | convex | 32×32×32 | – | -0.00951680 ≤ inf E ≤ -0.00947745 | 0.416 % | 2.4 s | -0.00949862 |
| u² − u | the same | convex | 48×48×48 | – | -0.00950668 ≤ inf E ≤ -0.00948922 | 0.184 % | 6.1 s | -0.00949862 |
| u⁴/4 + u²/2 − u | u = 0 on ∂(0,1)² | convex | 64×64 | – | -0.01676427 ≤ inf E ≤ -0.01675442 | 0.059 % | 1.4 s | none known |
| u⁴/4 + u²/2 − u | u = 0 on ∂(0,1)³ | convex | 32×32×32 | – | -0.00980025 ≤ inf E ≤ -0.00976023 | 0.410 % | 2.7 s | none known |
| (u² − 1)²/4 | (−3,3)×(0,1): u = −1 and +1 on the ends, sides free | calibration | 96×4 | 32 | 0.944210 ≤ inf E ≤ 0.945112 | 0.096 % | 39 s | 0.94496743 |
| (u² − 1)²/4 | (−3,3)×(0,1)²: the same | calibration | 48×2×2 | 32 | 0.944198 ≤ inf E ≤ 0.945551 | 0.143 % | 54 s | 0.94496743 |
| 4(u² − 1)² | u = x on ∂(−1,1)² | calibration | 32×32 | 16 | 8.245472 ≤ inf E ≤ 8.305251 | 0.720 % | 43 s | none known |
| 4(u² − 1)² | the same | calibration | 64×64 | 32 | 8.275406 ≤ inf E ≤ 8.291702 | 0.197 % | 6 min | none known |
| 4(u² − 1)² | u = xy on ∂(−1,1)² | calibration | 64×64 | 32 | 10.679179 ≤ inf E ≤ 10.952917 | 2.500 % | 6 min | none known |
| (u² − 1)² | u = x on ∂(−1,1)³ | calibration | 16×16×16 | 16 | 8.005335 ≤ inf E ≤ 8.040144 | 0.433 % | 7 min | none known |
| (u² − 1)² | the same | calibration | 32×32×32 | 16 | 8.015122 ≤ inf E ≤ 8.027319 | 0.152 % | 18 min | none known |

Times are for the whole certificate (candidate, field, and both certified evaluations) with Wolfram Engine 15.0 on
a 16-core x86-64 processor, each example in a fresh kernel; the kernel's array operations use several threads, so
these times are not one-core times like those of Section 8.1. `VerifyCertificate` recomputes both bounds from a
certificate in 1 to 15 seconds. Every reference value lies inside its certified interval.

The upper bounds agree with those of Section 8.1 to all digits shown, and so do the lower bounds of the convex
rows (the fields solve the same finite-volume problem) and of the two walls. In the other calibration rows the
fields come from two different optimisers and the flux bounds differ (Lemma 5, last paragraph), so the lower bounds
agree less closely: 8.275406 against 8.275401 for the 2D ramp, 10.679179 against 10.679217 for the saddle, and
8.015122 against 8.014875 for the 3D ramp. On the reference's own fields the port's lower bounds equal the
reference's where the datum is constant on each face and are slightly higher where it varies; for the stored field
of the 3D ramp the port also gives 8.015122.

## 9. What is not covered

* **Not claimed:** that the candidate is the minimiser, or that the minimiser is unique. The gap bounds how far below
  the candidate any admissible field can be.
* Potentials that depend on x (for example a source term f(x)u). Theorem 1 covers them; the evaluation in Section 3
  assumes W is the same on a cell. A source term can be reduced to cell-wise constants with the Poincaré inequality
  on each cell at second-order cost; this is not implemented.
* Potentials without superlinear growth in convex mode, in particular the pure Poisson problem W(u) = −f u, which
  needs exactly equilibrated fluxes.
* Domains other than boxes; curved boundary data other than through their multilinear interpolant.
* The Python reference takes polynomial W. The method needs W at sample points, bounds for W″ on intervals, and the
  tail suprema of Lemma 4.
* Vector-valued fields.
* Calibration mode is costly in three dimensions and converges slowly in the number of knot intervals when the
  minimiser is not unique (see the saddle example).

## References

1. G. Alberti, G. Bouchitté, G. Dal Maso, *The calibration method for the Mumford–Shah functional and
   free-discontinuity problems*, Calc. Var. Partial Differential Equations 16 (2003), 299–333.
2. G. Bouchitté, I. Fragalà, *A duality theory for non-convex problems in the calculus of variations*, Arch. Ration.
   Mech. Anal. (2018), doi:10.1007/s00205-018-1219-3; announced in C. R. Math. Acad. Sci. Paris 353 (2015), 375–379.
3. T. Pock, D. Cremers, H. Bischof, A. Chambolle, *Global solutions of variational models with convex
   regularization*, SIAM J. Imaging Sci. 3 (2010), doi:10.1137/090757617. (Functional lifting: the numerical
   counterpart of the same duality, without certification.)
4. W. Prager, J. L. Synge, *Approximations in elasticity based on the concept of function space*, Quart. Appl. Math.
   5 (1947), 241–269. (Complementary energy, the hypercircle.)
5. D. Gilbarg, N. S. Trudinger, *Elliptic Partial Differential Equations of Second Order*, Section 7.4 (chain rule
   for piecewise C¹ functions of Sobolev functions).
6. T. F. Russell, M. F. Wheeler, *Finite element and finite difference methods for continuous flows in porous
   media*, in The Mathematics of Reservoir Simulation, SIAM (1983). (Mixed elements with quadrature as cell-centred
   differences; cited from memory, attribution confirmed through later papers of Wheeler and Yotov.)

A brief web search (November 2026) found the theory above and uncertified numerical implementations of it, and did
not find a packaged implementation that certifies the bound with rigorous arithmetic. That search is not a
literature review.
