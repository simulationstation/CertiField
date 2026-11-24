# Scalar fields in two and three dimensions, part 2: position-dependent potentials, constraints, periodicity

This note extends `docs/SCALAR_ND.md` (*part 1*). It keeps the numbering of part 1: Theorem 1, Proposition 2,
Lemmas 3 to 5, Propositions 6 and 7 are the statements there; the statements here are numbered from 8.

**Status (23 November 2026).** The statements and proofs below are complete and were written before any code. They
have been checked by their author only. Everything in them is classical or elementary (references at the end) and
no novelty is claimed. The implementation is `Kernel/ScalarFieldND2.wl`, reached through the model `ScalarFieldND`
when one of the keys "Source", "Mean", "PeriodicDirections", "ValueRange", "DataTerm" is given; its tests are
`Tests/ScalarFieldND2.wlt`. Problems without these keys are handled exactly as in version 0.3.0.

## 1. The extended problem

* Ω, the cells and the edge lengths hᵢ are as in part 1; Lᵢ = bᵢ − aᵢ.
* Each direction i is either *periodic* (the two faces orthogonal to eᵢ are identified) or its two faces are each a
  Dirichlet face or a natural face.
* The stiffness is a constant symmetric positive definite matrix K in the statements. The implementation takes
  K = k·I as in part 1.
* The potential depends on position: on a cell c,

      W(x, t) = W_c(t) + a(x) t ,

  where a is multilinear on c (it need not be continuous across cell faces), and W_c is finite and continuous on a
  closed interval 𝕀 (the *value range*, the same for all cells; 𝕀 = ℝ or 𝕀 = [τ₋, τ₊]) and W_c = +∞ outside 𝕀.
* Optionally the mean of u is prescribed: ∫_Ω u = m|Ω|.

Admissible class and energy:

    𝒜 = { u ∈ H¹(Ω) : u = g on Γ_D, equal traces on opposite periodic faces [, ∫_Ω u = m|Ω| ] },
    E[u] = ∫_Ω ∇u·K∇u + W(x, u) dx,      E[u] := +∞ if u ∉ 𝕀 on a set of positive measure.

Kinds of W_c the implementation is to take (the theory does not depend on this list):

| Kind | W_c(t) | 𝕀 |
|---|---|---|
| polynomial | Σ wⱼ tʲ, rational coefficients | ℝ, or an interval [τ₋, τ₊] |
| polynomial with entropy terms | P(t) + c₁ (t − τ₋) ln(t − τ₋) + c₂ (τ₊ − t) ln(τ₊ − t), c₁, c₂ ≥ 0 rational, P a polynomial | [τ₋, τ₊] |
| table | a polynomial plus, for each cell, its own continuous piecewise linear function on a uniform grid of [τ₋, τ₊] | [τ₋, τ₊] |
| linear | w₀ + w₁ t (no growth; Section 6) | ℝ |

The regular-solution free energy t ln t + (1 − t) ln(1 − t) + χ t(1 − t) on [0, 1] is the second kind with
c₁ = c₂ = 1.

## 2. The inequality

**Theorem 8 (Theorem 1 for the extended problem).** Let S(x,t) = Σₘ aₘ(x) βₘ(t) be as in Section 2 of part 1 and,
in addition, let aₘ·ν take the same value at corresponding points of opposite periodic faces. For x in the interior
of a cell and t ∈ 𝕀 ∖ Z put

    Φ(x,t) = div_x S(x,t) + ¼ ∂_t S(x,t)·K⁻¹ ∂_t S(x,t) − W(x,t),          M(x) = sup_{t ∈ 𝕀∖Z} Φ(x,t).

If ∫_Ω M⁺ < ∞, then for every u ∈ 𝒜 (without the constraint on the mean)

    E[u]  ≥  ∫_{Γ_D} S(x, g(x))·ν dℋⁿ⁻¹  −  ∫_Ω M(x) dx.

*Proof.* If u ∉ 𝕀 on a set of positive measure there is nothing to prove. Otherwise u(x) ∈ 𝕀 almost everywhere,
and the proof of Theorem 1 is repeated with three changes.

*Step 2 (pointwise inequality).* For p, q ∈ ℝⁿ, q·p ≤ ¼ q·K⁻¹q + p·Kp, because
|K^{1/2}p − ½K^{−1/2}q|² ≥ 0. With q = ∂_tS(x,u), p = ∇u this gives, on {u ∉ Z},
div F ≤ Φ(x,u) + W(x,u) + ∇u·K∇u ≤ M(x) + W(x,u) + ∇u·K∇u. On {u = t₀}, t₀ ∈ Z: 𝕀 is a non-degenerate interval and
Z is finite, so there are t ∈ 𝕀∖Z with t → t₀; W(x,·) is continuous on 𝕀 and ¼q·K⁻¹q ≥ 0, so
M(x) ≥ (div_xS)(x,t₀) − W(x,t₀) as before. The dependence of W on x plays no role: the inequality is pointwise in x.

*Step 3 (divergence theorem).* The two faces of a periodic pair carry the same values of aₘ·eᵢ and of the trace of
u, with opposite outward normals, so their contributions cancel like those of an interior face.

*Step 4.* M is bounded below on each cell by Φ(·,t) for one fixed t ∈ 𝕀∖Z, which is bounded there because a is. ∎

**Corollary 9 (prescribed mean).** For every λ ∈ ℝ,

    inf { E[u] : u ∈ 𝒜, ∫u = m|Ω| }  ≥  λ m |Ω| + ∫_{Γ_D} S(x,g)·ν − ∫_Ω M_λ ,

where M_λ is formed with W(x,t) − λt in place of W(x,t).

*Proof.* For u with ∫u = m|Ω|, E[u] = E_λ[u] + λ m|Ω|, where E_λ has the potential W − λt; apply Theorem 8 to
E_λ on the larger class without the constraint. ∎

The right-hand side is concave in (S, λ) jointly (the flux is linear, M_λ is a supremum of functions that are
convex in S and affine in λ), so λ can be handed to the optimiser as one more unknown. A constant λ is the term
a(x)t with a ≡ −λ: nothing new has to be evaluated.

## 3. Rigorous evaluation of the lower bound

The discrete fields are those of Section 3 of part 1. With a periodic direction the face values on the two
identified faces are one unknown. Notation of part 1: sᵢ^±, D_c, a dot for d/dt. Let a_v be the value of a at the
vertex v of the cell c (the cell's own value, if a jumps).

**Proposition 10 (cells to vertices).** With

    Φ_c(v,t) = D_c(t) − W_c(t) − a_v t + (1/4k) Σᵢ ( ṡᵢ^{vᵢ}(t) )²,          M_c(v) = sup_{t ∈ 𝕀∖Z} Φ_c(v,t),

∫_c M dx ≤ |c| 2⁻ⁿ Σ_v M_c(v). For a constant positive definite K the last term is ¼ q_v·K⁻¹q_v with
q_v = (ṡᵢ^{vᵢ}(t))ᵢ and the statement is the same.

*Proof.* As for Proposition 2: for fixed t, −a(x)t is affine in each local coordinate ξᵢ, and the quadratic term is
convex in each ξᵢ (only the i-th component of ∂_tS depends on ξᵢ, affinely, and K⁻¹ is positive definite). So
Φ(·,t) is convex in each coordinate separately and lies below the multilinear interpolant of its vertex values,
which are the Φ_c(v,t). ∎

*Lemma 3 is unchanged*: the term −a_v t is affine in t and does not enter the second derivative. For a table, W_c″
vanishes inside every table interval, so the subintervals of Lemma 3 must not straddle a table node (the number of
table intervals divides K·Q_s) and the allowance is (w²/8)(−D_c″)⁺.

*Lemma 4 with a linear term.* For t > t_hi, Φ_c(v,t) is a constant plus t(d_c − a_v) − W_c(t), so the tail value is
T⁺(d_c − a_v) (with the multiplier of Corollary 9, T⁺(d_c − a_v + λ)); likewise T⁻. The conjugate tables of part 1
are used at these shifted arguments. If 𝕀 = [τ₋, τ₊] the spline range is taken to be 𝕀 itself and there are no
tails. *Proposition 6 with a linear term:* for S = tσ₀, M(x) = |σ₀(x)|²/(4k) + W_c*(d_c − a(x)); the first term is
integrated exactly as before, and x ↦ W_c*(d_c − a(x)) is convex in each coordinate (a convex function of a
function that is affine in ξᵢ), so its integral over c is at most |c| 2⁻ⁿ Σ_v W_c*(d_c − a_v).

**Lemma 11 (a subinterval that ends at a singular point).** Let J = [α, β] of width w lie in one knot interval,
and let Φ_c(v,·) = C + P₂ + P₁ on J with C convex, P₂ of class C², and P₁ concave and continuous on J and
differentiable at β. Then

    sup_J Φ_c(v,·) ≤ max( Φ_c(v,β),  Φ_c(v,α) + [ P₁(β) − P₁(α) − P₁′(β) w ] ) + (w²/8) sup_J (−P₂″)⁺.

In particular, if α = τ₋ and P₁(t) = −c₁ (t − τ₋) ln(t − τ₋), the bracket equals c₁w:

    sup_J Φ_c(v,·) ≤ max( Φ_c(v,β),  Φ_c(v,α) + c₁ w ) + (w²/8) sup_J (−P₂″)⁺,

and symmetrically max( Φ_c(v,α), Φ_c(v,β) + c₂ w ) on the subinterval that ends at τ₊. No logarithm or derivative
of W is evaluated at the singular end; P₂″ = D_c″ minus the second derivative of the rest of W_c, which is bounded
on J.

*Proof.* By the proof of Lemma 3, C + P₂ lies below its chord plus (w²/8) sup_J(−P₂″)⁺. A concave function lies
below its tangent at β. The sum of the chord and the tangent is affine, so its maximum over J is at an end. At β
it is (C+P₂)(β) + P₁(β) = Φ_c(v,β). At α it is (C+P₂)(α) + P₁(β) − P₁′(β)w = Φ_c(v,α) − P₁(α) + P₁(β) − P₁′(β)w.
For the entropy term, P₁(α) = 0, P₁(β) = −c₁ w ln w and P₁′(β) = −c₁(ln w + 1). ∎

(The version with both tangents: for every θ ∈ [0,1], sup_J(C + P) ≤ max( θΦ(α) + (1−θ)(Φ(β) − m_β w),
θ(Φ(α) + m_α w) + (1−θ)Φ(β) ) when P is concave and differentiable at both ends, with m_α, m_β the slope of the
chord of C plus P′(α), P′(β). It is not needed for the kinds of Section 1.)

**Lemma 12 (enclosing the logarithm with the four operations).** Let x > 0 be a machine number, x = μ·2ᵉ with
e an integer and μ ∈ [0.70, 1.42] (an exact operation: the binary exponent is split off, and μ is doubled if it is
below 0.7071). With z = (μ − 1)/(μ + 1), |z| ≤ ζ = 9/50,

    ln x = e ln 2 + 2 Σ_{j=0}^{n} z^{2j+1}/(2j+1) + R,       |R| ≤ 2 ζ^{2n+3} / ( (2n+3)(1 − ζ²) ),

and ln 2 = 2 Σ_{j≥0} 3^{−(2j+1)}/(2j+1), whose partial sums are lower bounds with remainder below
2·3^{−(2n+3)}/((2n+3)(1 − 1/9)).

*Proof.* ln μ = 2 artanh z = 2 Σ z^{2j+1}/(2j+1) for |z| < 1, and the tail is bounded by a geometric series. ∎

The implementation evaluates z and the sum in outward-rounded interval arithmetic and adds the rational bound for
R (n = 12 gives |R| < 10⁻²¹); ln 2 is enclosed once in exact rational arithmetic. The logarithm is increasing, so
an interval is mapped end by end. **Nothing new is trusted**: no library logarithm enters a certified number.

## 4. The competitor

u_h is the Q1 function of Proposition 7 (with equal nodal values on identified periodic nodes).

**Proposition 13.**

1. *Linear term.* ∫_c a u_h = |c| Σ_{v,v′} a_v U_{v′} Π_j m(v_j, v′_j) with m as in Proposition 7 (the integral of
   the product of two multilinear functions).
2. *Entropy terms and other potentials that are not polynomials.* The midpoint bound of Proposition 7 holds
   whenever W_c″ is bounded on [min_v U_v, max_v U_v]; W_c″ = P″ + c₁/(t − τ₋) + c₂/(τ₊ − t) is enclosed by
   interval arithmetic. If W_c is convex on that interval, also
   |c| W_c(ū_c) ≤ ∫_c W_c(u_h) ≤ |c| 2⁻ⁿ Σ_v W_c(U_v).
3. *Tables.* Write the piecewise linear W_c = W⁺ − W⁻ with W⁺, W⁻ convex piecewise linear (W⁻ collects the kinks
   with a negative jump of the slope). Then
   |c| [ W_c(ū_c) − e⁻ ] ≤ ∫_c W_c(u_h) ≤ |c| [ 2⁻ⁿ Σ_v W_c(U_v) + e⁻ ],  e⁻ = 2⁻ⁿ Σ_v W⁻(U_v) − W⁻(ū_c) ≥ 0.
4. *Exact mean.* Let φ_h be the Q1 function with nodal value 1 at every node that is not on Γ_D and 0 on Γ_D, and
   δ = ( m|Ω| − ∫u_h ) / ∫φ_h. Then u_h + δφ_h ∈ 𝒜 has mean exactly m. Its nodal values lie in
   [U + δ₋φ, U + δ₊φ] for any enclosure [δ₋, δ₊] of δ, and the certified upper bound is the upper end of the
   energy enclosure evaluated on these nodal intervals.

*Proof.* 1 is a direct computation. 2: if W_c is convex on the range of u_h over c, then W_c∘u_h is convex in each
coordinate (u_h is affine in each), hence below the multilinear interpolant of its vertex values; the lower bound
is Jensen's inequality, the mean of u_h over c being ū_c. 3: apply 2 to W⁺ (upper bound by vertex values) and to
W⁻ (lower bound by the centre value) and subtract; the other inequality likewise. 4: ∫(u_h + δφ_h) = m|Ω| by the
choice of δ, the boundary values on Γ_D are unchanged, and the energy of a Q1 function is a continuous function of
its nodal values, enclosed over the box of nodal values by interval arithmetic. ∫u_h and ∫φ_h are |c| times sums
of vertex means, enclosed in the same arithmetic; ∫φ_h > 0 unless every node lies on Γ_D. ∎

## 5. Linear potentials (no superlinear growth)

For W(x,t) = w₀ + A(x) t (A = w₁ + a, multilinear on each cell) the supremum in Theorem 8 is finite only if
div σ₀ = A everywhere, which stored machine numbers cannot satisfy exactly, and which the fields of part 1
(divergence constant on each cell) cannot satisfy at all when A varies inside a cell. Two routes exist: make the
equilibrium exact in rational arithmetic, or bound the residual and pay for it out of the gradient term. **The
second is taken**, for three reasons. Exact equilibrium needs a parametrisation of the divergence-free part of the
flux space (stream functions in two dimensions, edge potentials and harmonic fields in three, with cases for the
boundary conditions) and a particular field of higher degree for a source that varies in a cell. The residual
route costs only terms of the size of the rounding errors when the source is constant on each cell, and O(h²)
otherwise. And its constants are elementary for boxes.

**Lemma 14 (Poincaré inequality on a box).** Let c be a box with edge lengths hᵢ and w ∈ H¹(c) with ∫_c w = 0.
Then ‖w‖_{L²(c)} ≤ (h_max/π) ‖∇w‖_{L²(c)}.

*Proof.* The functions φ_μ(x) = Πᵢ cos(μᵢπ(xᵢ − aᵢ)/hᵢ), μ ∈ ℕ₀ⁿ, are a complete orthogonal system of L²(c),
−Δφ_μ = λ_μφ_μ with λ_μ = π² Σ μᵢ²/hᵢ², and ∂_νφ_μ = 0 on ∂c. Hence ∫∇w·∇φ_μ = λ_μ∫wφ_μ for w ∈ H¹(c), and
∫∇φ_μ·∇φ_μ′ = λ_μ∫φ_μφ_μ′ = 0 for μ ≠ μ′. Bessel's inequality in L²(c; ℝⁿ) for the orthogonal system {∇φ_μ}_{μ≠0}
gives ‖∇w‖² ≥ Σ_{μ≠0} (∫∇w·∇φ_μ)²/‖∇φ_μ‖² = Σ_{μ≠0} λ_μ (∫wφ_μ)²/‖φ_μ‖² ≥ (π/h_max)² ‖w‖², by Parseval's identity
for w, whose coefficient at μ = 0 vanishes. ∎

**Lemma 15 (along lines).** Let F = {xᵢ = aᵢ} or {xᵢ = bᵢ} be a face of Ω, π_F the orthogonal projection onto F,
u ∈ H¹(Ω) and γ its trace on F. Then ‖u − γ∘π_F‖_{L²(Ω)} ≤ (Lᵢ/√2) ‖∂ᵢu‖_{L²(Ω)}.

*Proof.* For u ∈ C¹(Ω̄), |u(x) − u(π_Fx)|² ≤ |xᵢ − aᵢ| ∫_{aᵢ}^{bᵢ} |∂ᵢu|² dxᵢ by the fundamental theorem of
calculus and the Cauchy–Schwarz inequality; integrate over Ω. Both sides are continuous on H¹(Ω) (the trace is
continuous into L²(F)) and C¹(Ω̄) is dense. ∎

**Proposition 16 (the bound for linear potentials).** Let σ₀ be a lowest-order Raviart–Thomas field with σ₀·ν = 0
on the natural faces and equal values on identified periodic faces, d_c = div σ₀ on c, Ā_c the mean of A over c
(the mean of its vertex values), ρ_c = d_c − Ā_c, and

    ‖σ₀‖² = ∫_Ω |σ₀|²,      ‖ρ‖² = Σ_c |c| ρ_c²,      R₁² = Σ_c ∫_c (A − Ā_c)²,      Fl = ∫_{Γ_D} g σ₀·ν.

(a) If the face F ⊂ Γ_D is orthogonal to eᵢ, then for every u ∈ 𝒜

    E[u] ≥ w₀|Ω| + Fl − ( ‖σ₀‖ + (h_max/π) R₁ + (Lᵢ/√2) ‖ρ‖ )² / (4k) − ‖ρ‖ √|Ω| max_F |g|.

(b) If the mean is prescribed, then for every u ∈ 𝒜 and every λ ∈ ℝ, with A replaced by A − λ throughout,

    E[u] ≥ w₀|Ω| + λ m|Ω| + Fl − ( ‖σ₀‖ + (h_max/π) R₁ + (L_max/π) ‖ρ‖ )² / (4k) − m Σ_c |c| ρ_c.

*Proof.* Let X = ‖∇u‖_{L²(Ω)}. As in step 3 of Theorem 1 with S = tσ₀, ∫(div σ₀)u + ∫σ₀·∇u = Fl. Since A − Ā_c
has mean zero on c,

    ∫_Ω A u = Σ_c ∫_c Ā_c u + Σ_c ∫_c (A − Ā_c)(u − ū_c) = Fl − ∫σ₀·∇u − ∫ρu + Σ_c ∫_c (A − Ā_c)(u − ū_c).

By Lemma 14 and the Cauchy–Schwarz inequality the last sum is at least −(h_max/π) R₁ X, and ∫σ₀·∇u ≤ ‖σ₀‖X.
(a) By Lemma 15, ‖u‖ ≤ ‖g∘π_F‖ + (Lᵢ/√2)X and ‖g∘π_F‖_{L²(Ω)} = √Lᵢ ‖g‖_{L²(F)} ≤ √|Ω| max_F|g|, so
|∫ρu| ≤ ‖ρ‖√|Ω| max_F|g| + (Lᵢ/√2)‖ρ‖X. (b) u − m has mean zero on Ω, so by Lemma 14 on Ω,
−∫ρu = −m∫ρ − ∫ρ(u − m) ≥ −mΣ|c|ρ_c − (L_max/π)‖ρ‖X; the multiplier enters as in Corollary 9. In both cases
E[u] − w₀|Ω| = kX² + ∫Au ≥ kX² − BX + (terms without X) with B the bracket of the statement, and
kX² − BX ≥ −B²/(4k). ∎

*Evaluation.* ‖σ₀‖² is the exact formula of Proposition 6; ∫_c(A − Ā_c)² = |c| ( Σ_{v,v′} A_vA_{v′}Π_j m(v_j,v′_j)
− Ā_c² ), and it is zero when A is constant on c; ρ_c, Fl and the sums are enclosed by interval arithmetic. An
upper bound s for a square root is a machine number s whose square, rounded down, is verified to be at least the
upper bound of the radicand. π enters through an exact rational lower bound, obtained from Machin's formula
π = 16 arctan(1/5) − 4 arctan(1/239) with alternating series (partial sums bracket the limit). max_F|g| is the
largest modulus of the nodal values on F, because g is multilinear between them.

*Tightness.* When A is constant on each cell and the flux of the minimiser lies in the flux space, the bound
equals the minimum up to ‖ρ‖, which is the residual of a linear solve. A problem with neither a Dirichlet face nor
a prescribed mean is refused: its energy is not bounded below unless ∫A = 0, and then it does not have a unique
scale.

## 6. Choosing the fields (heuristic)

Nothing here is trusted.

* *Convex mode with a linear term or a prescribed mean.* The finite-volume problem of part 1 with the cell means
  Ā_c added to W′, and with the constraint Σ_c p_c = m·(number of cells) when the mean is prescribed; λ is the
  multiplier of that constraint.
* *Linear potentials.* The same finite-volume system is linear and is solved directly, so that ρ is of the size
  of the rounding errors.
* *Calibration mode.* The smoothed functional of part 1 with the per-vertex term −a_v t, the shifted tails, the
  end rule of Lemma 11 in place of the curvature allowance on the two end subintervals, and λ as one more unknown
  (the gradient with respect to λ is m|Ω| minus the smoothed mean of t).
* *Candidate with a prescribed mean.* The Q1 energy is minimised over nodal values with the weighted mean fixed
  (a uniform shift of the free nodes restores it after every step); for potentials on an interval the potential is
  continued outside by its second-order Taylor polynomial in the solver only, and the result is clipped into the
  interval.

## 7. Known answers that the tests are to use

All are elementary; the proofs are one line each.

1. *Source, constant* (E = ∫k|∇u|² − f u on (0,1) × cross-section, u = 0 on x₁ = 0, 1, natural elsewhere, f = 4k):
   u* = x₁(1 − x₁), min E = −k/3 per unit cross-section. The flux 2k(1 − 2x₁) lies in the flux space, so the lower
   bound must equal the minimum to rounding.
2. *Source, linear* (f = c x₁): u* = c x₁(1 − x₁²)/(12k), min E = −c²/(180k) per unit cross-section.
3. *Convex potential with a prescribed mean*, natural or periodic faces: min E = |Ω| W(m) (Jensen).
4. *Non-convex potential with a prescribed mean on a small box.* If W″ ≥ −2κ on 𝕀, all faces are natural or
   periodic and k ≥ κ (L_max/π)², then u ≡ m is a minimiser and min E = |Ω| W(m): W(t) + κ(t − m)² is convex, and
   κ∫(u − m)² ≤ k∫|∇u|² by Lemma 14. For the regular-solution energy, W″ ≥ 4 − 2χ, so κ = χ − 2.
5. *A table that is the same in every cell*, natural faces: min E = |Ω| min W_c.
6. *One-dimensional minimisers*: a strip with data that do not depend on x₂ and periodic or natural faces in x₂
   has the energy of the one-dimensional problem times the width, to be compared with `ScalarField1D` or with a
   one-dimensional reference value (floating point, labelled as such).

## 8. What is not covered

* **Not claimed:** anything about a problem class before its known-answer and falsification tests have passed.
* A matrix K that is not a multiple of the identity (Theorem 8 and Proposition 10 cover it; the competitor energy
  and the heuristics were not written for it).
* Potentials on an interval and tables in convex mode (they are evaluated in calibration mode, which is valid for
  every potential); linear potentials in calibration mode.
* A position-dependent term that is not affine in u, apart from the kinds of Section 1.
* Total-variation and other non-quadratic gradient terms.

## References

1. L. E. Payne, H. F. Weinberger, *An optimal Poincaré inequality for convex domains*, Arch. Rational Mech. Anal.
   5 (1960), 286–292; M. Bebendorf, *A note on the Poincaré inequality for convex domains*, Z. Anal. Anwend. 22
   (2003), 751–756. (The constant diam/π for convex domains; for a box Lemma 14 is elementary. Both cited from
   memory; the bibliographic data are still to be verified.)
2. J. W. Cahn, J. E. Hilliard, *Free energy of a nonuniform system. I. Interfacial free energy*, J. Chem. Phys. 28
   (1958), 258–267. (Gradient energy with the regular-solution free energy. Cited from memory; to be verified.)
3. The references of part 1, in particular [1]–[3] there for the inequality and [4] for the complementary energy.

The literature search for this part has not been done yet; it is listed as open in the log.
