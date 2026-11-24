# The level-set lower bound for Skyrme fields of degree two: audit and extension

This note audits the level-set lower bound behind the model `SkyrmeDegreeOne` for its use of the degree, states
the bound for maps of degree ±B, proves it, and records what the certified programme gives for B = 2.

**Status.** The analytic reduction extends to every degree B ≠ 0 with one change: the flux through a level
surface, c(τ) = 4π sin²τ, becomes |B| c(τ). There is no obstruction. Theorems 1 and 2 below are proved in this
note; the proofs follow those of the degree-one note step by step and have been checked by their author only. The
source remarks, without proof, that the bound for other degrees is obtained "by replacing c by |B|c" and quotes a
floating-point value 2.209 for degree two; no novelty is claimed for the extension. Numbers called *certified*
come from the outward-rounded programme `Kernel/SkyrmeLevelSetDegree.wl`; numbers called *floating point* are
ordinary numerics. The state of the computations is in the log at the end.

**Sources read for the audit.** (a) `Kernel/SkyrmeLevelSet.wl`, `Kernel/Minorants.wl` and the model text in
`Kernel/Models.wl`. (b) The level-set note, "A level-set lower bound for the Skyrme energy" (working note dated
12 November 2026; reference [2] of the README), Sections 1 to 6, in a public copy, and the listings of its
programmes that accompany it. (c) The note on degree two (reference [1] of the README) in its public version of
12 November 2026. The record that the README cites for [1] could not be retrieved on 23 November 2026 (the DOI
answers "record deleted"), so the version of [1] that the README refers to was not read. Every step used below is
written out here, so this note can be checked without those sources.

Units as in the paclet: E(U) = ∫_{ℝ³} |dU|² + |Λ²dU|² + 2μ²(1 − cos f) dx, f = dist_{S³}(U, 1), reported as
E/12π². Faddeev's bound is E ≥ 12π²|B|.

## 1. Statements

Notation. c(τ) = 4π sin²τ; g(v) = (36π)^{1/3} v^{2/3}; V(τ) = 2μ²(1 + cos τ);
h_κ(D) = 2√((1 + 2D)(D² + 2D + κ)) for κ, D ≥ 0. A *minorant* is a triple m = (α, β, κ) ∈ [0, ∞)³ with
α + βD ≤ h_κ(D) for all D ≥ 0.

The class. For an integer B ≠ 0, 𝒞_B is the set of measurable U : ℝ³ → S³ ⊂ ℝ⁴ with weak derivative,
E(U) < ∞, U − 1 ∈ L⁶(ℝ³), and deg U := (1/2π²) ∫ det dU dx = B. (If dU ∈ L² there is a constant U_∞ with
U − U_∞ ∈ L⁶; for μ > 0 finite energy forces U_∞ = 1, for μ = 0 it is a normalisation by a rotation of S³.)

The one-dimensional problem. Let b > 0, τ ∈ (0, π), v ≥ 0. A pair (A, M) of measurable functions
[τ, π) → [0, ∞), with M integrable on [τ, s] for every s < π, is *admissible from v* if A(s) ≥ g(v̄(s)) for
almost every s, where v̄(s) = v + ∫_τ^s M. Put

    Q_b(s) = sup_m [ αA(s) + β b c(s) − κM(s) ] + V(s)M(s),        W_b(τ, v) = inf ∫_τ^π Q_b(s) ds,

the supremum over all minorants, the infimum over the pairs admissible from v. Q_b ≥ 0 because (0, 0, 0) is a
minorant. W_b(τ, ·) is nondecreasing: g is increasing, so a pair admissible from v₂ is admissible from every
v₁ < v₂, with the same Q_b.

**Theorem 1 (level-set inequality in degree B).** Let B ≠ 0 and U ∈ 𝒞_B. Then for every τ₁ ∈ (0, π)

    E(U)  ≥  inf_{v ≥ 0} W_{|B|}(τ₁, v).

**Theorem 2 (the certified programme).** For every positive integer B, every pion mass μ ∈ ℚ, μ ≥ 0, every grid
(levels 0 < τ₀ < … < τ_N < π and increasing volume nodes at each level) and every certified table of minorants,
the number returned by `levelSetLowerBoundDegree[B, μ, grid, table]`, multiplied by 12π², is at most
inf_{v ≥ 0} W_B(τ₀, v). Hence E(U)/12π² is at least that number for every U ∈ 𝒞_B ∪ 𝒞_{−B}.

For B = 1 these are the statements behind `SkyrmeDegreeOne`.

**Relations between degrees (proved in Section 4).** W_B ≥ 6B ∫_τ^π c, so the bound is never below Faddeev's
up to the part of the integral that the grid leaves out. W_B(τ, Bv) ≤ B·W_1(τ, v), so the bound for degree B is at
most B times the value of the degree-one problem. It is strictly smaller in the cases computed: the bound is not
linear in B. For B ≥ 2 the extremal of the relaxed problem starts with *more* area than the isoperimetric
inequality requires (the level surface of B separate cores), and the isoperimetric constraint becomes active only
later (near τ = 1.05 for B = 2, μ = 0).

## 2. The audit: where the degree-one proof uses the degree

| Step of the degree-one argument | Uses the degree? | In degree B |
|---|---|---|
| Pointwise inequality on a level surface (Lemma A) | no | unchanged |
| Coarea formula on the level sets of θ = dist(U, −1) | no | unchanged |
| Flux: ∫_{Σ_τ} D dℋ² ≥ c(τ), from the degree formula | **yes, the only place** | ∫_{Σ_τ} D dℋ² ≥ \|B\| c(τ) |
| Volume: \|{θ < τ}\| ≥ ∫₀^τ M | no | unchanged |
| Isoperimetric inequality A(τ) ≥ g(v(τ)) | no | unchanged; valid, though level surfaces of a map of higher degree need not be near a sphere, which costs sharpness and not validity |
| Minorants (M1), (M2) and the table | no (properties of h_κ) | unchanged; the same table is used |
| (M3): level energy ≥ αA + βc − κM | through the flux | αA + β\|B\|c − κM |
| Monotonicity of W in the volume | no | unchanged |
| Floor Bog = 6∫c (minorant (0, 6, 0)) | through the flux | 6\|B\|∫c |
| Floors Ξ, Ξ_μ (minorants (2√κ, 0, κ)) | no: β = 0, the flux does not enter | unchanged; they do not improve with B |
| Second-order step bound (Lemma D) | through the flux integrals C_a, C_b | C_a, C_b multiplied by \|B\| |
| Backward induction: chords, windows, far targets, monotone hull | through C₀, C_a, C_b only | unchanged otherwise |
| Regularity of level sets (Sobolev maps) | no | unchanged; the degree formula (Lemma B) holds for every degree |
| Normalisation at infinity, U − 1 ∈ L⁶ | no | unchanged |

Three remarks from the audit.

1. *The closed-form floors.* The written proof of the floors in the source integrates a differential inequality
   up to τ = π using Y(π, ·) = 0. For a map whose set {θ < π} has infinite volume (the hedgehog, for example) the
   volume v̄(s) tends to infinity as s → π and the end value needs an argument. The statement is correct; Lemma C
   below proves it for pairs that are admissible on [τ, π) with locally finite volume, which is what a map
   provides. This concerns degree one as much as degree B.
2. *Sub-additivity.* The inequality I₂ ≤ 2I₁ between the infima is not used anywhere below. According to
   Esteban's erratum [E] it is not known whether I_k ≤ I_ℓ + I_{k−ℓ} holds, because a finite-energy map cannot yet
   be cut off to a map that is constant near infinity. The upper bound used for degree two is the energy of an
   explicit map: two disjoint copies of a trial map that is the vacuum outside a ball.
3. *Jensen's inequality.* Theorem 1 of the source passes through the convex envelope ψ of the level energy and
   Jensen's inequality. The certified statement does not need it: (M3) uses one affine minorant at a time. The
   proof below takes that route.

## 3. Proof of Theorem 1

Throughout, U ∈ 𝒞_B, θ = dist_{S³}(U, −1) ∈ [0, π] in its precise representative, and b = |B|.

**Facts quoted from the literature.** (F1) θ ∈ W^{1,2}_loc; almost everywhere on {U ≠ ±1}, ∇θ = dUᵀ ν(U) with ν
the unit normal of the geodesic spheres about −1, and ∇θ = 0 almost everywhere on {U = ±1}. (F2) The coarea
formula for Sobolev functions [MSZ]: ∫ f |∇θ| dx = ∫₀^π ∫_{θ⁻¹(τ)} f dℋ² dτ for measurable f ≥ 0. (F3) The
coarea formula for functions of bounded variation, De Giorgi's structure theorem and the isoperimetric
inequality P(Ω) ≥ g(|Ω|) for sets of finite perimeter and finite volume [EG]. These are used exactly as in the
source; they were not re-derived here.

Write Σ*_τ = θ⁻¹(τ) ∩ {∇θ ≠ 0}. By (F2) with f = 1_{∇θ = 0}, ℋ²(θ⁻¹(τ) ∖ Σ*_τ) = 0 for almost every τ, and a
property that holds almost everywhere in {∇θ ≠ 0} holds ℋ²-almost everywhere on Σ*_τ for almost every τ.

**Lemma A (pointwise).** Let x be a point of differentiability with U(x) ≠ ±1 and a = ∇θ(x) ≠ 0. Let T = a^⊥,
let H be the tangent plane at U(x) of the geodesic sphere {dist(·, −1) = θ(x)}, and F = dU(x). Then F maps T into
H, and with D = |det(F|_T : T → H)|

    |F|² + |Λ²F|²  ≥  |a|²(1 + 2D) + 2D + D².

*Proof.* dθ = ⟨ν, F·⟩, so ⟨ν, Ft⟩ = a·t = 0 for t ∈ T: F(T) ⊂ H. Let n = a/|a|, τ₁, τ₂ an orthonormal basis of T
and G = F|_T. Then Fn = |a|ν + h with h ∈ H, and |F|² = |G|² + |a|² + |h|². In the basis τ₁∧τ₂, n∧τ₁, n∧τ₂ of
Λ²ℝ³, |Λ²F|² = |Gτ₁∧Gτ₂|² + Σᵢ |(|a|ν + h)∧Gτᵢ|² = (det G)² + |a|²|G|² + Σᵢ|h∧Gτᵢ|², because ν∧Gτᵢ is orthogonal
to Λ²H ∋ h∧Gτᵢ. Drop the terms with h and use |G|² ≥ 2|det G| = 2D. ∎

**Step 1 (coarea).** With ν̂ := 1/|∇θ| on Σ*_τ, φ(D, ν̂) := (1 + 2D)/ν̂ + (2D + D²)ν̂ and
M(τ) := ∫_{Σ*_τ} ν̂ dℋ²,

    E(U)  ≥  ∫₀^π [ ∫_{Σ*_τ} φ(D, ν̂) dℋ² + V(τ) M(τ) ] dτ.                                  (1)

Indeed E(U) ≥ ∫_{∇θ ≠ 0} (W(dU) + V(θ)) dx, by (F2) this is ∫₀^π ∫_{Σ*_τ} (W(dU) + V(τ))/|∇θ| dℋ² dτ, and by
Lemma A, W(dU)/|∇θ| ≥ |a|(1 + 2D) + (2D + D²)/|a| = φ(D, ν̂). (The mass term is 2μ²(1 − cos f) = V(θ).)

**Lemma B (degree formula).** Let U be weakly differentiable with dU, Λ²dU ∈ L² and U − 1 ∈ L⁶. Then for every
continuous ψ on S³, ∫ ψ(U) det dU dx = deg U · ∫_{S³} ψ ω, with ω the volume form of S³.

*Proof.* (As in the source; the degree is arbitrary.) By density take ψ smooth. Since H³(S³) = ℝ,
ψω − ψ̄ω = dβ with ψ̄ the mean of ψ and β a smooth 2-form, which after subtracting a closed form satisfies
|β(u)| ≤ C|u − 1|. For U ∈ W^{1,2}_loc ∩ L^∞ the forms dUⁱ∧dUʲ are closed in distributions and lie in L², so
d(U*β) = U*(dβ) in distributions. Testing with a cut-off χ_R (χ_R = 1 on B_R, |∇χ_R| ≤ 2/R, support in B_{2R})
and Hölder's inequality with exponents 6, 2, 3 on the annulus B_{2R}∖B_R, of volume CR³,

    |∫ χ_R U*(dβ)| = |∫ dχ_R ∧ U*β| ≤ (C/R) ∫_{B_{2R}∖B_R} |U − 1| |Λ²dU| ≤ C ‖U − 1‖_{L⁶(ℝ³∖B_R)} ‖Λ²dU‖_{L²(ℝ³∖B_R)} → 0.

Hence ∫ U*(ψω) = ψ̄ ∫ U*ω = deg U ∫_{S³} ψω. ∎

**Step 2 (three facts, for almost every τ ∈ (0, π)).**

(i) *Flux.* ∫_{Σ*_τ} D dℋ² ≥ b c(τ). At almost every point with ∇θ ≠ 0, in the frames of Lemma A the matrix of F
is block triangular (Fτᵢ ∈ H, Fn = h + |a|ν), so det dU = |∇θ| det(F|_T) with the signed determinant; where
∇θ = 0 and U ≠ ±1 the image of dU lies in H and det dU = 0. For χ ∈ C_c(0, π), Lemma B with
ψ = χ(dist(·, −1)) and (F2) give

    ∫₀^π χ(τ) ∫_{Σ*_τ} det(F|_T) dℋ² dτ = ∫ χ(θ) det dU dx = B ∫₀^π χ(τ) c(τ) dτ,

since the geodesic sphere of radius τ in S³ has area c(τ). So ∫_{Σ*_τ} det(F|_T) dℋ² = B c(τ) for almost every τ,
and ∫ D ≥ |∫ det(F|_T)| = b c(τ). **This is the only place where the degree enters.**

(ii) *Volume.* v(τ) := |{θ < τ}| is finite for τ < π because U − 1 ∈ L⁶, and by (F2) with f = 1/|∇θ| on
{∇θ ≠ 0, θ < τ}: v(τ) ≥ ∫₀^τ M.

(iii) *Isoperimetry.* A(τ) := ℋ²(Σ*_τ) ≥ g(v(τ)): by (F3), for almost every τ the set {θ < τ} has finite
perimeter, its reduced boundary lies in θ⁻¹(τ) up to an ℋ²-null set, and P({θ < τ}) ≤ ℋ²(θ⁻¹(τ)) = A(τ).

**Step 3 (one minorant at a time).** Let m = (α, β, κ) be a minorant. By the inequality between the arithmetic
and geometric means, φ(D, ν̂) + κν̂ = (1 + 2D)/ν̂ + (D² + 2D + κ)ν̂ ≥ h_κ(D) ≥ α + βD. Integrating over Σ*_τ and
using (i) and β ≥ 0,

    ∫_{Σ*_τ} φ(D, ν̂) dℋ²  ≥  αA(τ) + β b c(τ) − κM(τ)                                         (2)

for every τ at which (i)–(iii) hold and A(τ), M(τ) are finite, simultaneously for all minorants. Both are finite
for almost every τ ∈ (0, π): M by (ii), and A because by (F2), for τ′ < τ″ < π,
∫_{τ′}^{τ″} A dτ = ∫_{τ′<θ<τ″} |∇θ| dx ≤ |{θ < τ″}|^{1/2} ‖dU‖_{L²}. (If ∫ D = ∞ the left side of (2) is infinite,
since φ ≥ 4D.)

**Step 4.** By (1) and (2), E(U) ≥ ∫₀^π Q_b(s) ds with the level data (A, M) of U. For τ₁ ∈ (0, π) put
v₁ = ∫₀^{τ₁} M. On [τ₁, π), v̄(s) = v₁ + ∫_{τ₁}^s M = ∫₀^s M ≤ v(s) < ∞ by (ii), and A(s) ≥ g(v(s)) ≥ g(v̄(s)) by
(iii) and the monotonicity of g. So (A, M) is admissible from v₁, and since Q_b ≥ 0,

    E(U) ≥ ∫_{τ₁}^π Q_b ≥ W_b(τ₁, v₁) ≥ inf_{v ≥ 0} W_b(τ₁, v).  ∎

## 4. The one-dimensional problem: what the programme uses

**(M1)** Both sides being nonnegative, m is a minorant if and only if
P_m(D) = 8D³ + (20 − β²)D² + (8 + 8κ − 2αβ)D + 4κ − α² ≥ 0 on [0, ∞). **(M2)** κ ↦ h_κ(D) is concave, so convex
combinations of minorants are minorants. Neither contains the flux; the certified table of
`Kernel/Minorants.wl` is used unchanged.

**Lemma C (floors).** With r = (3v/4π)^{1/3}, each of

    Bog_b(τ) = 6b ∫_τ^π c,        Ξ(τ, v) = 4π(π − τ)² r,        Ξ_μ(τ, v) = 16π cos²(τ/2)(r + μr²)   (μ > 0)

is ≤ W_b(τ, v). Ξ and Ξ_μ are concave and nondecreasing in v.

*Proof.* (0, 6, 0) is a minorant, since h₀(D)² − 36D² = 8D(D − 1)², so Q_b ≥ 6bc: this is Bog_b. For the other
two let (A, M) be admissible from v with finite cost. For every κ ≥ 0, (2√κ, 0, κ) is a minorant (h_κ is
increasing in D), so Q_b ≥ sup_κ [2√κ A − κM] + VM; in particular M > 0 almost everywhere on {A > 0}, and
Q_b ≥ A²/M + VM ≥ g(v̄)²/M + VM. In terms of the volume radius r(s) of v̄(s), g(v̄) = 4πr² and M = 4πr²r′.

Ξ: if r(τ) = 0 there is nothing to prove. Otherwise for τ < s₁ < π, by the Cauchy–Schwarz inequality,
(s₁ − τ)² ≤ ∫_τ^{s₁} (r²/r′) ds · ∫_τ^{s₁} (r′/r²) ds ≤ ∫_τ^{s₁} (r²/r′) ds / r(τ), hence
∫_τ^π Q_b ≥ 4π ∫_τ^{s₁} r²/r′ ≥ 4π r(τ)(s₁ − τ)², and s₁ → π gives Ξ(τ, v). No assumption on v̄ near π is needed.

Ξ_μ: put L = 2cos(τ/2), so V = μ²L², Ξ_μ = 4πL²ϕ(r) with ϕ(r) = r + μr², which satisfies r√(μ²r² + ϕ′) = ϕ.
Taking κ = V + ∂_vΞ_μ = μ²L² + L²ϕ′/r² in the bound above, evaluated along the path,

    Q_b ≥ 2√κ g(v̄) − ∂_vΞ_μ M = 8πLϕ − ∂_vΞ_μ v̄′ ≥ 8πL sin(s/2)ϕ − ∂_vΞ_μ v̄′ = −(d/ds) Ξ_μ(s, v̄(s)),

so Ξ_μ(τ, v) ≤ Ξ_μ(s₁, v̄(s₁)) + ∫_τ^{s₁} Q_b for every s₁ < π (the function s ↦ Ξ_μ(s, v̄(s)) is absolutely
continuous on [τ, s₁]). It remains to see that lim inf_{s₁→π} Ξ_μ(s₁, v̄(s₁)) = 0. Suppose Ξ_μ(s, v̄(s)) ≥ ε > 0
for all s ≥ s₂. Since L² ≤ (π − s)², this forces r(s) ≥ c₁/(π − s) and v̄(s) ≥ c₂/(π − s)³ near π. On the other
hand Q_b ≥ VM and V is decreasing, so V(s′)(v̄(s′) − v̄(s)) ≤ ∫_s^{s′} VM ≤ ∫_τ^π Q_b < ∞, while
V(s′) ≥ (4μ²/π²)(π − s′)², so V(s′)v̄(s′) ≥ c₃/(π − s′) → ∞: a contradiction. ∎

**Lemma D (second-order bound for one step).** Let 0 < τ_a < τ_b < π, Δ = τ_b − τ_a, t = (s − τ_a)/Δ, let m_a, m_b
be minorants, (A, M) admissible from v at τ_a, w = v̄(τ_b), w̄ ≥ w and γ = g′(w̄). If

    κ_a − κ_b ≤ Δ (γ min(α_a, α_b) + 2μ² min(sin τ_a, sin τ_b)),                              (3)

then

    ∫_{τ_a}^{τ_b} Q_b  ≥  g(v) Δ (α_a + α_b)/2 + b (β_a C_a + β_b C_b) − (κ_b − V(τ_b))(w − v),

with C_a = ∫_{τ_a}^{τ_b} (1 − t) c and C_b = ∫_{τ_a}^{τ_b} t c.

*Proof.* By (M2), m(s) = (1 − t)m_a + t m_b is a minorant for every s, so Q_b ≥ α(s)A + β(s) b c − K(s)M with
K = κ(s) − V(s). Since g is concave and v ≤ v̄(s) ≤ w ≤ w̄, A ≥ g(v̄(s)) ≥ g(v) + γ(v̄(s) − v). Integrating
∫ K v̄′ by parts,

    ∫ Q_b ≥ g(v) ∫α + b ∫βc − K(τ_b)(w − v) + ∫ [γα(s) + K′(s)] (v̄(s) − v) ds.

K′ = (κ_b − κ_a)/Δ + 2μ² sin s and sin is concave on [0, π], so (3) makes the bracket nonnegative, and
v̄ − v ≥ 0. ∎

The flux enters Lemma D through b C_a and b C_b, and the far-target bounds through b C₀ with C₀ = ∫ c over the
step: ∫ Q_b ≥ ρ b C₀ for a minorant (0, ρ, κ) with κ ≤ V(τ_b), and ∫ Q_b ≥ 6 b C₀ + V(τ_b)(w − v).

**The backward induction (proof of Theorem 2).** Levels τ₀ < … < τ_N; at level i nodes
0 = v_{i,0} < … < v_{i,J_i}. The programme constructs functions Λ_i that are affine on each cell
[v_{i,j}, v_{i,j+1}] (jumps at the nodes are allowed) and constants T_i with

    Λ_i ≤ W_b(τ_i, ·) on every cell,        T_i ≤ W_b(τ_i, v) for all v ≥ v_{i,J_i}.            (H_i)

*Start.* On each cell Λ_N is the chord of Ξ(τ_N, ·) or of Ξ_μ(τ_N, ·), which lies below these concave functions,
and T_N is the larger of their values at the last node (they are nondecreasing in v). (H_N) holds by Lemma C.

*Step.* Assume (H_{i+1}). Fix a cell I = [v_{i,j}, v_{i,j+1}], a point v ∈ I and a pair admissible from v with
finite cost; let w = v̄(τ_{i+1}) ≥ v. The restriction of the pair to [τ_{i+1}, π) is admissible from w, so its cost
is at least ∫_{τ_i}^{τ_{i+1}} Q_b + W_b(τ_{i+1}, w). For a target cell k of level i+1 inside a finite window, the
programme has chosen one pair of minorants satisfying (3) with w̄ = v_{i+1,k+1}, and one affine function
L_k ≤ W_b(τ_{i+1}, ·) on cell k: the affine piece of Λ_{i+1}, or the constant sup of Λ_{i+1} over the cells
below, which is a lower bound because W_b is nondecreasing in v. If w lies in cell k, then
w ∈ R_k = [max(v_{i+1,k}, v_{i,j}), v_{i+1,k+1}], and by Lemma D the cost is at least an expression that is affine
in w; its minimum over R_k is attained at an end of R_k and is a concave function b_k(v) of v (g is concave,
α ≥ 0, the rest is affine in v). If w lies above the window or beyond the last node, one of the two far-target
bounds, added to the infimum of the corresponding lower bounds for W_b(τ_{i+1}, ·) over those targets, gives a
function b_∞(v) that is affine in v. Cells below the one containing v_{i,j} cannot be reached, since w ≥ v. Hence
W_b(τ_i, v) ≥ min(min_k b_k(v), b_∞(v)) on I. The right-hand side is concave, so it lies above its chord; Λ_i on
I is that chord, or the chord of Ξ(τ_i, ·) or Ξ_μ(τ_i, ·), or the constant 6b ∫_{τ_i}^{τ_N} c, each of which is a
lower bound by Lemma C. T_i is the largest of: the values of Ξ and Ξ_μ at the last node, that constant, a
far-target bound over all targets above the last node, and the largest value of Λ_i (monotonicity). This is
(H_i).

*End.* By (H₀) and monotonicity, inf_v W_b(τ₀, v) is at least the minimum over the cells of the running maximum
of the cell minima of Λ₀, and at least T₀; the programme returns the smaller of the two, divided by an upper bound
of 12π², rounded down. Every quantity is evaluated with the arithmetic of `SkyrmeLevelSet.wl`; the only new
rigorous quantities are the lower bounds of b C₀, b C_a, b C_b, obtained by multiplying exact rational lower
bounds of C₀, C_a, C_b by the integer b before rounding down. The choices of grids, windows, tangent columns and
multipliers affect only how large the returned number is. ∎

**Proposition E (relations between degrees).** (a) W_b(τ, v) ≥ 6b ∫_τ^π c. (b) For every integer b ≥ 1,
W_b(τ, bv) ≤ b W_1(τ, v).

*Proof.* (a) is the first floor of Lemma C. (b) If (A, M) is admissible from v for the problem with b = 1, then
(bA, bM) is admissible from bv: its volume is b v̄ and bA ≥ b g(v̄) ≥ b^{2/3} g(v̄) = g(b v̄). Its integrand is
sup_m [αbA + βbc − κbM] + VbM = b Q_1. ∎

By (b), the certified number for degree two can never exceed twice the value of the degree-one relaxed problem,
which in turn is at most twice the energy of any degree-one map. This is one of the checks of Section 7 and does
not rely on sub-additivity of the infima.

## 5. Implementation

* `Kernel/SkyrmeLevelSetDegree.wl`: `levelSetLowerBoundDegree[B, μ, grid, table]`. The body is derived from
  `levelSetLowerBound` by textual substitution: the factor B in the three flux integrals of a step
  (`C0`, `Ca`, `Cb`) and in the two floating-point numbers `ca`, `cb` that steer heuristic choices. The Bogomolny
  floor is accumulated from `C0` and so carries the factor automatically. Nothing else differs. For B = 1 the
  two functions perform the same operations.
* `Reference/SkyrmeDegree/relaxed.py`: floating-point extremals of the relaxed problem in degree B (not
  rigorous). `make_grids.py`: the grids, laid out as described in the level-set note, around the volume path of
  the extremal. For B = 1 and μ = 0 the generator reproduces the stored Quick grid of `Data/Skyrme` (same
  number of levels and nodes, levels to 10⁻⁹, reference volumes to 5×10⁻⁵).
* `Data/SkyrmeDegree`: grids for B = 2, and the stored claims.
* The model `SkyrmeDegreeTwo`. Lower bound: the programme with B = 2. Upper bound: the hedgehog trial map of
  `SkyrmeDegreeOne` is the vacuum outside a ball, so two translated copies with disjoint supports form one
  continuous map of degree two with twice the energy; the certified upper bound is twice the certified energy of
  the trial. Stored and reported claims are the computed bounds weakened by one part in 10⁹.

## 6. Numbers

Floating point (extremals of the relaxed problem, `relaxed.py`): R₂ = 2.20844 (μ = 0) and 2.61111 (μ = 1), in
units of 12π²; for comparison R₁ = 1.19210 and 1.41439, which agree with the level-set note. The source quotes
2.209 for degree two at μ = 0.

Certified (`levelSetLowerBoundDegree`, Wolfram Engine 15.0, one core). Every map of degree ±2 has E/12π² at least

| μ | Grid | Certified lower bound | Float value R₂ | Levels | Transitions | Fallbacks | Time |
|---|---|---|---|---|---|---|---|
| 0 | Quick | 2.207981 | 2.208440 | 1011 | 4,643,316 | 1,401,645 | 51 s |
| 1 | Quick | 2.572283 | 2.611114 (see below) | 887 | 5,081,334 | 66,696 | 77 s |

Both numbers are certified lower bounds and both passed the path test. The Standard grids exist but have not been
run. For μ = 0 the certified number is 0.02 % below the floating-point value of the relaxed problem, which is the
loss of the Quick grid (for degree one the same resolution loses 0.03 %); it is 10.4 % above Faddeev's bound 2
and 7.4 % below twice the degree-one bound (2 × 1.191748).

**Open point for μ = 1.** The certified number is 1.49 % below the cost 2.611114 of the extremal found by
`relaxed.py`. It is not yet known whether the Quick grid loses this much at μ = 1, or whether the relaxed problem,
which is not convex, has a cheaper admissible path than that extremal (the shooting finds one extremal, not
necessarily the minimiser). The certified bound is valid in either case; what is unsettled is the value of the
relaxed problem, and with it how much a finer grid can gain. Until this is settled, 2.611 must not be quoted as
the value of the relaxed problem for μ = 1. See the log for the two runs that were started to decide it.

Floating-point values of the relaxed problem for other degrees (μ = 0, extremals of `relaxed.py`, not certified):
R_B = 1.19210, 2.20844, 3.21814, 4.22504, 8.24100 for B = 1, 2, 3, 4, 8, that is R_B/B = 1.192, 1.104, 1.073,
1.056, 1.030. The bound per unit of degree tends to Faddeev's as B grows.

For orientation (not results of this note): Faddeev's bound for degree two is 2; the degree-two note certifies
I₂ ≤ 2.358835 (μ = 0) and I₂(1) ≤ 2.721257 with an axially symmetric trial map; the numerical value of
Battye and Sutcliffe [BS] for μ = 0 is 2.3582 ± 0.002. The bound of this note would leave a gap of about 6.4 %
(μ = 0) and 4.1 % (μ = 1) to those upper bounds.

## 7. Checks

Carried out where the log says so:

* B = 1: the new function returns the same lower bound, the same number of transitions and the same number of
  fallbacks as `levelSetLowerBound` on the stored grids.
* Path test: at every level the certified function, evaluated at the volume of the floating-point extremal, must
  not exceed the cost of that extremal beyond the level (an admissible path, so its cost is an upper bound of
  W_B there). At the first level this says that the certified number is below the floating-point value R_B.
* Sanity constraints: the certified number for B = 2 is at most 2R₁ (Proposition E), at most twice the certified
  hedgehog upper bound, at most the certified upper bounds for I₂ quoted above, and at least 2 minus the part
  of Faddeev's integral below τ₀ and above τ_N.
* Stored claims are recomputed from the stored grids and compared.

Not available: a second, independently written certified programme for degree two. The degree-one programme has
one (the Python verifier of the level-set note, reproduced to 3×10⁻¹² by the paclet); the degree-B function
differs from the degree-one function in the five places listed in Section 5. The same five changes in that
verifier (the factor B in its three step integrals of c and in the two values of c used for heuristics), run on
the grids of `Data/SkyrmeDegree`, would give the independent check; it reads the same grid files.

## 8. What an existence argument for degree three would need

By Esteban's theorem [E, Theorem 1], I_k is attained if I_k < Σ I_{d_i} for every decomposition of k into at
least two nonzero integers. For k = 3 the decompositions with a negative entry have Σ|d_i| ≥ 5, so Faddeev's
bound disposes of them as soon as I₃ < 5·12π². There remain 1 + 2 and 1 + 1 + 1, and since sub-additivity is not
known both must be checked: I₃ < I₁ + I₂ and I₃ < 3I₁. With certified lower bounds this needs a certified upper
bound for I₃ below min(LB(I₁) + LB(I₂), 3 LB(I₁)).

In units of 12π², for μ = 0: LB(I₁) = 1.192040 (Publication grid), so 3 LB(I₁) = 3.5761. With the value of this
note, LB(I₁) + LB(I₂) is about 1.1920 + 2.2084 = 3.4004. The numerical minimum for degree three is 3.4386 with
a quoted uncertainty of ±0.003 [BS, Table: E/B = 1.1462 ± 0.001]. So the second condition is within reach of a
good trial map (margin 0.137), and the first is not: the sum of the lower bounds is about 0.038 *below* the
numerical value of I₃, before any slack of a trial map. To close it the two lower bounds together would have to
improve by more than 0.038, for example the bound for I₂ from 2.208 to above 2.247 (it would then be within
4.7 % of I₂ instead of 6.4 %). Even the exact value of I₁ in place of its lower bound (1.2315 if the hedgehog is
the minimiser) would leave 1.2315 + 2.2084 − 3.4386 = 0.0013, which no realistic trial map can use. The level-set
method as it stands does not give existence in degree three.

## 9. The calibration method of the vector-field programme for Skyrme fields

That method bounds the energy from below by a null Lagrangian plus a pointwise supremum. For Skyrme fields on ℝ³
its reach is known and it is short: the companion note of the level-set note on topological bounds proves that
every lower bound obtained from a first-order null Lagrangian with continuous coefficients (a calibration on
ℝ³ × S³, which includes polyconvexity plus divergence terms and the occupation-measure and sum-of-squares
relaxations) is at most |B|·B_V, where B_V = 12π² for μ = 0 (Faddeev's bound) and 1.17003·12π² for μ = 1; the
mechanism is that translated and rotated copies of one map visit every point of the graph space with all
orientations, and no function that is affine in the minors can touch the energy on a whole orbit. The convexified
problem is therefore linear in B and sees no binding energy. I read the two-page proof and did not re-derive it.
If it is right, the calibration method cannot exceed 2·12π² for degree two when μ = 0, which is the value the
stop rule of the task treats as a negative result; the level-set bound is above it because the isoperimetric
inequality is not a pointwise inequality between the energy density and a null Lagrangian.

## 10. Not claimed

* That the bound is close to I₂. It is about 6 % below the numerical value (μ = 0).
* That 2 × (hedgehog energy) is a competitor for the degree-two minimiser. It is the trivial upper bound from two
  separated Skyrmions, used because it can be certified inside the paclet.
* Existence or any property of a minimiser of degree two or three.
* Any inequality between I₂ and 2I₁. (The degree-two note proves I₂ < 2I₁ with a trial map that is not part of
  this repository.)
* That the facts (F1)–(F3) were re-derived; they are quoted.

## References

* [E] M. J. Esteban, *Existence of 3D Skyrmions. Erratum, complete version*, arXiv:math-ph/0401042 (2004);
  Comm. Math. Phys. 251 (2004) 209–210. Theorem 1 (the existence criterion over all decompositions) and the remark
  that the non-strict inequalities I_k ≤ I_ℓ + I_{k−ℓ} are not known were read in the arXiv text.
* [BS] R. A. Battye, P. M. Sutcliffe, *Skyrmions, fullerenes and rational maps*, Rev. Math. Phys. 14 (2002)
  29–85, arXiv:hep-th/0103026. The table of energies was read in the arXiv text: E/B = 1.2322, 1.1791, 1.1462
  for B = 1, 2, 3, quoted error ±0.001 in E/B.
* [MSZ] J. Malý, D. Swanson, W. P. Ziemer, *The co-area formula for Sobolev mappings*, Trans. Amer. Math. Soc. 355
  (2003) 477–492. (Quoted through the level-set note; not read.)
* [EG] L. C. Evans, R. F. Gariepy, *Measure Theory and Fine Properties of Functions*, CRC Press, 1992. (Textbook
  facts; quoted.)
* The level-set note and the degree-two note: references [2] and [1] of the README.

**What was searched (23 November 2026).** Three web searches: for the energy table of [BS] and the statement of
[E] (both then read in the arXiv texts), and for rigorous lower bounds for Skyrme fields of degree two beyond
Faddeev's. The last one returned only Faddeev's bound and Harland's bound for massive pions (Phys. Lett. B 728
(2014) 518), which are linear in the degree. This is not a literature review.

## Log

* 23 November 2026. Sources located and read (see above). Audit written: the degree enters through the flux
  only; Theorems 1 and 2 with proofs; the end-point argument for the floors added (Lemma C). Floating-point
  extremals: R₂ = 2.20844 (μ = 0), 2.61111 (μ = 1); the stop rule (within 1 % of 2) does not fire. Grid generator
  written and compared with the stored degree-one grid. `Kernel/SkyrmeLevelSetDegree.wl` written.
* 23 November 2026, later. First certified run. B = 1 on the stored grid `Data/Skyrme/mu0_quick`: the new function
  and `levelSetLowerBound` both return 1.1917480092669916 with 4,917,322 transitions and 1,457,421 fallbacks
  (identical). B = 2, μ = 0, Quick grid: certified 2.207981318679515, floating-point value of the relaxed problem
  2.2084398; path test passed at all 1011 levels (smallest slack 3.4×10⁻⁶·12π², at the last level; largest
  4.5×10⁻⁴·12π², at the first). Model `SkyrmeDegreeTwo` registered, test suite `Tests/SkyrmeDegree.wlt` and
  `tools/skyrme_degree_claims.wls` written; grids for μ = 1 (Quick) and μ = 0 (Standard) generated.
  **Done:** audit and proofs; float values; generator; kernel file; B = 1 identity; first certified number.
  **In progress:** the run that stores the claims for the two Quick grids (μ = 0 and μ = 1).
  **Not yet run:** `Tests/SkyrmeDegree.wlt` as a whole; the existing suites on this branch; the Standard grids.
  **Next:** fill the table of Section 6 and the README paragraph from the stored claims.
* 23 November 2026, at the stop (work in progress).
  **State of the audit: complete.** Every step of the degree-one argument has been read in the source and
  classified in the table of Section 2: the pointwise lemma, the coarea step, the flux, the volume, the
  isoperimetric step, the minorants (M1)–(M3), monotonicity in the volume, the three floors, the second-order
  step bound, the backward induction, and the extension to finite-energy Sobolev maps with the degree formula.
  The degree is used in the flux only. **Verdict: the theorem extends to every degree; no obstruction**
  (Theorems 1 and 2, proved in Sections 3 and 4). Not re-derived: the facts (F1)–(F3), which are quoted. Of the
  code, only the lines that contain the flux were compared with the proof; the rest is the degree-one function.
  **Done since the last entry:** stored claims for the two Quick grids (`Data/SkyrmeDegree/claims.json`):
  μ = 0: computed 2.207981318679515, claim 2.207981316; μ = 1: computed 2.572283249253187, claim 2.572283246;
  path test passed for both (smallest slack 3.4×10⁻⁶ and 3.2×10⁻⁶ in units of 12π²). Standard grids generated for
  both masses. Floating-point values for B = 3, 4, 8.
  **Stopped unfinished:** two runs to settle the open point for μ = 1 (Section 6): the level-by-level comparison
  of the certified function with the cost of the extremal (run `levelSetLowerBoundDegree` with `"Trace" -> True`
  on `b2_mu1_quick` and print `ctg` minus the trace against the level, to see where the 1.49 % arise), and a
  direct minimisation over discretised two-phase paths started from several initial paths.
  **Not run at all:** `Tests/SkyrmeDegree.wlt` as a whole, and in particular `CertifiedVariationalBound` and
  `VerifyCertificate` for the model `SkyrmeDegreeTwo` (the pieces run so far are `levelSetLowerBoundDegree` for
  B = 1 on `Data/Skyrme/mu0_quick` and for B = 2 on the two Quick grids, and the claims tool with `--write`);
  the existing suites (fast, `--nd`, `--skyrme`) on this branch after the edits to `Kernel/Models.wl` and
  `Kernel/CertiField.wl`; the Standard grids. The static syntax check of all sources passes.
  **Not written:** the README paragraph.
  **Next, in this order:**
  1. `wolframscript -file tools/run_tests.wls --skyrme-degree`, and repair what fails. One failure is expected:
     the test `degree-two-mu1-below-and-near-the-relaxed-value` asks for the certified value within 0.002 of
     2.6111 and the certified value is 2.5723; do not change the test before step 2 is decided.
  2. Settle the open point for μ = 1: the trace by level; the Standard grid
     (`tools/skyrme_degree_claims.wls --write b2_mu1_standard`); the direct minimisation. If a cheaper path
     exists, `relaxed.py` needs the correct extremal for μ > 0 and the grids for μ = 1 must be regenerated around
     it.
  3. The existing suites: fast, `--nd`, `--skyrme`.
  4. Claims for the Standard grids; in `tools/skyrme_degree_claims.wls` the numbers are printed inside
     `NumberForm[…]` wrappers, which should be replaced by `ToString[…, InputForm]`.
  5. Complete the table of Section 6, then the README paragraph (plain product text, numbers from `claims.json`).
