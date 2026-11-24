# Magnetisation in thin films: certified bounds for unit-vector fields

This note states and proves what the magnetisation track certifies and describes how the two bounds are computed.
It is the analytic reduction behind `Reference/Magnetisation` (the Python reference).

**State of the work (23 November 2026).** Sections 1 to 5 (statements and proofs) are written. The reference
implementation and the tests are in progress; Sections 7 and 8 are filled in as they are completed.
`Reference/Magnetisation/LOG.md` is the running log.

**Status.** The method is a calibration (null-Lagrangian) argument. Its ingredients are classical: the degree of
maps into the sphere, the Belavin–Polyakov inequality, the helical-derivative factorisation used for chiral
skyrmions, and Lin's null Lagrangian for x/|x|. No novelty is claimed for any of them, nor for their combination.
The short proofs below have been checked by their author only. Numbers quoted as *certified* come from exact
rational arithmetic or outward-rounded interval arithmetic on exact inputs; numbers marked [FLOAT] are ordinary
floating point.

## 1. The problems

**Thin film.** Ω ⊂ ℝ² is a bounded open set with Lipschitz boundary (a rectangle in the implementation), A > 0, and
κ, K, h are real numbers (rational in the implementation). For m : Ω → S² ⊂ ℝ³,

    E[m] = ∫_Ω A|∇m|² + κ m·(∇×m) + K(1 − m₃²) + h(1 − m₃) dx,      ∇×m = (∂₂m₃, −∂₁m₃, ∂₁m₂ − ∂₂m₁),

    𝒜 = { m ∈ H¹(Ω; ℝ³) : |m| = 1 a.e., m = e₃ on ∂Ω (trace) },      𝒜_Q = { m ∈ 𝒜 : Q(m) = Q },

    Q(m) = (1/4π) ∫_Ω m·(∂₁m × ∂₂m) dx.

Q(m) is an integer (Lemma 1). The certified statement is about inf { E[m] : m ∈ 𝒜_Q } for a given integer Q. With
p = ∂₁m, q = ∂₂m the integrand is f(m, p, q) = A(|p|² + |q|²) + d₁(m)·p + d₂(m)·q + V(m₃), where
d_i(y) = κ y × e_i (a tangent vector at y) and V(t) = K(1 − t²) + h(1 − t).

*Other Dzyaloshinskii–Moriya terms.* If e₁, e₂ in d_i are replaced by a_i = R e_i, R a rotation about e₃ (interfacial
coupling is the rotation by a right angle), then E_R[m] = E[R⁻¹m], and m ↦ R⁻¹m maps 𝒜_Q onto itself. So the infimum
is the same and every bound below applies unchanged.

**Strip (domain wall).** Ω = (−L, L) × (0, w), E as above with κ = 0, m = −e₃ on {x₁ = −L}, m = e₃ on {x₁ = L},
nothing prescribed on the other two sides.

**Cube (hedgehog).** Ω = (−1, 1)³, E[m] = A ∫_Ω |∇m|², m ∈ H¹(Ω; S²) with the trace of x/|x| on ∂Ω.

## 2. Three lemmas

**Lemma 1 (degree).** Let m ∈ 𝒜. Then Q(m) is an integer, and for every bounded Borel function ρ : S² → ℝ

    ∫_Ω ρ(m) m·(∂₁m × ∂₂m) dx = Q(m) ∫_{S²} ρ dℋ².

*Proof.* Extend m by e₃ to ℝ²; the extension has locally square-integrable gradient because m − e₃ ∈ H¹₀(Ω). Let
σ : S² ∖ {N} → ℝ² be an orientation-preserving stereographic projection and ũ = m ∘ σ, ũ(N) = e₃. Then ũ = e₃ near
N, ũ ∈ H¹(S²; S²) (the Dirichlet integral is conformally invariant in two dimensions), and by the change of
variables ∫_Ω ρ(m) m·(∂₁m × ∂₂m) dx = ∫_{S²} ρ(ũ) ũ*ω, where ω is the area form and ũ*ω = ũ·(∂₁ũ × ∂₂ũ) in
oriented coordinates. Smooth maps are dense in H¹(S²; S²) (Schoen and Uhlenbeck [9]). For a smooth u and
continuous ρ, ∫ ρ(u) u*ω = ∫ u*(ρω) = deg(u) ∫ ρω with deg(u) ∈ ℤ. Let u_k be smooth, u_k → ũ in H¹ and almost
everywhere. Then ∂₁u_k × ∂₂u_k → ∂₁ũ × ∂₂ũ in L¹, and multiplying by the bounded, almost everywhere convergent
factors u_k and ρ(u_k) keeps convergence in L¹. So ∫ ρ(u_k) u_k*ω → ∫ ρ(ũ) ũ*ω. With ρ = 1 the integers deg(u_k)
converge to Q(m), hence are eventually equal to it, and the formula follows for continuous ρ. Both sides are
integrals of ρ against finite signed Borel measures on S²; measures that agree on continuous functions are equal. ∎

For maps that are constant on ∂Ω this lemma contains the families (T) and (Q) of the general method with
coefficients that do not depend on x: adding c(m) m·(∂₁m × ∂₂m) to the integrand changes the energy of every
m ∈ 𝒜_Q by the same number Q ∫ c.

**Lemma 2 (divergences).** Let S : S² → ℝ² be the restriction of a C¹ map and m ∈ 𝒜. Then S(m) ∈ H¹(Ω; ℝ²),
div[S(m)] = ∇_τS₁(m)·∂₁m + ∇_τS₂(m)·∂₂m almost everywhere (∇_τ is the tangential gradient on the sphere), and
∫_Ω div[S(m)] dx = 0.

*Proof.* Extend S to a C¹ map on ℝ³ with bounded derivative. The chain rule for a C¹ function with bounded
derivative and an H¹ map gives S(m) ∈ H¹ and ∂_i[S_j(m)] = ∇S_j(m)·∂_i m. Since |m|² = 1, m·∂_i m = 0 almost
everywhere, so only the tangential part of ∇S_j enters. S(m) − S(e₃) is a Lipschitz function, vanishing at 0, of
m − e₃ ∈ H¹₀(Ω), hence lies in H¹₀(Ω), and the integral of the divergence of an H¹₀ field is zero. ∎

**Lemma 3 (the supremum over the gradient).** Let y ∈ S², T = y^⊥, J v = y × v for v ∈ T, A > 0, 0 ≤ u ≤ 4A,
ℓ₁, ℓ₂ ∈ T, and w₊ = ℓ₁ + Jℓ₂, w₋ = ℓ₁ − Jℓ₂. For all p, q ∈ T

    A(|p|² + |q|²) − (u − 2A) y·(p × q) − ℓ₁·p − ℓ₂·q
        = (u/4)|p + Jq|² − ½ w₊·(p + Jq) + ((4A − u)/4)|p − Jq|² − ½ w₋·(p − Jq).

Hence the infimum of the left side over p, q ∈ T is

    −|w₊|²/(4u) − |w₋|²/(4(4A − u))        if 0 < u < 4A,
    −|w₋|²/(16A) if u = 0 and w₊ = 0;   −|w₊|²/(16A) if u = 4A and w₋ = 0;   −∞ otherwise.

*Proof.* J is a rotation of T by a right angle and p·Jq = −y·(p × q). So |p ± Jq|² = |p|² + |q|² ∓ 2 y·(p × q),
and (u/4)|p + Jq|² + ((4A − u)/4)|p − Jq|² = A(|p|² + |q|²) − (u − 2A) y·(p × q). Also
½ w₊·(p + Jq) + ½ w₋·(p − Jq) = ℓ₁·p + Jℓ₂·Jq = ℓ₁·p + ℓ₂·q. The map (p, q) ↦ (p + Jq, p − Jq) is a bijection of
T × T, and the two squares are minimised separately. ∎

The condition |c| ≤ 2A on the coefficient c = u − 2A of the minor is the Belavin–Polyakov inequality
A(|p|² + |q|²) ≥ 2A |y·(p × q)|.

## 3. The lower bound

**Theorem 4.** Let S : S² → ℝ² be C¹, u : S² → [0, 4A] bounded and Borel, Γ ∈ ℝ. Put ℓ_i = ∇_τS_i − d_i,
w± = ℓ₁ ± Jℓ₂ and assume that at every y ∈ S² the infimum in Lemma 3 is at least −V(y₃) − Γ. Then for every
m ∈ 𝒜_Q

    E[m]  ≥  Q ∫_{S²} (u − 2A) dℋ²  −  Γ |Ω|.

*Proof.* At almost every x the vectors p = ∂₁m, q = ∂₂m are tangent at y = m(x), and by Lemma 3 and the
assumption

    f(y, p, q) − [ (u(y) − 2A) y·(p × q) + ∇_τS₁(y)·p + ∇_τS₂(y)·q ] ≥ −Γ.

All terms are integrable. Integrate over Ω and use Lemma 1 with ρ = u − 2A and Lemma 2. ∎

Any S, u, Γ that satisfy the assumption give a valid bound. If Γ = 0 the bound does not depend on Ω: it holds for
every bounded Lipschitz domain.

**Proposition 5 (fields that depend on y₃ and are equivariant).** Let φ ∈ C¹[−1, 1], S(y) = φ(y₃)(y₂, −y₁), and
u = u(y₃). With t = y₃,

    |w₊|² = 4P(t),   P = ¼ (1 − t)² [ ((1 + t)(φ + κ))′ ]²,        |w₋|² = 4R(t),   R = ¼ (1 + t)² [ ((1 − t)(φ − κ))′ ]².

Consequently the assumption of Theorem 4 holds if, on each interval of a finite partition of [−1, 1], u is a
polynomial with 0 ≤ u ≤ 4A and one of the following holds on that interval:

    (a) u is not identically 0 or 4A, and   Π_Γ := P(4A − u) + R u − (V + Γ) u (4A − u) ≤ 0;
    (b) u ≡ 0,  P ≡ 0  and  R ≤ 4A(V + Γ);
    (c) u ≡ 4A,  R ≡ 0  and  P ≤ 4A(V + Γ).

Then E[m] ≥ 2πQ ( ∫₋₁¹ u dt − 4A ) − Γ|Ω| for every m ∈ 𝒜_Q.

*Proof.* Let ρ² = 1 − t² > 0 and e_φ = e₃ × y = (−y₂, y₁, 0), τ = ∇_τ y₃ = e₃ − t y. These are orthogonal tangent
vectors of length ρ with J e_φ = τ, Jτ = −e_φ. One finds ∇_τ y₁ = (−y₂ e_φ − t y₁ τ)/ρ², ∇_τ y₂ = (y₁ e_φ − t y₂ τ)/ρ²,
and d₁ = κ J∇_τ y₁, d₂ = κ J∇_τ y₂. With α = (φ − κt)/ρ² and β = (φ′ρ² − φt + κ)/ρ²,

    ℓ₁ = α y₁ e_φ + β y₂ τ,      ℓ₂ = α y₂ e_φ − β y₁ τ,      Jℓ₂ = β y₁ e_φ + α y₂ τ,

so w₊ = (α + β)(y₁ e_φ + y₂ τ) and w₋ = (α − β)(y₁ e_φ − y₂ τ), with |w±|² = (α ± β)² ρ⁴. Finally
(α + β)ρ² = (1 − t)((1 + t)(φ + κ))′ and (α − β)ρ² = −(1 + t)((1 − t)(φ − κ))′. At the poles the formulas hold by
continuity. For the second statement: where 0 < u < 4A, (a) divided by u(4A − u) is the condition of Theorem 4. At
a zero t₀ of u in an interval of type (a), Π_Γ(t₀) = 4A P(t₀) ≤ 0 gives w₊ = 0, and R/(4A − u) ≤ V + Γ holds at
nearby points of the interval where u > 0, hence at t₀ by continuity; this is the condition for u = 0. Zeros of
4A − u are treated in the same way, and (b), (c) are the conditions of Lemma 3 directly. On a common endpoint of
two intervals either polynomial may be used to define u. The value of the bound is Theorem 4 with
∫_{S²} g(y₃) dℋ² = 2π ∫₋₁¹ g(t) dt. ∎

For piecewise polynomial φ and u with rational coefficients, conditions (a) to (c) are polynomial inequalities
with rational coefficients, and ∫ u dt is a rational number. The certified lower bound is therefore an exact
rational multiple of π minus Γ|Ω|.

**Three checks of Proposition 5 against known results** (Γ = 0):

* κ = 0, V = 0, φ = 0, u ≡ 4A (case (c)) or u ≡ 0 (case (b)): E ≥ 8πA|Q| for Q > 0 and Q < 0 respectively. This is
  the Belavin–Polyakov bound [3]; on a bounded domain it is the infimum and is not attained.
* A = ½, κ = r ∈ (0, 1], V = ½(1 − t)² (K = −½, h = 1), φ ≡ r, u ≡ 2r²: P = r²(1 − t)², R = 0,
  Π₀ = (2 − u)(1 − t)²(r² − u/2) = 0. Hence E ≥ 4π(1 − 2r²) in 𝒜₋₁. On the plane this is the exact infimum, attained
  by a Belavin–Polyakov map of radius 2r (Döring and Melcher [5]; all degrees in [6]). For r = 1 it is the
  Bogomolny bound of the critical-coupling model [4].
* K = 0, h ≥ κ²/(2A), φ ≡ κ, u = κ²(1 − t)/h: P = κ²(1 − t)², R = 0, Π₀ = 0. Hence E ≥ 8πA − 4πκ²/h in 𝒜₋₁.

The optimal φ, u for the Zeeman model are not of this simple form; Section 6 describes how they are found, and
Section 8 the resulting bounds.

**Why fields independent of x.** On ℝ² the right-hand side of the general bound is concave in the dual fields and
invariant under translations of x and under simultaneous rotations of x and m about e₃. Averaging shows that
fields independent of x and equivariant are optimal there. On a bounded domain with m = e₃ on the boundary,
x-dependent fields could gain only through the boundary. They are not used here.

**Remark (no exact calibration).** For the axisymmetric skyrmion m₀ of the Zeeman model, equality in Theorem 4
along m₀ would require u − 2A = −(A|∇m₀|² − V)/|m₀·(∂₁m₀ × ∂₂m₀)|, which is below −2A in the tail of the
skyrmion. So the bound of Proposition 5 cannot be attained by m₀, and whether m₀ minimises E in 𝒜₋₁ remains the
open question stated in [1].

## 4. The strip and the cube

**Lemma 6 (polar angle).** For m ∈ H¹(Ω; S²), ϑ = arccos(m₃) belongs to H¹(Ω) and |∂_iϑ| ≤ |∂_i m| almost everywhere.

*Proof.* Θ(y) = arccos(y₃) is the geodesic distance on S² from e₃, so |Θ(y) − Θ(y′)| ≤ d_g(y, y′) ≤ (π/2)|y − y′|
and Θ extends to a Lipschitz function on ℝ³; hence ϑ = Θ(m) ∈ H¹(Ω). On almost every line parallel to a coordinate
axis m is absolutely continuous, and for two points s < s′ on it |ϑ(s) − ϑ(s′)| ≤ d_g(m(s), m(s′)) ≤ ∫_s^{s′} |∂_i m|,
because the geodesic distance is at most the length of the curve m traces on the sphere. ∎

**Proposition 7 (strip).** For the strip problem with V ≥ 0, inf E = w · e(L), where
e(L) = inf { ∫_{−L}^{L} Aθ′² + V(cos θ) dx : θ(−L) = π, θ(L) = 0 }. Moreover, if s is a C¹ piecewise polynomial
on [−1, 1] with (1 − t²) s′(t)² ≤ 4A V(t), then E[m] ≥ w (s(1) − s(−1)) for every admissible m.

*Proof.* By Lemma 6, E[m] ≥ ∫_Ω A|∂₁ϑ|² + V(cos ϑ) dx, and for almost every x₂ the function ϑ(·, x₂) is admissible
for e(L); this gives ≥. The map (sin θ(x₁), 0, cos θ(x₁)) gives ≤. For the second statement,
∂₁[s(m₃)] = s′(m₃) ∂₁m₃ and |∂₁m₃| ≤ (1 − m₃²)^{1/2} |∂₁m|, so
∂₁[s(m₃)] − A|∂₁m|² ≤ s′(m₃)²(1 − m₃²)/(4A) ≤ V(m₃); integrate, using
∫_Ω ∂₁[s(m₃)] dx = w (s(1) − s(−1)) (Fubini and the boundary values). ∎

For V = K(1 − t²) and s(t) = 2σt with σ² ≤ AK this gives E ≥ 4σw, the wall energy of the infinite line; it is below
w · e(L) by a factor that is exponentially close to one in L(K/A)^{1/2}. The number e(L) is what the existing
one-dimensional model `MagneticDomainWall` certifies.

**Lemma 8 (cube; Lin's null Lagrangian).** Let Ω ⊂ ℝ³ be a bounded Lipschitz domain and
N₀(F) = (tr F)² − tr(F²). (a) If u, v ∈ H¹(Ω; ℝ³) and u − v ∈ H¹₀, then ∫ N₀(∇u) = ∫ N₀(∇v). (b) If |u| = 1 almost
everywhere, then N₀(∇u) ≤ |∇u|² almost everywhere. (c) u₀ = x/|x| ∈ H¹(Ω; S²) and N₀(∇u₀) = |∇u₀|² = 2/|x|².
Consequently u₀ minimises ∫|∇m|² among maps in H¹(Ω; S²) with its boundary values, and for Ω = (−1, 1)³

    min E = 2A ∫_Ω |x|⁻² dx = 48 A J,      J = ∫₀¹∫₀¹ ds dt / (1 + s² + t²).

*Proof.* (a) N₀(∇u) − N₀(∇v) = B(∇(u − v), ∇(u + v)) with B(F, G) = tr F tr G − tr(FG). For smooth φ with compact
support and smooth w, B(∇φ, ∇w) = ∂_i(φ_i ∂_j w_j) − ∂_j(φ_i ∂_i w_j), whose integral vanishes; B is continuous on
L² × L², and smooth functions are dense. (b) At almost every point uᵀ∇u = 0. In an orthonormal basis (n₁, n₂, u)
the matrix of ∇u has third row zero; with upper left block ((a, b), (c, d)), N₀ = 2(ad − bc) ≤ a² + b² + c² + d².
(c) Direct computation; 2/|x|² is integrable in three dimensions. Then A∫|∇m|² ≥ A∫N₀(∇m) = A∫N₀(∇u₀) = A∫|∇u₀|².
The divergence theorem for x/|x|² on Ω minus a small ball gives ∫_Ω |x|⁻² = ∫_{∂Ω} x·ν/|x|² = 24 J. ∎

The minimality of x/|x| is due to Brezis, Coron and Lieb [7] (for the ball); the proof by this null Lagrangian is
Lin's [8]. In the terms of the general method, N₀ is the family (T) with T_a = u × e_a at the critical size of the
coefficient. The certified statement for the cube is an enclosure of 48AJ by rigorous quadrature of the rational
function 1/(1 + s² + t²); both bounds equal it.

## 5. The upper bound: a normalised interpolant

Let M_h be the continuous function that is bilinear on every cell of a tensor-product grid on the rectangle, with
nodal values M_ij ∈ ℝ³ equal to (0, 0, 1) at the nodes of ∂Ω (for the strip: (0, 0, ∓1) on the two ends), and
suppose D = |M_h|² ≥ D_min > 0 on Ω. Then m_h = M_h/|M_h| is Lipschitz, has the required boundary values exactly,
and

    |∂_i m_h|² = |M_h × ∂_i M_h|²/D²,      (m_h × a)·∂_i m_h = −a·(M_h × ∂_i M_h)/D,
    m_h·(∂₁m_h × ∂₂m_h) = det[M_h, ∂₁M_h, ∂₂M_h]/D^{3/2},      m_{h,3} = M_{h,3}/D^{1/2}.

On a cell with corner values B₀₀, B₁₀, B₀₁, B₁₁ and local coordinates (ξ, η) ∈ [0, 1]²:

* M_h × ∂_ξM_h does not depend on ξ: it is v₀(1 − η)² + 2v₁η(1 − η) + v₂η² with v₀ = B₀₀ × B₁₀,
  v₁ = ½(B₀₀ × B₁₁ + B₀₁ × B₁₀), v₂ = B₀₁ × B₁₁; likewise M_h × ∂_ηM_h depends on ξ only.
* D is a polynomial of degree two in each variable whose nine Bernstein coefficients are B₀₀·B₀₀, B₀₀·B₁₀,
  B₁₀·B₁₀, B₀₀·B₀₁, ½(B₀₀·B₁₁ + B₁₀·B₀₁), …; their minimum a and maximum b bound D on the cell.
* det[M_h, ∂_ξM_h, ∂_ηM_h] is bilinear, with corner values det[B₀₀, B₁₀, B₀₁], det[B₀₀, B₁₀, B₁₁],
  det[B₀₀, B₁₁, B₀₁], det[B₁₀, B₁₁, B₀₁].

The weights 1/D², 1/D, 1/D^{1/2} are convex functions of D. On [a, b] they are bounded by linear functions of D
with rational coefficients:

    1/D² ≤ (a² + ab + b² − (a + b)D)/(a²b²),     1/D ≤ (a + b − D)/(ab),     1/D ≥ 2/d₀ − D/d₀²  (any d₀ > 0),
    (3 − D)/2 ≤ 1/D^{1/2} ≤ ½(1 + (a + b − D)/(ab)).

With these, each term of the energy on a cell is bounded above by the integral of a polynomial, which is an exact
combination of the corner values (products of Bernstein polynomials integrate to ratios of binomial
coefficients), plus a remainder of second order in b − a for the terms without a sign. `upper_bound` evaluates
these expressions with outward-rounded interval arithmetic. The error of the bound is of fourth order in the
angle between neighbouring nodal vectors.

**Degree of the competitor.** By Lemma 1, Q(m_h) is an integer. `degree` encloses ∫ det[…]/D^{3/2} cell by cell
(the range of the bilinear numerator is the hull of its corner values; 2/(b(1 + b)) ≤ D^{−3/2} ≤ (1 + a)/(2a²)) and
accepts the competitor only if the enclosure, divided by an enclosure of 4π, contains exactly one integer.

For the strip there is no degree, and the same energy bound is used with the boundary values of the strip.

## 6. Choosing the fields (heuristic)

Nothing in this section is trusted. The pair (φ, u) is found by minimising ½∫u dt subject to the pointwise
condition of Proposition 5 in floating point (a dynamic programme in the variables ln(1 − t) and (1 + t)(φ + κ),
followed by a fit of a C¹ piecewise cubic φ and of a piecewise polynomial u slightly above the smaller root of
Π₀ = 0), and is then rounded to rational data and checked exactly. The competitor is the interpolant of a
numerically computed axisymmetric profile on a graded grid.

## 7. Checks

To be completed with the implementation.

## 8. Results

To be completed with the implementation. Floating-point values of the optimised dual of Proposition 5 for the
Zeeman model are in `Reference/Magnetisation/LOG.md`.

## 9. What is not claimed, and what is not covered

* **Not claimed:** that the competitor is a minimiser, or that the axisymmetric skyrmion is the minimiser in its
  degree class. The certified interval is for the infimum of E over 𝒜_Q.
* The stray field (magnetostatic energy) is not included.
* Dual fields that depend on x are not used. For domains that are small compared with the skyrmion the bound of
  Proposition 5 is valid but not tight.
* Boundary values other than e₃ (thin film) and ±e₃ on two opposite sides (strip); free boundaries for the thin
  film with a Dzyaloshinskii–Moriya term (the term then contributes on the boundary).
* Three dimensions beyond the exchange-only cube.

## References

1. X. Li, C. Melcher, *Stability of axisymmetric chiral skyrmions*, J. Funct. Anal. 275 (2018) 2817–2844,
   arXiv:1711.07717. (The Zeeman model; journal data from memory, the arXiv text was read.)
2. C. Melcher, *Chiral skyrmions in the plane*, Proc. R. Soc. A 470 (2014) 20140394. (Not read; cited through [1]
   and [6].)
3. A. A. Belavin, A. M. Polyakov, *Metastable states of two-dimensional isotropic ferromagnets*, JETP Lett. 22
   (1975) 245–247. (From memory.)
4. B. Barton-Singer, C. Ross, B. J. Schroers, *Magnetic skyrmions at critical coupling*, Comm. Math. Phys. 375
   (2020) 2259–2280, arXiv:1812.07268.
5. L. Döring, C. Melcher, *Compactness results for static and dynamic chiral skyrmions near the conformal limit*,
   Calc. Var. Partial Differential Equations 56 (2017) 60, arXiv:1611.01984. (The minimum 4π(1 − 2r²); quoted
   from the introduction of [6].)
6. S. Ibrahim, T. Miura, C. Román, I. Shimizu, *Phase transition thresholds and chiral magnetic fields of general
   degree*, arXiv:2512.24598. (Read: Section 1, Theorem 1.1.)
7. H. Brezis, J.-M. Coron, E. H. Lieb, *Harmonic maps with defects*, Comm. Math. Phys. 107 (1986) 649–705.
8. F.-H. Lin, *Une remarque sur l'application x/|x|*, C. R. Acad. Sci. Paris Sér. I 305 (1987) 529–531. (Not
   read; title and pages from memory; later papers found in the search credit it with the null-Lagrangian proof.)
9. R. Schoen, K. Uhlenbeck, *Boundary regularity and the Dirichlet problem for harmonic maps*, J. Differential
   Geom. 18 (1983) 253–268. (Density of smooth maps in H¹(M²; N); from memory, not confirmed by the search.)
10. S. Komineas, C. Melcher, S. Venakides, *The profile of chiral skyrmions of small radius*, Nonlinearity 33
    (2020) 3395–3408, arXiv:1904.01408.
11. A. Chernyavsky, J. J. Bramburger, G. Fantuzzi, D. Goluskin, *Convex relaxations of integral variational
    problems: pointwise dual relaxation and sum-of-squares optimization*, arXiv:2110.03079. (Not compared yet.)

**What was searched (23 November 2026) and what was found.** Web searches for the papers above, and reading of the
arXiv texts of [1] and [6]. Found: the model and the open question in [1]; the exact minimal energies of the
critical family in [6], which Proposition 5 reproduces with constant φ and u; the statement in [10] that a
rigorous lower bound matching the small-radius asymptotics was not available there. Not found: a computed lower
bound for the energy of chiral skyrmions in the Zeeman model at moderate fields. That search is not a literature
review, and the absence of a result in it is not evidence that none exists.
