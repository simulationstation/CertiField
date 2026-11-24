# Skyrme model, degree B: certified upper bounds from the rational-map ansatz

This note states and proves what the model `SkyrmeRationalMap` certifies: an upper bound for the infimum I_B of the
Skyrme energy over maps of degree B, obtained as the energy of one explicit map built from a rational map of the
Riemann sphere and a radial profile (the ansatz of Houghton, Manton and Sutcliffe [1]).

**Status.** The ansatz, its energy formula and the lower bound used here are classical; no novelty is claimed for
any of them. The proofs below are short and have been checked by their author only. What the package adds is a
certified enclosure of the two integrals. Numbers quoted as *certified* come from outward-rounded and exact
rational arithmetic; numbers marked (float) are ordinary floating point; numbers marked (published) are quoted
from the literature.

## 1. The problem

Units as elsewhere in the paclet: for U : ℝ³ → S³ and a pion mass μ ≥ 0,

    E_μ(U) = ∫_{ℝ³} |dU|² + |Λ²dU|² + 2μ² (1 − cos f) dx,        f = dist_{S³}(U, 1),

reported as E_μ/12π². For an integer B, 𝒞_μ^B is the class of weakly differentiable maps with E_μ(U) < ∞,
U − 1 ∈ L⁶(ℝ³) and degree (1/2π²)∫ det dU = B, and I_B(μ) = inf { E_μ(U) : U ∈ 𝒞_μ^B }.

## 2. The competitor

Let R = p/q be a rational map of degree B ≥ 1: p and q are polynomials in z without a common zero and
max(deg p, deg q) = B. Through stereographic projection, z = tan(θ/2) e^{iφ} is a point of the unit sphere, and R
defines a smooth map n_R : S² → S² of degree B. Let f : [0, ∞) → ℝ be Lipschitz with f(0) = π and f(r) = 0 for
r ≥ r_max. The competitor is

    U(x) = exp( i f(|x|) n_R(x/|x|)·σ ) = ( cos f, sin f · n_R ) ∈ S³ ⊂ ℝ⁴,        U(0) = −1.

**Theorem 1 (Houghton–Manton–Sutcliffe [1]; classical).** U ∈ 𝒞_μ^B for every μ ≥ 0, and

    E_μ(U) = 4π ∫₀^∞ ( r² f′² + 2B (f′² + 1) sin² f + ℐ sin⁴ f / r² + 2μ² r² (1 − cos f) ) dr,
    ℐ = ℐ(R) = (1/4π) ∫_{S²} J⁴ dA,        J(z) = (1 + |z|²) |R′(z)| / (1 + |R(z)|²).

J is the factor by which the differential of n_R stretches lengths (the same in every tangent direction, because R is
holomorphic).

*Proof.* n_R is smooth and f is Lipschitz, so U is locally Lipschitz away from the origin; near the origin
|sin f| ≤ C|x|, so U is Lipschitz on ℝ³, and U = 1 for |x| ≥ r_max. At a point x ≠ 0 where f is differentiable let
e_r be the radial unit vector and e₁, e₂ orthonormal tangent vectors of the sphere of radius r = |x|. Then

    dU(e_r) = f′ · (−sin f, cos f · n_R),        dU(e_j) = (0, sin f · dn_R(e_j)),    |dn_R(e_j)| = J/r.

The vector (−sin f, cos f·n_R) has length one and is orthogonal to (0, w) for every w ⊥ n_R; dn_R(e₁) ⊥ dn_R(e₂)
because n_R is conformal. So the singular values of dU are |f′|, sin f·J/r, sin f·J/r, and

    |dU|² = f′² + 2 sin² f J²/r²,        |Λ²dU|² = 2 f′² sin² f J²/r² + sin⁴ f J⁴/r⁴.

Integrate in polar coordinates. J² dA is the pull-back of the area form of the target sphere under the
orientation-preserving map n_R, so ∫_{S²} J² dA = 4πB; and ∫ J⁴ dA = 4π ℐ by definition. This gives the formula. All
terms are bounded and vanish for r ≥ r_max, so the energy is finite and U − 1 ∈ L⁶. The determinant of dU in
oriented frames is −f′ sin² f J²/r² up to the fixed sign of the orientation convention, so

    deg U = ± (1/2π²) ∫₀^∞ (−f′) sin² f dr · ∫_{S²} J² dA = ± (1/2π²) · (π/2) · 4πB = ± B,

where ∫₀^∞ (−f′) sin² f dr = [f/2 − sin(2f)/4] between f = 0 and f = π, whatever the shape of f. The sign is that of
the hedgehog (R = z, B = 1), which has degree one in the convention of the paclet. ∎

**Corollary.** I_B(μ) ≤ E_μ(U) for every such R and f.

## 3. The lower end of the interval

**Lemma 2 (Faddeev; classical).** E_μ(U) ≥ 12π² |B| for every U ∈ 𝒞_μ^B.

*Proof.* If λ₁, λ₂, λ₃ ≥ 0 are the singular values of dU at a point, the inequality of the arithmetic and geometric
means for six numbers gives |dU|² + |Λ²dU|² = Σλᵢ² + Σ_{i<j} λᵢ²λⱼ² ≥ 6 λ₁λ₂λ₃ = 6 |det dU|. Integrating,
E_μ(U) ≥ 6 ∫ |det dU| ≥ 6 |∫ det dU| = 12π² |B|. ∎

A certificate of this model therefore asserts B ≤ I_B(μ)/12π² ≤ u, with u the certified upper bound. **These
intervals are wide**: the lower end is the topological bound, which no map attains. For B = 1 the paclet has a much
better lower bound (the model `SkyrmeDegreeOne`); for B ≥ 2 it has none.

## 4. The angular integral ℐ

**As a plane integral.** With W = p′q − pq′, |R′| = |W|/|q|² and 1 + |R|² = (|p|² + |q|²)/|q|², and
dA = 4 dx dy/(1 + |z|²)²:

    ℐ = (1/π) ∫_{ℝ²} F_{p,q}(z) dx dy,        F_{p,q}(z) = (1 + |z|²)² |W(z)|⁴ / (|p(z)|² + |q(z)|²)⁴.

The denominator is positive because p and q have no common zero.

**Four pieces over one rectangle.** Put p̃(u) = u^B p(1/u), q̃(u) = u^B q(1/u). Then W̃ = −u^{2B−2} W(1/u), and a
direct computation gives F_{p̃,q̃}(u) du = F_{p,q}(z) dx dy under u = 1/z (the inversion is a rotation of the
sphere), so the integral over |z| ≥ 1 is the integral of F_{p̃,q̃} over the unit disk. The half-disk
{|z| ≤ 1, Re z ≥ 0} is the image of the rectangle 𝒬 = [0,1] × [−1,1] under

    (ρ, t) ↦ z = ρ (1 + it)/(1 − it),        dx dy = 2ρ/(1 + t²) dρ dt,

and the other half-disk is its image under z ↦ −z. Hence

    π ℐ = H(p, q) + H(p(−·), q(−·)) + H(p̃, q̃) + H(p̃(−·), q̃(−·)),        H(p, q) = ∫_𝒬 G_{p,q}(ρ, t) dρ dt,

    G_{p,q} = 2ρ (1 + ρ²)² (1 + t²)³ |𝒲|⁴ / (|P|² + |Q|²)⁴,
    P = Σ_k a_k ρ^k (1 + it)^k (1 − it)^{B−k},   Q likewise from q,   𝒲 = Σ_k w_k ρ^k (1 + it)^k (1 − it)^{2B−2−k},

where a_k, w_k are the coefficients of p and W. (P = (1 − it)^B p(z), |1 − it|² = 1 + t², |z| = ρ.) G is a
quotient of polynomials in (ρ, t) whose denominator is positive on 𝒬.

**Lemma 3 (midpoint rule with remainder).** Let G be C² on a rectangle c with centre m and sides h₁, h₂. Then

    h₁h₂ [G(m) + (h₁²/24) inf_c ∂₁₁G + (h₂²/24) inf_c ∂₂₂G] ≤ ∫_c G ≤ h₁h₂ [G(m) + (h₁²/24) sup_c ∂₁₁G + (h₂²/24) sup_c ∂₂₂G].

*Proof.* For a C² function φ on an interval of length h with midpoint m, φ(s) − φ(m) − φ′(m)(s − m) lies between
½(s − m)² inf φ″ and ½(s − m)² sup φ″; integrating, ∫φ − hφ(m) lies between (h³/24) inf φ″ and (h³/24) sup φ″.
Apply this to ρ ↦ G(ρ, t) for each t and then to t ↦ G(m₁, t). ∎

**How ℐ is enclosed.** 𝒬 is divided into cells with dyadic corners. On each cell the value of G at the centre and
enclosures of ∂_ρρG and ∂_ttG over the cell are computed by evaluating the expression for G in interval arithmetic
on second-order Taylor data (value, two first and three second derivatives, each an interval): sums and products
by the Leibniz rule, and one reciprocal of the enclosure of |P|² + |Q|², which must be positive on the cell (if the
enclosure is not positive the cell is refined, and the computation fails if it cannot be made positive). The
coefficients of p and q are Gaussian rationals, enclosed by machine intervals. Lemma 3 then gives an upper and a
lower bound on each cell; the cell bounds are added exactly (as rationals) and divided by a rational lower or
upper bound of π. Only the upper bound enters a certificate.

**Degree.** That p and q have no common zero and that max(deg p, deg q) = B is decided exactly: the resultant of p
and q, a Gaussian rational, is not zero.

## 5. The radial integral

With p(r) = π − f(r) the integrand of Theorem 1 is

    e(r) = r² p′² + 2B sin² p (1 + p′²) + ℐ sin⁴ p / r² + 2μ² r² (1 + cos p).

The coefficient of ℐ is non-negative, so e increases if ℐ is replaced by a certified upper bound. The trial profile
is of the kind used for the hedgehog (Kernel/SkyrmeHedgehog.wl): p is piecewise cubic with rational data, p(0) = 0,
p rises linearly to π on the last cell, f = 0 beyond. The integral is enclosed cell by cell exactly as there, with
B and the upper bound of ℐ in the two places where the hedgehog has 1:

* first cell: exact integral of the polynomial majorant from |sin p| ≤ |p| and 1 + cos p ≤ 2;
* middle cells: Simpson and first-order enclosures in interval arithmetic;
* last cell: exact integral of the polynomial majorant from |sin f| ≤ |f|, 1 − cos f ≤ f²/2, 1/r² ≤ 1/r²_{n−1}.

For B = 1 and R = z (ℐ = 1) the competitor is the hedgehog and the computation is that of `SkyrmeDegreeOne`.

## 6. Choosing the map and the profile (heuristic)

The maps are those of [1] for B = 2, …, 8 (Section 9 lists them). Where a map of [1] contains √3 the competitor
uses a rational number close to √3 instead; where [1] gives a parameter to two or three digits the competitor uses
that decimal number. Any rational map of degree B is admissible, so this affects only the value of ℐ. The profile
is the minimiser of a finite-difference version of the radial energy, found by Newton's method in floating point,
and rounded to rational data. Neither choice is trusted.

## 7. Arithmetic and trusted base

Rigorous code: +, −, ×, ÷ of IEEE-754 doubles with the outward rounding of Kernel/RigorousArithmetic.wl (the
angular integral); exact rational arithmetic (the resultant, the sums over cells, the polynomial majorants);
Wolfram Language `Interval` arithmetic at 30 digits (the radial quadrature, as for the hedgehog). No library
elementary function, linear solver or optimiser enters a certified number. Trusted in addition: this note.

## 8. Checks

To be filled in after the implementation.

## 9. Results

To be filled in after the implementation.

## 10. What is not claimed

* **Not claimed:** that these maps are minimisers or even critical points of the energy. The minimal energies found
  numerically in the literature are lower (for example 1.171 against 1.208 per baryon for B = 2, published [1]).
* Not claimed: that I_B is attained, or that the certified number is the best that the ansatz gives (the profile
  and, for B = 5, 6, 8, the parameters of the map are only approximately optimal).
* Not claimed: any lower bound beyond B for B ≥ 2.

## References

1. C. J. Houghton, N. S. Manton, P. M. Sutcliffe, *Rational maps, monopoles and Skyrmions*, Nucl. Phys. B 510
   (1998), 507–537, arXiv:hep-th/9705151. (The ansatz (2.3), the energy (2.10), the integral ℐ (2.11), the maps of
   Section 3 and the energies of Table 1.)
2. R. A. Battye, P. M. Sutcliffe, *Skyrmions, fullerenes and rational maps*, Rev. Math. Phys. 14 (2002), 29–86,
   arXiv:hep-th/0103026. (Table 1: minimal ℐ and energies per baryon for B ≤ 22.)
3. S. Krusch, *Quantization of Skyrmions*, arXiv:hep-th/0610176. (Equation (5): the energy with the pion mass
   term, in the normalisation used here.)
4. L. D. Faddeev, *Some comments on the many-dimensional solitons*, Lett. Math. Phys. 1 (1976), 289–293. (The bound
   of Lemma 2; cited from memory.)

What was searched and found: to be filled in.

## Log

* 23 November 2026. Note written before the implementation. The energy formula of the task description was
  compared with [1] (2.10)–(2.11) and with [3] (1), (5), (6), and with the density in Kernel/SkyrmeHedgehog.wl at
  B = 1, ℐ = 1: they agree, including the factor 2μ² of the mass term.
* 23 November 2026, work stopped here (the machine was shut down). State for whoever resumes:
  **Written, never run in a kernel:** Kernel/SkyrmeRationalMaps.wl (interval arithmetic on second-order Taylor
  data, the adaptive midpoint rule for ℐ, the radial bound, the maps, the Newton solver for the profile, and the
  model functions `certifySkyrmeRationalMap`, `verifySkyrmeRationalMap`, `skyrmeRationalMapDescription`,
  `exportRationalMapCertificate`, `verifyRationalMapCertificateFile`). Expect bugs. The file is not yet loaded by
  Kernel/CertiField.wl and the model is not yet registered in Kernel/Models.wl.
  **Run (floating point only):** Reference/SkyrmeRationalMaps/check_rational_maps.py, function `angular`, on the
  maps of the kernel file: ℐ = 1.000000, 5.808259, 13.577235, 20.649626, 35.753385, 50.761040, 60.868042 for
  B = 1, …, 7 (literature: 1, π + 8/3, 13.58, 20.65, 35.75, 50.76, 60.87) and 85.63 near a = 0.136 for B = 8, each
  with degree B to twelve digits. **Error found in [1]:** for (3.17) the paper prints a = −3; that value gives
  ℐ = 216.53 (float), and a = −1/7 (which satisfies the same Wronskian condition) gives the icosahedral map with
  ℐ = 60.868; the kernel file uses a = −1/7. For B = 8 the paper's a = 0.14 gives 85.66; 0.136 gives 85.63.
  **Next, in order:** (1) load the kernel file in a script and run ℐ for B = 1, 2 against 1 and π + 8/3, then
  B = 3, …, 8; (2) one radial bound (B = 2, μ = 0, published 1.208 per unit of degree), and
  `rationalMapUpperBound[trial, 1, 1, μ]` against `hedgehogUpperBound` on the shipped trial data (must be equal);
  (3) add "SkyrmeRationalMaps.wl" to Kernel/CertiField.wl, the name "SkyrmeRationalMap" to `$CertiFieldModels`, the
  branches in `assembleCertificate` (raw and description) and in `verifyCertificate` of Kernel/Models.wl;
  (4) Tests/SkyrmeRationalMaps.wlt and a flag in tools/run_tests.wls: soundness of the Taylor arithmetic against
  exact derivatives at rational points, ℐ(z) = 1 and ℐ(z²) = π + 8/3 inside the enclosures, B = 1 against the
  hedgehog, published energies to 0.001, tampering, failure modes; (5) store certificates as
  Reference/SkyrmeRationalMaps/certificates/*.json and run check_rational_maps.py on them; (6) Sections 8–10 of
  this note, README, docs/ARCHITECTURE.md, a CI step for the Python check.
