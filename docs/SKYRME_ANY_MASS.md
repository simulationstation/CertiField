# Skyrme model, degree one: certified bounds for any pion mass

This note states what the model `SkyrmeDegreeOne` certifies when the pion mass is not one of the four masses for
which grids are shipped, and why. It adds no new mathematics to the level-set lower bound: it records that the
bound and its certified evaluation do not depend on the four masses, proves the one-line monotonicity statement
that is used as a fallback, and describes the heuristic that places the grids for a new mass.

**Status.** The lower bound is the level-set theorem of [1, 2], which this note uses and does not re-prove; the
statements proved here (Section 3) are elementary and have been checked by their author only. The grid generator
(Section 5) and the trial profile (Section 6) are heuristics: they decide how tight a certificate is, never whether
it is valid. No novelty is claimed. Numbers quoted as *certified* come from the paclet's outward-rounded and exact
rational arithmetic; numbers marked (float) are ordinary floating point.

## 1. The problem

Units as in the paclet: for U : ℝ³ → S³ and a pion mass μ ≥ 0,

    E_μ(U) = ∫_{ℝ³} |dU|² + |Λ²dU|² + 2μ² (1 − cos f) dx,        f = dist_{S³}(U, 1),

reported as E_μ/12π². 𝒞_μ is the class of weakly differentiable maps U with E_μ(U) < ∞, U − 1 ∈ L⁶(ℝ³) and degree
one (the class of [1, 2], in which the degree is defined by the integral of the pulled-back volume form), and

    I₁(μ) = inf { E_μ(U) : U ∈ 𝒞_μ }.

## 2. The level-set bound for a general mass

**Theorem A ([1, 2]; used here, not re-proved).** Let μ ≥ 0 be rational, let 0 < τ₁ < … < τ_N < π be machine numbers
and let any volume nodes be given at each level. The programme `levelSetLowerBound` (Kernel/SkyrmeLevelSet.wl),
run with the exact rational μ and these grids, returns a machine number c with

    E_μ(U) ≥ 12π² c        for every U ∈ 𝒞_μ.

The analytic part (coarea formula on the level sets of the distance to the antivacuum, the flux identity, the
isoperimetric inequality, affine minorants of h_κ, the step lemma, the three closed-form floors, monotonicity of the
value function in the volume) is written for the potential V(τ) = 2μ²(1 + cos τ) with an arbitrary μ ≥ 0. The four
masses 0, 263/500, 13/20, 1 enter the published results only as the values for which the programme was run.

The mass enters the programme in four places, each of which is enclosed for an arbitrary rational μ:

| Quantity | Where | Enclosure |
|---|---|---|
| V(τ_i) = 2μ²(1 + cos τ_i) | the step bound, the far targets | exact rational μ², rational enclosure of cos τ_i (alternating Taylor series with remainder), rounded outward to machine numbers |
| 2μ² min(sin τ_i, sin τ_{i+1}) | condition (∗) of the step lemma | exact rational μ², rational lower bound of the sine, rounded down |
| μ in Ξ_μ = 16π cos²(τ/2)(r + μr²) | the floor Ξ_μ and the terminal data | the largest machine number ≤ μ |
| the grids | everywhere | heuristic; any grid gives a valid bound |

Nothing else depends on the mass. In particular the table of affine minorants of h_κ is the same for every mass;
where a multiplier falls outside the table (κ > 400) the programme uses the last certified row, which is a valid
minorant with a smaller κ.

## 3. Monotonicity in the mass

**Lemma B.** Let 0 ≤ μ₀ ≤ μ. Then 𝒞_μ ⊂ 𝒞_{μ₀}, E_μ(U) ≥ E_{μ₀}(U) for every U ∈ 𝒞_μ, and I₁(μ) ≥ I₁(μ₀). Hence a
lower bound for E_{μ₀} on 𝒞_{μ₀} is a lower bound for E_μ on 𝒞_μ.

*Proof.* 1 − cos f ≥ 0, so E_{μ₀}(U) ≤ E_μ(U) for every map U. If U ∈ 𝒞_μ this is finite, and the other two
conditions in the definition of the class do not involve the mass; so U ∈ 𝒞_{μ₀}, and
E_μ(U) ≥ E_{μ₀}(U) ≥ I₁(μ₀). ∎

**Lemma C (the trial map).** The trial maps of the paclet (Kernel/SkyrmeHedgehog.wl: U = exp(i f(r) x·σ/r) with
f = π − p, p piecewise cubic with rational data, p(0) = 0, p rising linearly to π on the last cell, f = 0 beyond it)
belong to 𝒞_μ for every μ ≥ 0, and μ ↦ E_μ(U) is non-decreasing.

*Proof.* The map is Lipschitz, equals 1 outside a ball and has degree one; its energy is finite for every μ because
the integrand vanishes outside the ball. Monotonicity is as in Lemma B. ∎

Lemma B is used as a fallback: a certified lower bound at a mass μ₀ holds at every mass μ ≥ μ₀, at no cost. It is
weaker than the bound computed at μ itself by about the amount that I₁ grows between μ₀ and μ.

## 4. What a certificate asserts

Let μ be the requested mass (a rational number; a machine number is replaced by a rational within 10⁻¹² of it, as
before, and the certificate names the rational that was used).

**Theorem D.** A `CertifiedBound` for `SkyrmeDegreeOne` with pion mass μ contains a rational μ_L with
0 ≤ μ_L ≤ μ, a grid, and rational trial data. Let c be the number returned by `levelSetLowerBound` at the mass μ_L
on that grid, and let u be the upper end of the enclosure of E_μ(H)/12π² computed by `hedgehogUpperBound` for the
trial map H. Then

    12π² c ≤ I₁(μ) ≤ 12π² u.

*Proof.* Theorem A at the mass μ_L and Lemma B give the left inequality. H ∈ 𝒞_μ by Lemma C, which gives the
right one. ∎

Which μ_L and which grid:

| Requested μ | `Method` | μ_L | Grid |
|---|---|---|---|
| 0, 263/500, 13/20 or 1 | `Automatic` | μ | the shipped grid of that mass and resolution (unchanged behaviour) |
| any other rational in the supported range | `Automatic` | μ | generated for μ at the requested resolution (Section 5) and stored in the certificate |
| any rational in the supported range | `"Monotonicity"` | the largest of 0, 263/500, 13/20, 1 that is ≤ μ | the shipped grid of μ_L |
| any rational in the supported range | `"GeneratedGrid"` | μ | generated for μ even if a grid is shipped (used to test the generator) |

If the generator fails for a requested mass, `Automatic` falls back to `"Monotonicity"` and the certificate says so.
The upper bound always refers to the requested mass μ itself.

`VerifyCertificate` recomputes c from (μ_L, the stored or shipped grid) and u from (μ, the stored trial data), checks
0 ≤ μ_L ≤ μ, and compares with the stored bounds. It does not run the generator or the profile solver again.

## 5. The grids (heuristic)

The grids are those of the original verifier named in the header of Kernel/SkyrmeLevelSet.wl; its construction is
reproduced in Kernel/SkyrmeAnyMass.wl. Nothing in this section is needed for validity.

1. *Reference path.* The relaxed level-set problem has an extremal described by the radius r(τ) of the ball with
   the volume of {dist(U, −1) < τ} and a multiplier κ(τ):

       r′ = ∂_κ h**_κ(D),     κ′ = −2μ² sin τ − (2/r) [h**_κ − D ∂_D h**_κ](D),     D = sin²τ / r²,
       r ≈ τ/a and κ = a⁴ − a² at the centre,

   with h_κ(D) = 2√((1 + 2D)(D² + 2D + κ)) and h**_κ its convex envelope (h_κ for D ≥ D_t(κ), the tangent line
   through (0, 2√κ) below; D_t = 0 when ½(3 − √5) ≤ κ ≤ ½(3 + √5)). Where the bracket is negative (D above the root
   D* of D*² − D* = κ) the area constraint is inactive and r′ = sin²τ ∂_κh_κ(D*)/(D* r²), κ′ = −2μ² sin τ. The slope a
   is found by bisection between paths whose multiplier κ − V reaches zero before τ = π and paths for which it does
   not (for μ = 0: by the sign of r√κ − (π − τ) at the end).
2. *Reference radius.* r_ref(τ) is the monotone cubic interpolant of ln r in the variable ln(τ/(π − τ)), continued
   linearly; v_ref = (4π/3) r_ref³ and λ_ref = d ln v_ref/dτ.
3. *Levels.* τ₁ = 0.03, τ_N = π − 0.02, steps min(0.004, δ/λ_ref(τ)).
4. *Volume nodes.* At level i: 0 and v_ref(τ_i) e^u for a fixed set of u with spacing δ on −4 ≤ u ≤ 2, growing by
   5 % per cell outside (at most 0.2) from −32 to 22; nodes below the radius 0.1 sin τ_i are dropped.

δ is 0.04, 0.02 and 0.005 for the resolutions `"Quick"`, `"Standard"` and `"Publication"`.

Differences from the original construction, all without effect on validity: the differential equations are
integrated by an embedded Runge–Kutta pair of orders 5 and 4 with the original tolerances (the original uses
LSODA); D_t(κ) is computed from the equation (D² − 1 − κ)² = κ(4D + 5), which is the tangency condition with the
square root removed, by bisection of the increasing function D² − 1 − κ + √(κ(4D + 5)).

## 6. The trial map for a general mass (heuristic)

With h_κ in place of h**_κ the system of Section 5 is the hedgehog equation in Hamiltonian form, with τ = π − f.
The same shooting code therefore gives the hedgehog profile as r(τ). It is inverted by monotone interpolation,
continued beyond the reliable radius by the linearised tail (f ∝ (1 + μr) e^{−μr}/r², or 1/r² for μ = 0), sampled at
the nodes of `hedgehogTrialData` and rounded to rationals. The energy of the resulting trial map is enclosed by
`hedgehogUpperBound` exactly as for the shipped trial maps. Any profile data give a valid upper bound.

## 7. Arithmetic and trusted base

Unchanged: the certified numbers are produced by `levelSetLowerBound`, the certified minorant table and
`hedgehogUpperBound`, which are not modified. Trusted: IEEE-754 arithmetic (+, −, ×, ÷) with the outward rounding of
Kernel/RigorousArithmetic.wl, exact rational arithmetic, Wolfram Language `Interval` arithmetic in the quadrature,
the analytic reduction of [1, 2], and Section 3 of this note. Not trusted and not used in a certified number:
the differential-equation solver, the interpolation, the bisection, the choice of grids.

## 8. Checks

To be filled in after the implementation.

## 9. Results

To be filled in after the implementation.

## 10. Supported range, and what is not claimed

* Supported input: to be fixed after the tests of Section 8.
* **Not claimed:** that the hedgehog is the minimiser for any mass. The gap bounds how far below the hedgehog trial
  map a degree-one map can be.
* Not claimed: that I₁(μ) is attained for μ > 1. The bounds are bounds for the infimum over 𝒞_μ whether or not it is
  attained.
* Not claimed: that the generated grids are the best ones. A different grid may give a larger certified number.
* Irrational masses are not accepted as input. By Lemma B and Lemma C a certificate at a rational μ₀ ≤ μ gives a
  lower bound, and one at a rational μ₁ ≥ μ an upper bound, for an irrational μ.

## References

1. The preprint and the manuscript cited as [1] and [2] in the README of this repository (the level-set lower
   bound and its computer-assisted evaluation; the header of Kernel/SkyrmeLevelSet.wl names the note and the
   original verifier).
2. M. J. Esteban, *A direct variational approach to Skyrme's model for meson fields*, Comm. Math. Phys. 105 (1986),
   571–591; *Existence of 3D Skyrmions. Complete version*, arXiv:math-ph/0401042 (2004). M. J. Esteban, S. Müller,
   *Sobolev maps with integer degree and applications to Skyrme's problem*, Proc. R. Soc. Lond. A 436 (1992),
   197–201. (The class of maps.)
3. D. Harland, *Topological energy bounds for the Skyrme and Faddeev models with massive pions*, Phys. Lett. B 728
   (2014), 518–523, arXiv:1311.2403. (A lower bound that depends on the pion mass; it is a topological bound and is
   further from the hedgehog energy than the level-set bound.)

What was searched and found: to be filled in.

## Log

* 23 November 2026. Note written before the implementation. The frozen files were read; the mass enters
  `levelSetLowerBound` only through the four quantities of Section 2.
