# Nonlinear elasticity: what the method could certify (feasibility note)

**Status.** An assessment, not a model: no code, no certified numbers. Statements marked *classical* are cited;
the proposition of Section 2 has a short proof, checked by its author only, and is of a kind long known in the
uniqueness theory (Section 3); statements marked *sketch* are not proved here. No novelty is claimed.

**Setting.** Ω ⊂ ℝⁿ a box, n = 2 or 3; deformations u : Ω → ℝⁿ; E[u] = ∫_Ω W(∇u) − b·u dx, u = g on Γ_D.
W is *polyconvex*: W(F) = Φ(M(F)) with Φ convex and M(F) the list of minors, (F, det F) for n = 2 and
(F, cof F, det F) for n = 3. W(F) ≥ c|F|^p − C with p ≥ n, and 𝒜 = { u ∈ W^{1,p}(Ω; ℝⁿ) : u = g on Γ_D }. For
p ≥ n the minors of ∇u are integrable and are null Lagrangians (their integrals over Ω depend only on the values
of u on ∂Ω; classical [1, 2]). For p < n the determinant is not: x/|x| has det ∇u = 0 almost everywhere and the
boundary values of the identity. Nothing below applies there. The method is the principle of `docs/SCALAR_ND.md` extended to vector
fields: if ∫_Ω N(x, u, ∇u) dx = B(g) for every u ∈ 𝒜, then E[u] ≥ B(g) − ∫_Ω G, with
G(x) = sup over (y, F) of N(x, y, F) − W(F) + b·y.

## 1. Affine data: the certificate

Let Γ_D = ∂Ω, g(x) = Ax, b = 0. Take Λ in the subdifferential of Φ at M(A) and N(F) = Λ·M(F). Then
∫_Ω N(∇u) = |Ω| N(A) for every u ∈ 𝒜, and N(F) − W(F) ≤ N(A) − W(A) for every F. Hence E[u] ≥ |Ω| W(A): the
homogeneous deformation is a minimiser. This is the classical proof that polyconvexity implies quasiconvexity [1].

A certificate is a triple (A, Λ, δ): a rational matrix, rational multipliers (any Λ gives a valid bound; it need
not be an exact subgradient), and a certified number

    δ ≥ sup_F [ Λ·M(F) − W(F) ] − ( Λ·M(A) − W(A) ),      certifying   |Ω| (W(A) − δ) ≤ inf E ≤ |Ω| W(A).

The upper bound is the energy of u = Ax. The work is the supremum over F ∈ ℝ^{n×n}:

* If Φ is a sum of convex functions of F, of cof F and of det F (compressible neo-Hookean, Mooney–Rivlin and
  Ogden-type energies), the supremum is at most the sum of the three convex conjugates at Λ, with equality when Λ
  is the subgradient. Each is a closed form or a one-dimensional supremum of the kind the conjugate tables of
  `ScalarFieldND` already enclose.
* For a polynomial W without that structure, a sum-of-squares decomposition of
  W(F) − W(A) − Λ·(M(F) − M(A)) + δ, rounded to exact rationals, or interval branch and bound in n² variables.

This is cheap and would be a good unit test, but it re-proves a textbook statement. With traction-free faces it
is not available: the boundary term of a minor contains tangential derivatives of u on Γ_N, which are not data,
and the homogeneous state need not be the minimiser there (buckling).

## 2. Non-affine data: what the method can use, and where it stops

**(a) Divergences.** The energy has no growth in u. A field that is affine in u, S = P(x)ᵀu, contributes
N = P : ∇u + (div P)·u, and the supremum over y is finite only if div P + b = 0 exactly (rows of P), as in linear
elasticity. That is easy to arrange (rows of P are curls of continuous piecewise polynomial stream functions; P
need not be symmetric), and it gives the complementary-energy bound E[u] ≥ ∫_{Γ_D} (Pν)·g − ∫_Ω W*(P): the bound
for the convexified density W**. A field that is nonlinear in u is admissible only if its own divergence supplies
negative growth in u (this is how the same principle proves Poincaré inequalities [7]); whether that can help
here was not examined.

**(b) Minors.** With constant coefficients they are null Lagrangians. With coefficients that depend on x they are
not: for n = 2,

    r(x) det ∇u = div( r u₁ (∂₂u₂, −∂₁u₂) ) − u₁ ( ∂₁r ∂₂u₂ − ∂₂r ∂₁u₂ ),

and the last term is linear in u and in ∇u, with nothing in the energy to absorb it. So the multipliers of the
minors are constants c̄, and the fields available are N = div(Pᵀu) + c̄·m(∇u), m the minors of order at least
two, with G(x) = sup_F [ P(x) : F + c̄·m(F) − W(F) ] (the terms with b·u cancel).

**(c) How far that reaches.** Model class: n = 2, Γ_D = ∂Ω, W(F) = μ|F|²/2 + h(det F), μ > 0, h convex.

**Proposition.** Let u* ∈ 𝒜 be a weak solution of div DW(∇u*) + b = 0 with h′(det ∇u*) bounded, and let c̄ be a
constant. Then for every u ∈ 𝒜, with H = ∇u − ∇u*,

    E[u] − E[u*]  ≥  ∫_Ω ½ ( μ − |h′(det ∇u*) − c̄| ) |H|² dx.

In particular, if |h′(det ∇u*) − c̄| ≤ μ almost everywhere, u* is a global minimiser, and the principle with
P = DW(∇u*) − c̄ cof ∇u* and the multiplier c̄ returns E[u*]: there is no gap. Conversely (h twice
differentiable), if u* is a minimiser and some pair (P, c̄) returns E[u*], then |h′(det ∇u*) − c̄| ≤ μ almost
everywhere.

*Proof.* In two dimensions det(F₀ + H) = det F₀ + cof F₀ : H + det H. With F₀ = ∇u*, d₀ = det F₀, d = det ∇u,

    W(∇u) − W(F₀) − DW(F₀) : H = ½ μ|H|² + [ h(d) − h(d₀) − h′(d₀)(d − d₀) ] + h′(d₀) det H.

The bracket is non-negative. ∫ DW(F₀) : H = ∫ b·(u − u*) by the equation, and ∫ det H = 0 because u − u*
vanishes on ∂Ω, so c̄ ∫ det H may be subtracted. Finally |det H| ≤ ½|H|². For the converse, equality in the
principle forces F₀(x) to maximise F ↦ P : F + c̄ det F − W(F) for almost every x; the second-order condition
there, tested on the conformal and the anticonformal matrices H with cof F₀ : H = 0 (for which det H = ±½|H|²),
gives μ ∓ (h′(d₀) − c̄) ≥ 0. ∎

So the constant-multiplier bound is exact as long as the pressure-like quantity h′(det ∇u*) varies by at most 2μ
over the body, a condition on stress, not on the smallness of strain. Beyond it no field of the family attains
the minimum; presumably (*sketch*) the bound then stays a fixed distance below it however fine the grid. With
traction-free faces the minor is lost (c̄ = 0) and the window is |h′(det ∇u*)| ≤ μ; for a material that is
stress-free in the reference state, h′(1) = −μ, so the bound can be exact only if det ∇u* ≥ 1 everywhere: the
body is nowhere compressed in area. Compression with free faces is where buckling is possible, and the method
cannot be expected to decide it. In three dimensions
the remainder also contains F₀ : cof H and the cubic det H; the corresponding condition involves |∇u*| and needs
growth of order three. It was not worked out.

**Summary of what stops the method.** (1) u is unbounded and no lower-order term absorbs residuals: fields
affine in u must be exactly equilibrated (solvable). (2) Minors enter with constant coefficients only. (3) Hence
the bound is exact only inside the window above; outside it the certified interval remains valid and does not
close. (4) The supremum over F is a global optimisation in n² variables at every vertex (it reduces to one
variable for the model class, through the singular values of P; *sketch*).

## 3. What the literature already does

* Sum of squares. The pointwise dual relaxation of Chernyavsky, Bramburger, Fantuzzi and Goluskin [7] is the
  principle above with divergences of polynomial fields φ(x, u) only, the pointwise inequality being imposed by
  sum-of-squares conditions and solved as a semidefinite programme in floating point. It is sharp for several
  classes (eigenvalues of Sturm–Liouville problems, Poincaré constants). Fantuzzi and Tobasco [8] show that it
  coincides with the occupation-measure relaxation under a coercivity hypothesis, that both are sharp for
  one-dimensional, scalar-valued or convex problems, and that they are not sharp in general when domain and
  target both have dimension at least two (a double-well example with incompatible wells); see also [9]. In
  non-convex cases the relaxation is at least as strong as the convex envelope [10]. This matches (a): without
  minors one gets the convexified problem.
* Polyconvexity. Fantuzzi, Henrion, Kružík, Murali and Weis [11] give sum-of-squares conditions that certify
  polyconvexity of a polynomial W and compute polyconvex envelopes pointwise. That is the pointwise inequality
  of Section 1 for general polynomial W; they do not bound minima of boundary value problems.
* Global minimality by analysis. Uniqueness and global minimality of equilibrium solutions under small-strain or
  small-stress conditions: John [3], Knops and Stuart [4] (affine data), Zhang [5], Sivaloganathan and Spector
  [6] (uniform polyconvexity and pointwise bounds on derivatives of the stored energy at the solution). The
  proposition of Section 2 is of this kind.

Searched (November 2026): the items above, and "certified" or "rigorous" two-sided bounds for boundary value
problems of nonlinear elasticity with non-affine data. None was found. This is not a literature review.

## 4. Recommendation

**No-go** for nonlinear elasticity in general: outside the window of Section 2 (large pressure variation,
compression with free faces, instabilities) the families of fields available cannot close the gap, and that is a
property of constant minor multipliers, not of the grid. Three dimensions need the analysis that was not done.

**Go, as a small separate model, for a restricted class**: n = 2, W(F) = μ|F|²/2 + h(det F) with polynomial h
and rational data, displacement prescribed on all of ∂Ω. Inside the window the bound is exact in the limit, and
every ingredient exists: exactly equilibrated P from stream functions; G convex in P, so the cell integral is
bounded by vertex values as in Proposition 2 of `docs/SCALAR_ND.md`; a polynomial competitor energy, evaluated
exactly; and the certified gap itself reports when the window has been left.

**The toy problem that decides it.** Ω = (0, 1)², W(F) = ½|F|² + ½(det F − 1)² − det F (polyconvex, not convex,
stress-free at the identity, μ = 1, h′(d) = d − 2), datum g(x) = (x₁ + t x₁x₂, x₂) on ∂Ω, for which det ∇g
ranges over [1, 1 + t]. Run t = 1/2, 2 and 4: by the proposition the gap should close under refinement for small
t and stall once the variation of det ∇u* exceeds 2. Go if the certified gap at t = 1/2 falls below 1 % in
minutes and the stall appears where predicted; no-go if the gap does not close even at t = 1/2. The cost is
that of the two-dimensional elasticity model.

## Not claimed

Nothing here is implemented or certified. The proposition is stated for the two-dimensional model class only.
Beyond the window only non-attainment is proved, not a gap of positive size. The reduction of the supremum over
F and everything said about three dimensions are sketches.

## References

1. J. M. Ball, *Convexity conditions and existence theorems in nonlinear elasticity*, Arch. Rational Mech. Anal.
   63 (1977), 337–403.
2. J. M. Ball, J. C. Currie, P. J. Olver, *Null Lagrangians, weak continuity, and variational problems of
   arbitrary order*, J. Funct. Anal. 41 (1981), 135–174.
3. F. John, *Uniqueness of non-linear elastic equilibrium for prescribed boundary displacements and sufficiently
   small strains*, Comm. Pure Appl. Math. 25 (1972), 617–634.
4. R. J. Knops, C. A. Stuart, *Quasiconvexity and uniqueness of equilibrium solutions in nonlinear elasticity*,
   Arch. Rational Mech. Anal. 86 (1984), 233–249.
5. K. Zhang, *Energy minimizers in nonlinear elastostatics and the implicit function theorem*, Arch. Rational
   Mech. Anal. 114 (1991), 95–117.
6. J. Sivaloganathan, S. J. Spector, *On the uniqueness of energy minimizers in finite elasticity*, J. Elasticity
   133 (2018), 73–103; arXiv:1804.06971.
7. A. Chernyavsky, J. J. Bramburger, G. Fantuzzi, D. Goluskin, *Convex relaxations of integral variational
   problems: pointwise dual relaxation and sum-of-squares optimization*, SIAM J. Optim. 33 (2023), 481–512;
   arXiv:2110.03079.
8. G. Fantuzzi, I. Tobasco, *Sharpness and non-sharpness of occupation measure bounds for integral variational
   problems*, arXiv:2207.13570.
9. M. Korda, R. Rios-Zertuche, *The gap between a variational problem and its occupation measure relaxation*,
   arXiv:2205.14132.
10. D. Henrion, M. Korda, M. Kružík, R. Rios-Zertuche, *Occupation measure relaxations in variational problems:
    the role of convexity*, arXiv:2303.02434.
11. G. Fantuzzi, D. Henrion, M. Kružík, A. Murali, S. Weis, *Polyconvexity with moments and sums of squares*,
    arXiv:2604.11124.

The bibliographic data of items 1–11 were checked by web search in November 2026. Items 7–11 were read in
abstract (item 7 also in its introduction); items 3–6 are cited for what their abstracts and later citations
say they prove, not from a reading of the full texts.
