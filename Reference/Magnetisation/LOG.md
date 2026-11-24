# Track M (magnetisation): log

Dated entries, newest last. Numbers marked [FLOAT] are ordinary floating point and are not certified.

## 23 November 2026

### State (updated at the end of the day's entries below)

* Done: the note `docs/MAGNETISATION.md`, Sections 1 to 6 (statements and proofs); the exact checker for the
  lower-bound certificate (`magnet.py`), with known answers, falsification and an independent exact test of the
  pointwise inequality (`check_dual.py`, `pointwise.py`: all pass); the heuristic search for certificates
  (`dual_opt.py`); the first certified lower bound for the skyrmion class (`results/zeeman_h2.json`).
* In progress: certificates for other fields (`run_dual.py`); the certified competitor (energy and degree of a
  normalised bilinear interpolant, Section 5 of the note) is designed and not yet written.
* Next: `upper_bound` and `degree` with interval arithmetic; M-T2 and the critical family from above; M-T3 (strip,
  cross-check with `MagneticDomainWall`); M-T1 (enclosure of 48 J); M1 intervals at a large and a moderate field;
  stored certificates and `verify_all.py`; Sections 7 and 8 of the note; then the paclet model.
* Open: the search starts from phi = kappa and repairs infeasible cells by damping towards phi = kappa, which is
  strictly feasible only for h > kappa^2/(2A); for lower fields (the dynamic programme gives bounds down to
  h = 0.9) another starting point and repair are needed. Degree Q > 0 uses the constant certificate only.

### The dual family that is being built

Model: m : Ω → S², Ω ⊂ ℝ² bounded, m = e₃ on ∂Ω, degree Q,
E[m] = ∫ A|∇m|² + κ m·(∇×m) + K(1 − m₃²) + h(1 − m₃).

1. For maps that are constant on ∂Ω the families (T) and (Q) of the prompt with x-independent coefficients are one
   statement: for every bounded Borel function c on S², ∫_Ω c(m) m·(∂₁m × ∂₂m) dx = Q ∫_{S²} c dA. (Area formula;
   the degree is the same integer at almost every point of the sphere.)
2. Supremum over the tangent variables p = ∂₁m, q = ∂₂m, with J = y × · and ℓ₁, ℓ₂ the tangential linear
   coefficients: sup = |ℓ₁ + Jℓ₂|²/(4(2A + c)) + |ℓ₁ − Jℓ₂|²/(4(2A − c)), finite for |c| < 2A.
3. For fields that are independent of x and equivariant under rotation about e₃, S(y) = φ(y₃)(y₂, −y₁), c = c(y₃),
   with t = y₃ and u = 2A + c:

       |ℓ₁ + Jℓ₂|² = (1 − t)² [((1 + t)(φ + κ))′]²,     |ℓ₁ − Jℓ₂|² = (1 + t)² [((1 − t)(φ − κ))′]²,

   and the pointwise inequality holds if and only if the polynomial inequality
   Π(t) = P(4A − u) + R u − V u (4A − u) ≤ 0 holds on [−1, 1] with 0 ≤ u ≤ 4A, where P, R are the two squares
   divided by 4 and V(t) = K(1 − t²) + h(1 − t). Then E[m] ≥ 2πQ(∫₋₁¹ u dt − 4A) on every domain. For piecewise
   polynomial φ, u with rational coefficients this is checked in exact rational arithmetic and the bound is an
   exact rational multiple of π.
4. Checks of the formulas against known results: φ ≡ κ, u ≡ 4A gives the Belavin–Polyakov bound in the class
   where the Dzyaloshinskii–Moriya term does not help; for V = (1 − t)²/2, A = 1/2, κ = r ≤ 1 the choice φ ≡ r,
   u ≡ 2r² gives E ≥ 4π(1 − 2r²) in degree −1, which is the exact infimum proved by Döring and Melcher (2017) and
   quoted in arXiv:2512.24598 (Theorem 1.1 there gives all degrees). φ ≡ κ, u = κ²(1 − t)/h gives
   E ≥ 8πA − 4πκ²/h for the Zeeman model.

### Feasibility study [FLOAT]

`explore_reduced_dual.py`: A = 1/2, κ = 1, K = 0 (the model of Li and Melcher, arXiv:1711.07717). Column 2: energy
of the axisymmetric skyrmion on the plane (Newton on the discrete radial energy). Column 4: the value of the
reduced dual of item 3, optimised by dynamic programming over (ln(1 − t), (1 + t)(φ + κ)) with 260 × 1500 grid
points. Energies in units of 4π.

| h | skyrmion energy | radius | simple bound 1 − 1/h | optimised dual | gap, simple | gap, optimised |
|---|---|---|---|---|---|---|
| 0.90 | 0.130321 | 1.2964 | −0.11111 | 0.11568 | 185 % | 11.2 % |
| 1.00 | 0.234262 | 1.1431 | 0.00000 | 0.21994 | 100 % | 6.1 % |
| 1.25 | 0.416274 | 0.8663 | 0.20000 | 0.39965 | 52 % | 4.0 % |
| 1.50 | 0.533166 | 0.6922 | 0.33333 | 0.52573 | 37 % | 1.4 % |
| 2.00 | 0.672755 | 0.4800 | 0.50000 | 0.66859 | 26 % | 0.62 % |
| 3.00 | 0.802490 | 0.2835 | 0.66667 | 0.79226 | 17 % | 1.27 % |
| 5.00 | 0.895823 | 0.1451 | 0.80000 | 0.89002 | 11 % | 0.65 % |
| 10.0 | 0.956203 | 0.0584 | 0.90000 | 0.95416 | 5.9 % | 0.21 % |
| 20.0 | 0.981403 | 0.0238 | 0.95000 | 0.98024 | 3.2 % | 0.12 % |
| 50.0 | 0.993868 | 0.0075 | 0.98000 | 0.99355 | 1.4 % | 0.03 % |

Reading: the optimised translation-invariant dual is within about 1 % of the skyrmion energy for h ≥ 1.5 and
within 12 % down to h = 0.9, where the simple bound is useless. The values are not monotone in h at the level of
0.5 % (grid effects of the dynamic programme); the optimum is at least as good as each value. These are floating
point values of an upper bound for the optimal dual functional, not certificates; the certificate is the exact
check of item 3 for a spline close to the optimiser's output.

An exact calibration of this kind does not exist for the skyrmion: equality in the pointwise inequality along the
skyrmion would need c = −(A|∇m|² − V)/|m·(∂₁m × ∂₂m)|, which is below −2A in the tail (there
A(θ′ + sin θ/r)² ≈ hθ² > h(1 − cos θ) ≈ hθ²/2). The dual is nevertheless close.

On the plane the dual functional is concave and invariant under translations and under simultaneous rotations, so
fields independent of x and equivariant are optimal there. On a bounded square with m = e₃ on the boundary fields
that depend on x can only gain through the boundary; this has not been quantified.

### Literature checked today

* X. Li, C. Melcher, arXiv:1711.07717: the model E = ∫ ½|∇m|² + m·(∇×m) + (h/2)|m − e₃|², degree −1, the
  axisymmetric form and its ODE; "a major open question is whether m₀ minimizes the energy even globally within
  the homotopy class Q = −1".
* S. Ibrahim, T. Miura, C. Román, I. Shimizu, arXiv:2512.24598: E_r = D + rH + V with V = ½∫(1 − n₃)²; for
  0 < r ≤ 1 the infimum in degree k is 4π|k|(1 − 2r²) (k < 0) and 4πk (k ≥ 0); the degree −1 case is attributed
  to Döring and Melcher, with the explicit minimiser h(x/2r), h the Belavin–Polyakov map.
* S. Komineas, C. Melcher, S. Venakides, Nonlinearity 33 (2020) 3395: small-radius asymptotics; the lower bound
  for the energy is stated there as resting on a heuristic.

### Exact certificates (evening)

* `check_dual.py` (20 s): constant certificates reproduce 8πA|Q| (exchange only), 4π(1 − 2r²) for the critical
  family at r = 1/4, 1/2, 3/4, 1 (and u lowered by 10⁻⁶ is rejected in each case), 4π(1 − 1/h) for the Zeeman
  model, and 4√(AK)w for the wall. Of 400 random certificates for the critical family 60 are accepted and none
  claims more than the exact infimum (smallest margin 0, attained by the exact one); of 300 random certificates
  for exchange only 107 are accepted, none above 8πA|Q|. For 200 random fields with a positive Γ the pointwise
  inequality f − N + Γ ≥ 0 holds exactly at random rational points (evaluated from the three-dimensional
  definitions, without the reduced formulas).
* Search, first attempts and what failed: (1) dynamic programme path smoothed by a spline: the value of the
  programme depends strongly on the grid in a (0.630 with 200 × 900 points against 0.669 with 260 × 1500 at h = 2)
  and 10 of 27 cells had no valid u; abandoned. (2) L-BFGS on the nodal values and slopes of a C¹ cubic φ on a
  mesh geometric in 1 − t, objective ½∫u₋: at h = 2 the value is 0.6710 with the mesh ending at 1 − t = 2·10⁻³,
  but the field still varies below that scale and the last cell failed. (3) The same with φ ≡ κ imposed on the
  last cell and the mesh extended in stages (2·10⁻³, 2·10⁻⁴, 2·10⁻⁵), each stage started from the previous field:
  one cell without a root of Π₀ between the quadrature points; repaired by damping the field towards φ ≡ κ by
  0.3 %.
* Result, h = 2 (A = 1/2, κ = 1, K = 0), degree −1: the certificate with 41 cells passes the exact check (0.2 s);

      E[m] ≥ 8.448613244 = 4π · 0.672319      for every m ∈ 𝒜₋₁, on every bounded Lipschitz domain.

  The energy of the axisymmetric skyrmion on the plane is 4π · 0.672755 [FLOAT], so the lower bound is within
  0.07 % of it. The simple bound is 4π · 0.5. The certified upper bound on a square is still to be computed.

### Stop of the session, 23 November 2026 (work paused here)

**Done and committed.**

* Note `docs/MAGNETISATION.md`: Sections 1 to 6 with complete proofs (Lemma 1 degree and weighted degree formula,
  Lemma 2 divergences, Lemma 3 supremum over the gradient, Theorem 4, Proposition 5 with the three checks against
  known results, Lemma 6 and Proposition 7 for the strip, Lemma 8 for the cube, Section 5 formulas for the
  competitor). Sections 7 (checks) and 8 (results) are placeholders.
* `magnet.py`: exact checker (`Certificate.check`, Bernstein coefficients with bisection in rational arithmetic),
  exact bound, `wall_bound`. `check_dual.py` passes (output in the entry above). `pointwise.py`: the independent
  exact test.
* `dual_opt.py`, `run_dual.py`: heuristic search and stored certificates. Certified lower bounds, Zeeman model
  (A = 1/2, κ = 1, K = 0), degree −1, valid on every bounded Lipschitz domain; each passed the exact check (0.2 s)
  and the independent pointwise test:

  | h | certified: E ≥ | in units of 4π | skyrmion on the plane [FLOAT] | simple bound | file |
  |---|---|---|---|---|---|
  | 3/2 | 6.697023218 | 0.532932 | 0.533166 | 0.333333 | `results/zeeman_h3o2.json` |
  | 2 | 8.448613244 | 0.672319 | 0.672755 | 0.500000 | `results/zeeman_h2.json` |
  | 3 | 10.079480953 | 0.802100 | 0.802490 | 0.666667 | `results/zeeman_h3.json` |

  The lower bounds are within 0.04 % to 0.07 % of the floating-point skyrmion energy. The runs for h = 5 and
  h = 10 were stopped before they finished.

**Half-written (exact state).**

* `competitor.py` (`upper_bound`, `degree`, `graded_axis`, `radial_field`, `belavin_polyakov`) is written from
  Section 5 of the note and HAS NOT BEEN RUN: its first test (cut-off Belavin–Polyakov maps for exchange only and
  for the critical family, which must give energies just above 4π and just above 2π with degree −1 or +1) was
  stopped before it printed anything. Treat every line as untested. Things to check first: the Bernstein
  coefficients of |M_h|² (`_D_coefficients`), the tables `G2`, `T2`, `G12`, the sign of the
  Dzyaloshinskii–Moriya term and the helicity convention of `radial_field` (the skyrmion of the note is
  (sin θ e^{i(ψ+π/2)}, cos θ) with θ(0) = π; its degree is −1), the orientation in `degree`.
* No test has been written for the competitor, no certified upper bound exists yet, hence no certified interval
  and no gap for M1 yet: the percentages above compare with a floating-point number.

**Tests of the prompt.**

* M-T1 (cube): proved in the note (Lemma 8); the enclosure of 48AJ by rigorous quadrature is not written.
* M-T2 (exchange only): the lower bound 8πA|Q| is reproduced exactly (`check_dual.py`); the competitors from
  above wait for `competitor.py`.
* M-T3 (strip): proved in the note (Proposition 7); `wall_bound` gives 4√(AK)w exactly; the competitor and the
  comparison with `MagneticDomainWall` are not done.
* M-T4 (falsification): done for the lower bound (random certificates, lowered u, random Γ, independent
  pointwise test); the pairing with competitors is not done.
* Extra known answer found: the critical family with exact infimum 4π(1 − 2r²), reproduced exactly.

**Next steps, in order.**

1. Run and debug `competitor.py` on the two known cases (exchange only; critical family r = 1/2 with a
   Belavin–Polyakov map of radius 1); add an independent high-accuracy quadrature of the same normalised
   interpolant as a cross-check; then `check_competitor.py`.
2. M1: skyrmion profile by Newton on the radial energy with θ(L) = 0 (code in `explore_reduced_dual.py`,
   `skyrmion`), nodal values on a graded grid, certified upper bound and degree on a square of half-width about
   6/√(2h); certified intervals at h = 2 or 3 (moderate) and h = 10 or 20 (large). Gate G2 is then a matter of
   the competitor's accuracy; the lower bounds are already within 0.1 %.
3. M-T3 competitor and cross-check with the one-dimensional model (needs the Wolfram kernel once); M-T1
   enclosure of J; falsification with competitors; `verify_all.py`; Sections 7 and 8 of the note.
4. Make the search robust: lower fields (h ≤ κ²/(2A), where φ ≡ κ is not feasible), K ≠ 0, degree Q > 0, and a
   faster optimiser (analytic gradients; the numerical-gradient L-BFGS takes about a minute).
5. Then the paclet model `ThinFilmMagnet` (exact check with rational arithmetic in the Wolfram Language,
   competitor with the paclet's interval helpers).
