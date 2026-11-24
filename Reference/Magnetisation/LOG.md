# Track M (magnetisation): log

Dated entries, newest last. Numbers marked [FLOAT] are ordinary floating point and are not certified.

## 23 November 2026

### State

* Done: analysis of the dual family for the thin-film model; a floating-point feasibility study
  (`explore_reduced_dual.py`); literature check of the items below.
* In progress: the note `docs/MAGNETISATION.md` (statements and proofs), then the reference implementation.
* Next: exact certificate for the lower bound, certified competitor energy and degree, tests M-T1 to M-T4, M1.

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
