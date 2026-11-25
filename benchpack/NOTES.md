# Benchmark pack: method, derivations, tests, limits

Labels as in the ground rules (PROVED / CERTIFIED / FLOAT / DIAGNOSTIC / OPEN). Everything here is checked by the
agent only. The certificate data are described in `certificates/FORMAT.md`; the problem data are in
`problems/definitions.py`; the builder code is in `builder/`.

## 1. Bounds on the energy (both problems)

Linear elasticity, plane stress, E = 1, ν = 3/10, unit thickness, δ = 1; `C` the stiffness, `S = C^-1` the
compliance (definitions.py). `V` = displacements in H¹ that satisfy the prescribed displacement conditions,
`V0` = the same with zero data. `a(u,v) = ∫ ε(u):C:ε(v)`. The exact solution `u ∈ V` minimises `(1/2) a(v,v)` over
`V`, so `a(u, w) = 0` for `w ∈ V0`. `K` is `a(u,u)` (P2, full domain) or `4 a_q(u,u)` (P1, quarter domain).

**Statically admissible stress.** A symmetric, cellwise polynomial `τ` with `div τ = 0` in every cell, continuous
traction across every cell edge, zero traction on the free boundaries, and zero shear traction on mirror lines and
frictionless grips. Integrating by parts cell by cell, for every `v ∈ V`

    ∫ τ:ε(v) = W(τ) := work of the traction of τ on the prescribed boundary displacement,

because the traction terms on interior edges cancel, vanish on free boundaries, and on mirror lines and grips only
the prescribed normal displacement does work. P1 quarter: `W(τ) = (1/2) F_q(τ)`, `F_q = ∫_0^1 τ_xx(1,y) dy`.
P2: `W(τ) = F(τ) = −∫_0^1 τ_xy(1,y) dy`.

**Upper bound (PROVED, classical).** `(1/2) a(u,u) <= (1/2) a(v,v)` for every `v ∈ V`.

**Lower bound (PROVED, classical).** For symmetric `σ`, `e`: `(1/2) e:C:e >= σ:e − (1/2) σ:S:σ`, because the
difference is `(1/2)(e − Sσ):C:(e − Sσ) >= 0`. With `σ = τ` statically admissible and `e = ε(v)`, integrating gives
`(1/2) a(v,v) >= W(τ) − (1/2) Q(τ)`, `Q(τ) = ∫ τ:S:τ`, for every `v ∈ V`. Scaling `τ` gives `K >= W'^2/Q` with the
force normalisation of each problem: P1 `K >= F_q(τ)^2 / Q_q(τ)`, P2 `K >= F(τ)^2 / Q(τ)`.

## 2. Bounds on the patch average (Q2): polarisation with a prestress

`l(v) = ∫ s:ε(v)` with `s = w C (e_x⊗e_x)`, `w = 1/|P|` on the patch `P` and 0 elsewhere, so `l(v)` is the average of
`σ_xx(v) = (Cε(v))_xx` over `P`, and `Q2 = l(u)`. `l` is bounded on H¹ because `w` is bounded.

Assumptions: (A1) `a` is coercive on `V0` (true: the displacement conditions exclude every rigid motion; Korn's
inequality); (A2) for P1 the quarter problem has the restriction of the full-plate solution as its solution (the
full problem is invariant under both reflections and its stress is unique, so its stress is symmetric; the free
vertical translation is fixed by `u_y = 0` on y = 0).

**Identity (PROVED).** Let `J_k(v) = (1/2) a(v,v) − k l(v)` and `m(k) = min_{v∈V} J_k(v)`. Let `z ∈ V0` solve
`a(z, w) = l(w)` for all `w ∈ V0` (the adjoint problem; Lax–Milgram by A1). The minimiser of `J_k` is `u + k z`
(its first variation `a(u + kz, w) − k l(w) = 0` for `w ∈ V0`), and with `a(u, z) = 0` and `l(z) = a(z, z)`:

    m(k) = (1/2) a(u,u) − k l(u) − (k^2/2) a(z,z),   hence   Q2 = l(u) = (m(−k) − m(k)) / (2k)   for every k > 0.

**Upper bound on m(k) (PROVED).** `m(k) <= J_k(v)` for every `v ∈ V`.

**Lower bound on m(k): how the adjoint load enters the complementary principle (PROVED).** The load `k l` acts like a
prestress `k s` (equivalently an eigenstrain `k w e_x⊗e_x`, since `S s = w e_x⊗e_x`). Apply the pointwise inequality
of section 1 with `σ = τ + k s`, `τ` statically admissible:

    (1/2) ε:C:ε − k s:ε  >=  (τ + k s):ε − (1/2)(τ + k s):S:(τ + k s) − k s:ε  =  τ:ε − (1/2)(τ + k s):S:(τ + k s).

Integrating and using `∫ τ:ε(v) = W(τ)`, and `s:S:τ = w τ_xx`, `s:S:s = w^2 C_1111` (`C_1111 = 100/91`):

    m(k)  >=  W(τ) − (1/2) Q(τ) − k l(τ) − (k^2/2) C_1111 / |P|,      l(τ) = (1/|P|) ∫_P τ_xx.

The prestress needs no equilibrium of its own: the admissible set for `τ` is unchanged, the prestress only shifts
the stress inside the complementary energy. Equality holds for `τ = C ε(u + kz) − k s`, which is statically
admissible because `a(u + kz, w) = k l(w)` for `w ∈ V0`.

**Two-sided bound (PROVED).** With `v_±` admissible displacements and `τ_±` admissible stresses,

    (m_low(−k) − m_up(+k)) / (2k)  <=  Q2  <=  (m_up(−k) − m_low(+k)) / (2k).

If `v_± ≈ u ± k z` and `τ_± ≈ σ(u ± kz) ∓ k s` with energy errors `η_u^2` (primal problem) and `k^2 η_z^2` (adjoint),
the width is about `η_u^2/k + k η_z^2`, smallest (`2 η_u η_z`) at `k = η_u/η_z`. The builder chooses `k` from float
estimates of the two gaps and rounds it to a dyadic rational; any `k > 0` gives a valid bound.

## 3. P1: plate with a hole

**Model.** Quarter `Omega_q = [0,1]^2 minus the disc of radius 1/2`, mirror lines x = 0, y = 0, grip x = 1, free edge
y = 1 and free arc (definitions.py, `P1_QUARTER`). Q2 is the average of `σ_xx` over the half patch
`P_h = {1/2 <= r <= 11/20, π/2 − 1/10 <= θ <= π/2}`, `|P_h| = 21/8000` (by A2 the same as over the full patch).

**Displacements.** `v = (x/2, 0) + tensor-product B-splines` (degree `pu`, open knot vectors on [0,1]) with the first
and last x-coefficients of `v_x` and the first y-coefficients of `v_y` zero. The spline is defined on the whole
square; only its restriction to `Omega_q` matters.

**Stresses.** `τ = g(y) e_x⊗e_x + Airy(ρ^2 Ψ)`:
- `ρ = (x^2 + y^2 − 1/4)(7/4 + y^2 − x^2) = 2x^2 − x^4 + y^4 + (3/2)y^2 − 7/16`. The second factor is >= 3/4 on the
  square, so `ρ = 0` exactly on the hole. `ρ^2 Ψ` and its gradient vanish there, so the Airy stress has zero traction
  on the hole for every `Ψ` (on `ρ = 0` the Hessian of `ρ^2 Ψ` is `2Ψ ∇ρ⊗∇ρ`, and the Airy stress of it is
  `2Ψ ∇ρ^⊥⊗∇ρ^⊥`, orthogonal to the normal `∇ρ`). No approximation of the circle occurs.
- `ρ_x = 4x(1 − x^2)` vanishes on x = 0 and x = 1, and `ρ_y = 4y^3 + 3y` on y = 0. Hence `Ψ_x = 0` on x = 0, 1 and
  `Ψ_y = 0` on y = 0 give `τ_xy = 0` there (mirror lines and frictionless grip); `Ψ = Ψ_y = 0` on y = 1 gives zero
  traction on the top edge. These are equalities and zeros between B-spline coefficients.
- `g(y) = Σ γ_k (y − 1/2)^{e_k}` for y >= 1/2, 0 below, carries the force (`F_q = Σ γ_k (1/2)^{e_k+1}/(e_k+1)`). It
  depends on y only (so it is divergence-free), it is zero near the hole, and across y = const lines it is tangential.
  The Airy part carries no force through the grip (`[∂_y(ρ^2Ψ)](1, y)` vanishes at y = 0 and y = 1).
- Why two or more free terms in `g` (affects only sharpness, not validity): in the quarter the Airy function on the
  hole is `A + B y`, where `B` is the force and `A` the moment of `σ_xx` over the ligament x = 0 about the hole. The
  term `ρ^2 Ψ` cannot change `A` or `B`; `g` with two free coefficients can, so the family contains the exact stress
  in its closure. Exponents near `ps` keep `g` smooth at y = 1/2, where `ρ ~ x^2` on the line y = 1/2.

**Integrals.** Every integrand is a polynomial on each cell; it is built with rational coefficients (python-flint
`fmpq_mpoly`; stored floats are exact binary rationals) and integrated exactly. On cells cut by the circle, the part
inside the disc is subtracted; it is integrated in y and then along the arc, which reduces to `∫ cos^n t dt` with
endpoints whose cosine and sine are rationals or square roots of rationals; the endpoint angles are enclosed in Arb
balls (320 bits). The patch integrals use the divergence theorem on `∂P_h` (two arcs, the ray at `π/2 − 1/10`, the
segment on x = 0): `∫_{P_h} σ_xx(v) = ∮ (C_11 v_x n_x + C_12 v_y n_y)` and `∫_{P_h} τ_xx = ∮ Φ_y n_y` with `Φ` the
full Airy function of `τ` (C¹; `Φ_y = 0` on the hole). The pieces are split where they cross cell lines; the order of
the break points and the cell of each piece are decided with balls (an undecidable case raises an error). The only
inputs that are not exact rationals are `sin(1/10)`, `cos(1/10)`, `π` and the crossing angles, all Arb balls.

**Finding the fields (not trusted).** Float assembly with Gauss rules (cut cells by a θ-parametrised rule), sparse LU
with a tiny regularisation and iterative refinement; one solve for the primal and the dual problem at `k = 0` and one
for the adjoint parts, then `v_± = v_0 ± k v_1`, `τ_± = τ_0 ± k τ_1`. Float errors only make a bound weaker.

**Results (`certificates/P1_frictionless_N16_p9`: 16 x 16 cells, stress degree 9, displacement degree 8,
`g` exponents 7, 8, 9, `k = 2^-21`; 15 s to build).**

| Quantity | Rigorous bounds (builder) | Relative width | Claim (weakened, FORMAT §6) | Label |
|---|---|---|---|---|
| Q1 = K | 0.591200987507 <= K <= 0.591200987579 | 1.2e-10 | [0.5912009869, 0.5912009882] | CERTIFIED (checked by the agent only; independent checker pending) |
| Q2 | 1.0327375805 <= Q2 <= 1.0327680269 | 2.9e-5 | [1.032737550, 1.032768058] | CERTIFIED (checked by the agent only; independent checker pending) |

Both meet the acceptance criterion (1e-3) and the stretch goal (1e-4). The width of Q2 is set by the adjoint
problem, whose solution has stress jumps across the curved boundary of the patch that the smooth spline spaces do
not follow (adjoint energy gap about 25 against an adjoint energy of about 125-150); the primal energy gap is about
6e-12.

**Tests (`builder/tests_p1.py`, all pass, 90 s).**
1. Known answers, plate without hole: certified `1 <= K <= 1` to 1e-11 and `0.4999988 <= Q2 <= 0.5000012`.
2. Exact identities for random admissible stress fields (both grip variants): equilibrium, zero traction on the
   hole (each component of `τ∇ρ` divisible by `ρ`), on y = 1, zero shear on mirror lines and the frictionless grip,
   traction continuity across cell edges, zero grip force of the Airy part.
3. The exact evaluator and the float assembly agree on random fields (two independent code paths).
4. Falsification: perturbations of the certified fields of relative energy 0 to 1e-2, and other values of `k`, never
   give a K or Q2 bound that crosses the claimed interval.
5. Re-evaluation of the stored certificate reproduces the claims.

**What is guaranteed and what is not (P1).**
- Guaranteed, if the argument of sections 1–3 is right and python-flint's rational and ball arithmetic is correct:
  the K and Q2 of P1 (frictionless grips, as defined in problems/P1.md) lie in the claimed intervals. No
  floating-point rounding enters the certified numbers; the circle is not approximated.
- The argument is classical but has been checked by the agent only; the independent checker (a separate program,
  written from FORMAT.md) had not reported when these notes were written.
- Assumptions A1, A2 above; linear plane-stress elasticity; the weak (H¹) solution.
- Nothing is claimed about pointwise stresses, other patches, or the clamped variant (separate certificate).

## 4. P2: L-bracket

OPEN (next).
