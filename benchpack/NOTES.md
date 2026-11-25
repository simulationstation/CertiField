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

### 3.1 P1c: clamped grips

Same model, with `u_y = 0` added on the grip (the last x-coefficients of `v_y` are zero) and without the condition
`Ψ_x = 0` on x = 1 (the clamped grip may carry shear). The grip work is still `(1/2) F_q`, so all formulas of
sections 2 and 3 hold unchanged. Where the clamped grip meets the free edge y = 1 the stress is singular; the knots
are graded towards x = 1 and y = 1 (3 halvings) and towards the ray and the outer arc of the patch (7 halvings on both
sides of x = 21/400 and y = 11/20), where the adjoint stress jumps.

**Result (`certificates/P1c_N32_p7`: base 32 x 32 cells, 49 x 49 after grading, stress degree 7, displacement degree
6, `k = 5/8192`; 131 s).**

| Quantity | Rigorous bounds (builder) | Relative width | Claim (weakened) | Label |
|---|---|---|---|---|
| Q1 = K | 0.595719616964 <= K <= 0.595722280058 | 4.5e-6 | [0.5957196143, 0.5957222828] | CERTIFIED (checked by the agent only; independent checker pending) |
| Q2 | 1.0450481089 <= Q2 <= 1.0460953984 | 1.0e-3 | [1.045047061, 1.046096446] | CERTIFIED (checked by the agent only; independent checker pending) |

P1c has no acceptance criterion of its own; measured against those of P1, K meets 1e-4 and Q2 sits at 1.0e-3.
Clamping raises K by 0.77 % and the patch average by 1.2 % against the frictionless grips. Deeper grading towards the
clamped corner (4 or more halvings) made the float stress solve lose accuracy (elongated cells, as for P2), so the
corner limits the width of Q2 here. What is guaranteed and what is not: as for P1, for the clamped variant.

## 4. P2: L-bracket

**Model.** The full L-shaped domain (no symmetry). Both spaces are tensor-product B-splines on [−1,1]^2; cells in
the removed square are not integrated, so only the restriction to `Omega` matters.

**Displacements.** `v = (0, −(x+1)/2) + splines` with the first and last x-coefficients of both components zero:
`v = 0` on x = −1 and `v = (0, −1)` on the grip x = 1.

**Stresses.** `τ = Airy(φ_p + Ψ)`:
- `Ψ` has the knot 0 with multiplicity `ps − 1` in both directions (C¹ across x = 0 and y = 0), and every B-spline
  whose value or gradient is non-zero somewhere on a free edge has coefficient zero (index classes of FORMAT.md
  4.2). Then `Ψ` and its gradient vanish on all four free edges, so its traction is zero there, and `Ψ_x = 0` on
  y = 1 and on y = 0 (x >= 0), so it carries no vertical force through the grip. With the C¹ knot at 0 the zero
  conditions remove exactly the functions that touch the free edges (the two straddling B-splines are the only ones
  with value or slope on the line), not more.
- `φ_p = χ(y)(α + βx + γy)` with `χ = 1` for y <= 0 and `χ = 1 − 3y² + 2y³` on [0,1] (C¹; `χ(1) = χ'(1) = χ'(0) = 0`).
  It is linear for y <= 0 (zero stress there), its traction vanishes on y = 0 and y = 1, and its vertical grip
  force is `F = [φ_x(1,y)]_{y=0}^{1} = −β`. Why three parameters (sharpness only): for any admissible stress the
  Airy function can be normalised to vanish with its gradient on the top edge; on the other free boundary piece
  (y = 0 for x > 0, x = 0 for y < 0, y = −1 for x < 0, which is connected) it is then a linear function whose three
  coefficients are the two force components and the moment transmitted through the grip; `φ_p` provides all three.
- Knots: geometric layers towards x, y ∈ {−1, 0, 1} (the re-entrant corner and the four corners where a clamped edge
  meets a free edge) and towards 1/10 (the patch edges, where the adjoint stress jumps); multiplicity `ps − 1` at
  1/10 in the stress space and `pu` at 1/10 in the displacement space, so the adjoint fields may jump there.

**Integrals.** Every cell is a rectangle inside or outside `Omega`, and the patch `[0, 1/10]^2` is a union of cells.
All integrals are exact rationals (no transcendental number occurs in P2), computed cell by cell with python-flint,
in chunks that are saved (`benchpack/work/STEM/`) so that an interrupted build resumes.

**Results (`certificates/P2_h8_L8_p6`: 89 x 89 cells, base spacing 1/8, 8 dyadic layers with 2 cells each towards
the re-entrant corner lines, 6 towards the clamped corners, 4 towards 1/10; stress degree 6, displacement degree 5;
`k = 3/256`; 55 s to build on 2 processes).**

| Quantity | Rigorous bounds (builder, exact rationals) | Relative width | Claim (weakened) | Label |
|---|---|---|---|---|
| Q1 = K | 0.122603934395 <= K <= 0.122610210236 | 5.1e-5 | [0.1226039281, 0.1226102166] | CERTIFIED (checked by the agent only; independent checker pending) |
| Q2 | −0.3097348053 <= Q2 <= −0.3092223558 | 1.7e-3 | [−0.3097353178, −0.3092218433] | CERTIFIED (checked by the agent only; independent checker pending) |

Both meet the acceptance criteria (K: 1e-3, Q2: 1e-2). Q2 is negative: the patch at the re-entrant corner is in
compression.

**Second certificate, stretch goal (`certificates/P2_h32_L7_p6`: 129 x 129 cells, base spacing 1/32, 7 layers
towards the corner lines, 5 towards the clamped corners, 4 towards 1/10; degrees 6/5; `k = 9/1024`; 168 s).**

| Quantity | Rigorous bounds (builder, exact rationals) | Relative width | Claim (weakened) | Label |
|---|---|---|---|---|
| Q1 = K | 0.122604442179 <= K <= 0.122607381577 | 2.4e-5 | [0.1226044392, 0.1226073846] | CERTIFIED (checked by the agent only; independent checker pending) |
| Q2 | −0.3095829569 <= Q2 <= −0.3093255745 | 8.3e-4 | [−0.3095832144, −0.3093253171] | CERTIFIED (checked by the agent only; independent checker pending) |

This meets the stretch goal for Q2 (1e-3). A finer base grid with fewer layers keeps the cells less elongated for
the same smallest cell, which the float solve needs (see below). The two certificates are independent; their
intervals overlap, as they must.

**Limit of the tensor-product grading (DIAGNOSTIC).** Deeper grading (10 or more layers) makes the cells far from
the corner very elongated (aspect ratio 2^L); the float solve of the stress-function problem then fails (in one run
the float lower bound of K exceeded the float upper bound, and the exact evaluation of such a field gives a poor but
valid bound). Validity is never affected, only sharpness.

**Tests (`builder/tests_p2.py STEM`).**
1. Known answers of the exact evaluator on the L-shaped domain: the lift alone has `a = 15/52` and `l = 0`; the
   particular stress alone has `Q`, `l`, `F` equal to the closed forms computed by sympy.
2. Exact identities for random admissible stress fields: equilibrium, zero traction on the four free edges,
   traction continuity across cell edges, vertical grip force `−β`.
3. The exact evaluator agrees with the float assembly on random fields.
4. Falsification: perturbations of the certified fields of relative energy 0 to 1e-2, at other values of `k`, never
   cross the claimed intervals.
5. Re-evaluation of the stored certificate reproduces the claims.

**What is guaranteed and what is not (P2).**
- Guaranteed, if the argument of sections 1, 2 and 4 is right and python-flint's rational arithmetic is correct:
  the K and Q2 of P2 (problems/P2.md) lie in the claimed intervals. All certified numbers are exact rationals.
- The argument has been checked by the agent only; the independent checker had not reported when this was written.
- Assumptions: A1 (the clamped edge and the grip exclude rigid motions), linear plane-stress elasticity, the weak
  (H¹) solution, which exists although the stress is singular at the re-entrant corner.
- Nothing is claimed about the stress at the corner itself (it is unbounded) or about other patches.

## 5. P3: block with a spherical cavity (new certificates; the stored curved3d certificates are unchanged)

**Model** (FORMAT.md section 9): the octant `[0,1]^3` minus the ball of radius 1/2, mirror planes x, y, z = 0, grip
x = 1, free faces y = 1, z = 1 and the sphere; `K = 8 a_o(u,u)`; bounds `2 F_o(tau)^2 / Q_o(tau) <= K <= 8 a_o(v)`
(section 1 with the octant's grip work `F_o / 2`).

**What changes against curved3d** (whose lower bound was the weak side: 1.745105 at 24 cells per half-edge while
its upper bound was 1.758277):
- Degree. The stress is a second derivative of the potential, so with the same spline degree the stress converges
  one order slower than the strain of the displacement. The potential now has degree 4 (curved3d: 2).
- Particular field. curved3d's step column (axial stress outside the square column `|y|, |z| < 1/2`) jumps at the
  column faces, and the potential `rho^2 Psi` cannot cancel the jump near the points where the column touches the
  sphere (`rho^2` vanishes there to second order); in 2D (P1) step-like particular fields stalled the lower bound
  in the same way. It is replaced by `g(y, z) = sum gamma_k (y − 1/2)_+^{e_k} + delta_k (z − 1/2)_+^{e_k}` with
  exponents 4, 5, 6, whose coefficients are optimised together with the potential.
- Weight and ties as in curved3d (`rho = eta(x) + y^2 + z^2 − 1/4` with the flattened `eta`, parity ties at the
  mirror planes and the frictionless grip, potential and normal derivative zero on the free faces); the knot 3/4,
  where `eta` is only C¹, has multiplicity `ps − 1`.
- The y ↔ z symmetry of the problem is imposed on both fields (half the unknowns); dense direct solves.

**Integration.** Both integrands are non-negative. Mode `superset` integrates over curved3d's octree boxes (a cut
box is kept whole), which overestimates both energies by an amount proportional to the box size near the sphere;
for the displacement energy this bias is large (DIAGNOSTIC: at 8 cells per half-edge the float upper bound is
1.7732, 1.7653, 1.7615 at depths 3, 4, 5, halving with each level). Mode `exact` integrates the cells cut by the
sphere over their part outside the ball: the polynomial is integrated in z; the remaining integrals over a
rectangle intersected with a disc are reduced to arc integrals (`cos^n`, Arb endpoints) for the polynomial part
and, for the part with `sqrt(a^2 − x^2 − y^2)`, integrated exactly in y (closed forms with sqrt and asin) and in x
by Arb's rigorous integration (acb.integral), with the branch cuts checked. Full cells use the interval Gauss rules
of curved3d's kind (exact for the degree).

*Cut-cell reduction (PROVED, elementary).* For a cell `C = R × [z0, z1]` in the first octant, `R = [x0,x1]×[y0,y1]`,
and a polynomial `p`, let `Q(x,y,z) = ∫_{z0}^{z} p dz'` and `Z = (a² − x² − y²)^{1/2}`. Then
`∫_{C ∩ ball} p = ∫_{R ∩ D(r1)} Q(x,y,z1) + ∫_{R ∩ (D(r0) \ D(r1))} Q(x,y,Z)`, `D(r)` the disc of radius `r` about
the origin, `r_i² = a² − z_i²` (`D(r1)` empty if `z1 >= a`), because for `(x,y)` in `R` the segment of the cell
above `(x,y)` meets the ball in `z0 <= z < min(z1, Z)`. Writing `Q(x,y,Z) = P + Z W` (using `Z² = a² − x² − y²`), the
polynomial parts are integrated over `R ∩ D(r)` in y and then along the arc `x = r cos t`, which gives the integrals
`∫ cos^n t dt` of section 3 between endpoints whose cosines and sines are rationals or square roots of rationals.
For `∫_{R ∩ D(r)} Z W`, the y-integral is exact: with `b² = a² − x²`, `s = (b² − y²)^{1/2}`,
`K_0 = (y s + b² asin(y/b))/2`, `K_1 = (b³ − s³)/3`, `K_j = −y^{j−1} s³/(j+2) + (j−1) b²/(j+2) K_{j−2}` (integration
by parts) is `∫_0^y t^j (b² − t²)^{1/2} dt`; the remaining x-integral is enclosed by `acb.integral`, which returns
a ball containing the integral when the integrand reports non-analytic points (here: the branch cuts of the square
root and of asin, checked on every complex ball) — Arb's guarantee. Irrational interval ends are handled by
integrating to a nearby exact point and adding `|end − point| · sup|f|` over a ball around the end.

**Results (frictionless grips; 8 cells per half-edge, potential degree 4, displacement degree 3, exponents of `g`
4, 5, 6).**

| Certificate | Integration | Rigorous bounds (builder) | Relative width | Claim (weakened) | Label |
|---|---|---|---|---|---|
| `P3_frictionless_N8_p4_exact` | exact | 1.75716672 <= K <= 1.75762799 | 2.6e-4 | [1.757166257, 1.757628452] | CERTIFIED (checked by the agent only; independent checker pending) |
| `P3_frictionless_N8_p4` | superset, depths 4 / 6 | 1.75451382 <= K <= 1.75956610 | 2.9e-3 | [1.754508766, 1.759571150] | CERTIFIED (checked by the agent only; independent checker pending) |
| curved3d N = 24 (stored, unchanged) | superset | 1.74510549 <= K <= 1.75827651 | 7.5e-3 | | CERTIFIED (stored earlier) |

The exact certificate meets the target (0.1 %) by a factor of about four; build time 564 s on one process. The
three intervals are consistent (the exact one lies inside the other two). The superset certificate shows how
much of curved3d's width came from the integration domain: the same displacement space gives 1.75957 over the
superset at depth 6 and 1.75763 over the exact domain.

**Tests (`builder/tests_p3.py`, pass).** 1. Known answer: the lift alone has energy `(λ + 2μ)(1 − π/48)/4` over the
octant minus the ball; the computed upper bound exceeds it by 2e-13. 2. Exact identities of random admissible
stress fields (both grip variants): equilibrium, `tau ∇rho` divisible by `rho` on every cut cell, free faces, mirror
and grip shear, traction continuity, grip force formula. 3. Interval evaluator against the float assembly. The cut
cell integration agrees with mpmath on four cut cells to 1e-18. 4. With a certificate: perturbed stress fields
(ties restored) never give a lower bound above the claimed upper bound; 5. re-evaluation of the certificate.

**What is guaranteed and what is not (P3, new certificates).**
- Guaranteed, if the argument of section 1 (octant form) is right, python-flint's rational and ball arithmetic and
  Arb's rigorous integration are correct, and the outward-rounded float interval arithmetic of curved3d's `ia.py`
  is correct: the stiffness K of P3 (frictionless grips) lies in the claimed interval.
- The independent checker had not reported on these certificates when this was written; the argument has been
  checked by the agent only.
- Assumptions: linear 3D elasticity, the weak (H¹) solution, the octant symmetry (the problem is invariant under
  the three reflections; the stress is unique).
- The clamped variant is not covered by these new certificates (see below).
