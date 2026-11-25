# Certificate format (P1, P2)

Version: `benchpack-certificate-1`, revision 2 (general base point `b` of `g`, the clamped variant P1c, exact
integration in section 7). This file describes the certificate data completely, so that the stress and
displacement fields can be rebuilt from the data alone and every identity and bound checked in exact arithmetic,
without reading the builder's code. The problem data (domain, material, grips, patch, the weight `rho` of P1, the
symmetry reduction of P1) are in `benchpack/problems/definitions.py`, the only file the checker may share with the
builder.

All quantities are dimensionless: E = 1, ν = 3/10, plane stress, unit thickness, δ = 1.

## 1. Files

One certificate = two files with the same stem in `benchpack/certificates/`:

- `STEM.json`: the problem, the spaces (degrees and knot vectors), the roles of the fields, the claims, and the
  builder's own enclosures of every intermediate number (diagnostic, for comparison only);
- `STEM.npz` (NumPy, `numpy.load`): the coefficient arrays, all `float64`.

Numbers:
- strings in JSON of the form `"p/q"` or `"p"` are exact rationals (Python `fractions.Fraction(s)` reads them);
- every `float64` stored in the `.npz` file or in JSON is taken as the exact binary rational it represents;
- claim strings are exact decimal strings.

## 2. Splines

A space is given by a degree `p` and two knot vectors `knots_x`, `knots_y` (JSON lists of rational strings). Each knot
vector is non-decreasing, open (its first and last knot repeated `p+1` times) and spans the box of the model
(P1: [0,1]; P2: [−1,1]). With `t` a knot vector of length `n+p+1`, the `n` B-splines `B_0, …, B_{n−1}` of degree `p`
are defined by the Cox–de Boor recursion

    B_i^0(s) = 1 if t_i <= s < t_{i+1} else 0,
    B_i^k(s) = (s − t_i)/(t_{i+k} − t_i) B_i^{k−1}(s) + (t_{i+k+1} − s)/(t_{i+k+1} − t_{i+1}) B_{i+1}^{k−1}(s)

(terms with a zero denominator are omitted). On every non-empty knot span `[t_k, t_{k+1})` each `B_i` is a polynomial
with rational coefficients; only these polynomial pieces matter. A knot of multiplicity `m` makes the spline
`C^{p−m}` across that line.

A tensor-product spline with coefficient array `A` of shape `(n_x, n_y)` is

    s(x, y) = sum_{i,j} A[i, j] B_i(x) B'_j(y),

`B_i` on `knots_x`, `B'_j` on `knots_y`. The **cells** are the rectangles `[t_k, t_{k+1}] × [t'_l, t'_{l+1}]` of
non-empty spans; on each cell every field below is a polynomial with rational coefficients.

Each certificate has one **displacement space** (`displacement_space`: degree `pu`, knots) shared by all its
displacement fields, and one **stress space** (`stress_space`: degree `ps`, knots, and for P1 the exponents of `g`)
shared by all its stress fields.

npz array names: a displacement field called `NAME` is stored as `NAME_A` and `NAME_B` (shape `(n_x, n_y)` of the
displacement space); a stress field called `NAME` as `NAME_C` (shape `(n_x, n_y)` of the stress space) and
`NAME_par` (P1: the coefficients `gamma_k` of `g`; P2: `(alpha, beta, gamma)`).

## 3. P1 (quarter model; variant `frictionless`, and `clamped` = P1c)

Model: `Omega_q = {(x,y) in [0,1]^2 : x^2 + y^2 >= 1/4}` (`P1_QUARTER` in definitions.py). Mirror lines x = 0 and
y = 0, grip x = 1 (frictionless or clamped), free edge y = 1, free arc `x^2 + y^2 = 1/4`. Where the two variants
differ this is stated.

### 3.1 Displacement fields

    v_x(x, y) = x/2 + sum A[i,j] B_i(x) B'_j(y),       v_y(x, y) = sum B[i,j] B_i(x) B'_j(y).

Conditions (exact equalities of stored floats): `A[0, :] = 0`, `A[n_x−1, :] = 0`, `B[:, 0] = 0`; for `clamped` also
`B[n_x−1, :] = 0`. With open knot vectors this gives `v_x = 0` on x = 0, `v_x = 1/2 = δ/2` on x = 1, `v_y = 0` on
y = 0 (and `v_y = 0` on x = 1 for `clamped`): every such `v` is kinematically admissible for the quarter model.

Strain and energy: `eps = (dv_x/dx, dv_y/dy, dv_x/dy + dv_y/dx)`, `a_q(v, v) = ∫_{Omega_q} eps.C.eps dA` with `C` of
definitions.py (`C = [[100/91, 30/91, 0], [30/91, 100/91, 0], [0, 0, 5/13]]`). Stress of `v`: `C eps(v)`.

### 3.2 Stress fields

    tau = g(y) e_x⊗e_x + Airy(rho^2 Psi),
    Psi(x, y) = sum C[i,j] B_i(x) B'_j(y)                (stress space, degree ps),
    g(y)      = sum_k gamma_k (y − b)^{e_k} for y >= b,   g(y) = 0 for y < b,          (b = g_base >= 1/2)
    rho(x, y) = 2x^2 − x^4 + y^4 + (3/2)y^2 − 7/16 = (x^2 + y^2 − 1/4)(7/4 + y^2 − x^2)   (RHO in definitions.py),
    Airy(phi): tau_xx = phi_yy, tau_yy = phi_xx, tau_xy = −phi_xy.

The exponents `e_k` (non-negative integers) are `stress_space.g_exponents`, the base point `b` (a rational >= 1/2)
is `stress_space.g_base`; `gamma_k = NAME_par[k]`. So on each cell `tau` is a polynomial; `g` adds a polynomial in y
on cells with y >= b (`b` is always a knot of `knots_y` of the stress space, so `g` is polynomial on every cell).

Conditions on the coefficients (exact equalities of stored floats):

- `C[0, :] = C[1, :]` (`Psi_x = 0` on x = 0); for `frictionless` also `C[n_x−1, :] = C[n_x−2, :]` (`Psi_x = 0` on
  x = 1); for `clamped` no condition at x = 1;
- `C[:, 0] = C[:, 1]` (`Psi_y = 0` on y = 0);
- `C[:, n_y−1] = C[:, n_y−2] = 0` (`Psi = Psi_y = 0` on y = 1).

Identities that make `tau` statically admissible (all are polynomial identities, to be checked exactly):

1. `div tau = 0` on every cell (true for any Airy field and for `g(y) e_x⊗e_x`).
2. Traction continuity across every interior cell edge: across x = const, `(tau_xx, tau_xy)` continuous; across
   y = const, `(tau_xy, tau_yy)` continuous (follows from `rho^2 Psi` being C^1: knot multiplicities are at most
   `ps − 1`; `g` can jump only across y = b, where it is tangential).
3. Hole: on every cell meeting the circle, both components of `tau ∇rho` are divisible by `rho` (as polynomials). For
   the Airy part this holds because `rho^2 Psi` and its gradient vanish on `rho = 0`; the `g` part vanishes on all
   cells with y < b, and the circle has no points with y > 1/2. Since `rho = 0` on the circle exactly, the
   traction on the hole is zero.
4. Free edge y = 1: `tau_xy = tau_yy = 0`. Mirror lines: `tau_xy = 0` on x = 0 and on y = 0. Frictionless grip:
   `tau_xy = 0` on x = 1 (no condition at x = 1 for `clamped`).
5. Force: `F_q(tau) = ∫_0^1 tau_xx(1, y) dy = sum_k gamma_k (1 − b)^{e_k+1}/(e_k+1)` (the Airy part contributes
   `[d(rho^2 Psi)/dy](1, y)` from y = 0 to 1, which is 0).

Complementary energy: `Q_q(tau) = ∫_{Omega_q} tau.S.tau dA` with `S = [[1, −3/10, 0], [−3/10, 1, 0], [0, 0, 13/5]]`,
i.e. `tau_xx^2 + tau_yy^2 − (3/5) tau_xx tau_yy + (13/5) tau_xy^2`. Integrals are over `Omega_q` exactly (not over a
superset); on cells cut by the circle the polynomial integrand is integrated over the part with `x^2 + y^2 >= 1/4`.

### 3.3 Patch functional

`P_h = {(r cos t, r sin t) : 1/2 <= r <= 11/20, pi/2 − 1/10 <= t <= pi/2}`, `|P_h| = 21/8000` (half of the patch of
P1.md; by symmetry Q2 is the average of `sigma_xx` over `P_h`).

    l(v)   = (1/|P_h|) ∫_{P_h} (C eps(v))_xx dA = (8000/21) ∫_{P_h} (100/91)(dv_x/dx + (3/10) dv_y/dy) dA,
    l(tau) = (8000/21) ∫_{P_h} tau_xx dA.

### 3.4 Bounds

Q1 (role `fields.Q1`: displacement `v`, stress `tau`):

    F_q(tau)^2 / Q_q(tau)  <=  K  <=  4 a_q(v, v).

Q2 (role `fields.Q2`: `kappa` > 0 exact rational, displacements `v_plus`, `v_minus`, stresses `tau_plus`,
`tau_minus`). With `m(k) = min_v [ (1/2) a_q(v,v) − k l(v) ]` over admissible `v` (see benchpack/NOTES.md for the
argument):

    m(k)  <=  (1/2) a_q(v, v) − k l(v)                                               for every admissible v,
    m(k)  >=  (1/2) F_q(tau) − (1/2) Q_q(tau) − k l(tau) − (k^2/2) (100/91) / |P_h|    for every admissible tau,
    Q2 = (m(−kappa) − m(kappa)) / (2 kappa).

The certificate uses `v_plus`, `tau_plus` at `k = +kappa` and `v_minus`, `tau_minus` at `k = −kappa`:

    m_up(+k)  = (1/2) a_q(v_plus, v_plus)  − kappa l(v_plus),
    m_up(−k)  = (1/2) a_q(v_minus, v_minus) + kappa l(v_minus),
    m_low(+k) = (1/2) F_q(tau_plus)  − (1/2) Q_q(tau_plus)  − kappa l(tau_plus)  − (kappa^2/2)(100/91)(8000/21),
    m_low(−k) = (1/2) F_q(tau_minus) − (1/2) Q_q(tau_minus) + kappa l(tau_minus) − (kappa^2/2)(100/91)(8000/21),

    (m_low(−k) − m_up(+k)) / (2 kappa)  <=  Q2  <=  (m_up(−k) − m_low(+k)) / (2 kappa).

## 4. P2 (L-bracket; full domain, no symmetry)

Model: `Omega = [−1,1]^2` minus `(0,1] × [−1,0)` (closure of the L-bracket of definitions.py). Both spaces span
[−1,1]^2; cells inside `[0,1] × [−1,0]` are not part of `Omega` (0 is a knot in both directions, so every cell is
either in `Omega` or outside it).

### 4.1 Displacement fields

    v_x = sum A[i,j] B_i(x) B'_j(y),        v_y = −(x + 1)/2 + sum B[i,j] B_i(x) B'_j(y).

Conditions: `A[0, :] = A[n_x−1, :] = 0`, `B[0, :] = B[n_x−1, :] = 0`. So `v = 0` on x = −1 and `v = (0, −1)` on
x = 1. `a(v, v) = ∫_Omega eps.C.eps dA`.

### 4.2 Stress fields

    tau = Airy(phi_p + Psi),   phi_p = chi(y) (alpha + beta x + gamma y),
    chi(y) = 1 for y <= 0,   chi(y) = 1 − 3y^2 + 2y^3 for 0 <= y <= 1,
    Psi = sum C[i,j] B_i(x) B'_j(y)    (stress space), (alpha, beta, gamma) = NAME_par.

`chi` is C^1 with `chi(1) = chi'(1) = chi'(0) = 0`; 0 is a knot of `knots_y`, so `tau` is polynomial on every cell.

Index classes of the stress space: in x, `L_x = {i : support of B_i ⊂ [−1, 0]}`, `R_x = {i : support ⊂ [0, 1]}`,
`S_x` = the others (support contains 0 in its interior); `L_y`, `R_y`, `S_y` likewise in y. The knot 0 has
multiplicity `ps − 1` in both knot vectors (so `|S_x| = |S_y| = 2`).

Conditions on `C` (stored zeros):

- `C[i, j] = 0` if `j in {n_y−2, n_y−1}`                           (`Psi = Psi_y = 0` on y = 1);
- `C[i, j] = 0` if `i in L_x ∪ S_x` and `j in {0, 1}`               (`Psi = Psi_y = 0` on y = −1, x <= 0);
- `C[i, j] = 0` if `i in S_x` and `j in L_y ∪ S_y`                  (`Psi = Psi_x = 0` on x = 0, y <= 0);
- `C[i, j] = 0` if `i in R_x` and `j in S_y ∪ L_y`                  (`Psi = Psi_y = 0` on y = 0, x >= 0, and
  functions supported outside `Omega`).

Identities to check exactly: `div tau = 0` on every cell of `Omega`; traction continuity across interior cell edges;
zero traction on the four free edges of definitions.py (`P2['free']`); and the force

    F(tau) = −∫_0^1 tau_xy(1, y) dy = [phi_x(1, y)]_{y=0}^{y=1} = −beta.

`Q(tau) = ∫_Omega tau.S.tau dA`.

### 4.3 Patch functional

`S_p = [0, 1/10]^2`, `|S_p| = 1/100` (1/10 is a knot of every knot vector):
`l(v) = 100 ∫_{S_p} (C eps(v))_xx dA`, `l(tau) = 100 ∫_{S_p} tau_xx dA`.

### 4.4 Bounds

Q1: `F(tau)^2 / Q(tau) <= K <= a(v, v)`.

Q2: as in 3.4 with `m(k) = min_v [ (1/2) a(v,v) − k l(v) ]`,

    m_up(±k)  = (1/2) a(v_±, v_±) ∓ kappa l(v_±),
    m_low(±k) = F(tau_±) − (1/2) Q(tau_±) ∓ kappa l(tau_±) − (kappa^2/2)(100/91)(100),

and the same two quotients for Q2.

## 5. JSON keys

    format              "benchpack-certificate-1"
    problem, variant    "P1" "frictionless" | "P1" "clamped" (P1c) | "P2" "standard"
    model               "quarter" (P1) | "full" (P2)
    displacement_space  {"degree": pu, "knots_x": [...], "knots_y": [...]}
    stress_space        {"degree": ps, "knots_x": [...], "knots_y": [...], "g_exponents": [...], "g_base": "b"} (P1)
    fields              {"Q1": {"v": NAME, "tau": NAME},
                         "Q2": {"kappa": "p/q", "v_plus": NAME, "v_minus": NAME, "tau_plus": NAME, "tau_minus": NAME}}
    claims              {"Q1": {"lower": "...", "upper": "..."}, "Q2": {"lower": "...", "upper": "..."}}
    computed            the builder's rigorous enclosures [lo, hi] (floats) of a_q, Q_q, F_q, l for every field and of
                        the m bounds; diagnostic only
    label               "CERTIFIED (checked by the agent only)" or another label of the ground rules

A certificate may contain only `Q1` (then `fields.Q2` and `claims.Q2` are absent).

## 6. Claims

The builder computes rigorous bounds `lo <= Q <= hi` (exact rationals, and Arb ball arithmetic for the
transcendental constants of P1; see section 7). The claim is weakened before it is stored:
with `s = max(|lo|, |hi|)` and `w = max(1e-9 s, 1e-3 (hi − lo))`, the claimed lower number is `lo − w` rounded down
and the claimed upper number `hi + w` rounded up, both to 10 significant decimal digits. A checker confirms a claim
if its own rigorous bounds satisfy `lower_claim <= lo_checker` and `hi_checker <= upper_claim`.

## 7. How the builder evaluates the integrals (informational)

- On every cell the integrand is built as one polynomial with rational coefficients (python-flint `fmpq_mpoly`; the
  stored floats are exact binary rationals) and integrated over the cell exactly (rational result). No quadrature
  rule and no floating-point rounding are involved.
- P1 cells cut by the circle (`min x^2+y^2 over the cell < 1/4 < max`, decided exactly): the exact integral over the
  cell minus the integral over the part inside the disc. The latter is integrated in y first, then along the arc,
  which reduces to `∫ cos^n t dt` between endpoints whose cosines and sines are rationals or square roots of
  rationals; the only transcendental numbers are the endpoint angles. They are enclosed in Arb balls (python-flint,
  320 bits). Cells inside the disc are skipped.
- P1 patch integrals: divergence theorem on `P_h` (`∫ dv_x/dx = ∮ v_x n_x`, `∫ dv_y/dy = ∮ v_y n_y`,
  `∫ tau_xx = ∮ (phi_y) n_y` with `phi` the full Airy function of `tau`, which is C^1 and has `phi_y = 0` on the
  hole); the boundary pieces are the two arcs, the ray at angle `pi/2 − 1/10` and the segment on x = 0, split at the
  cell lines; trigonometric constants (`sin(1/10)`, `cos(1/10)`, crossing angles) as Arb balls.
- P2: every cell is inside or outside `Omega`; exact rational integrals as above.

## 8. What a checker has to do

1. Read the data; rebuild every field as one polynomial with rational coefficients per cell.
2. Check the coefficient conditions (3.1, 3.2, 4.1, 4.2) exactly, and the identities 1–5 of 3.2 (P1) or of 4.2 (P2).
3. Compute rigorous enclosures of `a`, `Q`, `F`, `l` for every field over the exact domain (or patch).
4. Form the bounds of 3.4 / 4.4 and compare with the claims (section 6).

## 9. P3 (block with a spherical cavity, 3D; octant model; only Q1 = K)

Problem data: `P3`, `P3_OCTANT`, `P3_LAME` in definitions.py. 3D isotropic elasticity, E = 1, ν = 3/10,
`λ = 15/26`, `μ = 5/13`; `C ε = λ tr(ε) I + 2μ ε`, `S τ = (1+ν) τ − ν tr(τ) I` (E = 1). Voigt order
`(xx, yy, zz, xy, yz, xz)`; tensor indices 0, 1, 2 = x, y, z.

**Model.** `Omega_o = {(x,y,z) in [0,1]^3 : x^2+y^2+z^2 >= 1/4}` (the whole octant; no wedge). Mirror planes x = 0,
y = 0, z = 0; grip x = 1; free faces y = 1, z = 1 and the sphere. All three coordinate directions use the same
cell size `1/N` (N in the JSON); cells `[e/N, (e+1)/N]` per axis, `e = 0..N−1`.

**Splines.** As in section 2, with three knot vectors (JSON `knots`: list of three lists, axes x, y, z). Tensor
product `s(x,y,z) = sum A[i,j,k] B_i(x) B'_j(y) B''_k(z)`.

**Displacement** (`V`, shape `(3, n_x, n_y, n_z)`, degree `pu`):

    v = (x/2) e_x + (s_0, s_1, s_2),   s_a = sum V[a, i, j, k] B_i(x) B'_j(y) B''_k(z).

Conditions (exact): `V[0, 0, :, :] = V[0, n_x−1, :, :] = 0` (v_x = 0 on x = 0, = 1/2 on x = 1);
`V[1, :, 0, :] = 0` (v_y = 0 on y = 0); `V[2, :, :, 0] = 0` (v_z = 0 on z = 0); clamped grips also
`V[1, n_x−1, :, :] = V[2, n_x−1, :, :] = 0`. The builder also imposes the y↔z symmetry of the problem:
`V[0, i, j, k] = V[0, i, k, j]` and `V[2, i, j, k] = V[1, i, k, j]` (this needs the same knots in y and z); a
checker need not rely on it (it is a restriction, not an admissibility condition).

**Stress** (`Psi`, shape `(6, n_x, n_y, n_z)`, degree `ps`; `gam`, length `2 n_g`):

    tau = g(y, z) e_x⊗e_x + inc(Phi),   Phi = rho^2 Psi  (a symmetric tensor field),
    Psi_c = sum Psi[c, i, j, k] B_i(x) B'_j(y) B''_k(z),  c = 0..5 for the components (xx, yy, zz, xy, yz, xz),
    inc(Phi)_ij = sum_{k,l,m,n} e_ikl e_jmn d_k d_m Phi_ln   (e the Levi-Civita symbol, d_k = d/dx_k),
    rho = eta(x) + y^2 + z^2 − 1/4,  eta(x) = x^2 for x <= 3/4,  eta(x) = x^2 − 4 (x − 3/4)^2 for x >= 3/4,
    g(y, z) = sum_k gamma_k (y − 1/2)^{e_k} [y >= 1/2] + sum_k delta_k (z − 1/2)^{e_k} [z >= 1/2],

with `(gamma_1..gamma_ng, delta_1..delta_ng) = gam` and `e_k = stress_space.g_exponents`. 1/2 is a knot in y and z and
3/4 a knot in x, so every field is one polynomial per cell. On the cube `rho` vanishes exactly on the sphere
(`eta >= 9/16 > 1/4` for x >= 3/4), is negative inside it and positive outside; `eta` is C¹ with `eta'(0) = eta'(1) = 0`.
The knot 3/4 has multiplicity `ps − 1` in the x-knots of the stress space (Psi is only C¹ there, like rho^2).

Conditions on `Psi` (exact). For component c = (a, b) and axis k let `s = (−1)^{[a=k] + [b=k]}` (+1: even, −1: odd):
- at the plane x_k = 0 (first index along axis k): even: index 0 equals index 1; odd: index 0 is zero;
- at x = 1 (last index along x), frictionless grips only: even: index n_x−1 equals index n_x−2; odd: index n_x−1 is
  zero; clamped grips: no condition at x = 1;
- at y = 1 and at z = 1 (last two indices along y, resp. z): zero, for all six components.
The builder also imposes the y↔z symmetry `Psi[swap(c), i, k, j] = Psi[c, i, j, k]` with swap = (xx→xx, yy→zz,
zz→yy, xy→xz, yz→yz, xz→xy), and `delta_k = gamma_k` (again a restriction a checker need not rely on).

Identities that make `tau` statically admissible (polynomial identities, cell by cell):
1. `div tau = 0` (true for any inc field and for `g(y,z) e_x⊗e_x`).
2. Traction continuity across cell faces: across x_k = const, `tau_kj` continuous for all j (Phi is C¹; the jumps of
   `g` across y = 1/2 or z = 1/2 are tangential).
3. Sphere: on every cell meeting it, each component of `tau ∇rho` is divisible by `rho`; `g` vanishes on every
   cell with y < 1/2 and z < 1/2, and the sphere has no point with y > 1/2 or z > 1/2.
4. Free faces y = 1, z = 1: zero traction. Mirror planes: x = 0: `tau_xy = tau_xz = 0`; y = 0: `tau_xy = tau_yz = 0`;
   z = 0: `tau_xz = tau_yz = 0`. Frictionless grip x = 1: `tau_xy = tau_xz = 0`.
5. Force through the grip: `F_o = ∫∫_{[0,1]^2} tau_xx(1, y, z) dy dz = sum_k (gamma_k + delta_k) (1/2)^{e_k+1}/(e_k+1)`
   (the inc part carries no force).

**Bounds.** `a_o(v) = ∫ eps(v):C:eps(v)`, `Q_o(tau) = ∫ tau:S:tau` over `Omega_o`;

    2 F_o(tau)^2 / Q_o(tau)  <=  K  <=  8 a_o(v).

(Solid cube, frictionless: `tau = (1/2) e_x⊗e_x`, `F_o = 1/2`, `Q_o = 1/4`, `v = (x/2, −νy/2, −νz/2)`, `a_o = 1/4`:
both give K = 2.)

**Integration domain** (JSON `integration`). Both integrands are non-negative, so a checker may integrate over the
exact domain or over any superset; the claims hold for the builder's choice, which is stated:
- `{"mode": "superset", "stress_depth": L, "displacement_depth": Lu}`: the union of the boxes of the octree rule
  below (depth L for `Q_o`, Lu for `a_o`), which contains `Omega_o`. Rule (curved3d's `build_boxes`): for each cell
  `(e_x, e_y, e_z)`, start with the cell (level 0). A box at level l has integer corner `lo = e 2^l + j` and side
  `1/(N 2^l)`; with `r_min = sum lo_k^2`, `r_max = sum (lo_k + 1)^2` and `t = (N 2^l)^2 / 4`: if `r_max <= t` the
  box is dropped (it lies in the closed ball); if `r_min >= t` it is kept whole; otherwise it is kept whole if
  l = L, else split into its 8 children at level l + 1.
- `{"mode": "exact"}`: `Omega_o` itself (cells cut by the sphere integrated over their part outside the ball).

JSON keys for P3: `problem "P3"`, `variant`, `model "octant"`, `N`, `displacement_space {"degree", "knots"}`,
`stress_space {"degree", "knots", "g_exponents", "g_base": "1/2", "rho", "eta_flat_from": "3/4"}`, `integration`,
`claims {"Q1": {"lower", "upper"}}` (section 6), `computed` (diagnostic), `label`; npz arrays `V`, `Psi`, `gam`.
