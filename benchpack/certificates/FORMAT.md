# Certificate format (P1, P2)

Version: `benchpack-certificate-1`. This file describes the certificate data completely, so that the stress and
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

## 3. P1 (frictionless grips; quarter model)

Model: `Omega_q = {(x,y) in [0,1]^2 : x^2 + y^2 >= 1/4}` (`P1_QUARTER` in definitions.py). Mirror lines x = 0 and
y = 0, grip x = 1 (frictionless), free edge y = 1, free arc `x^2 + y^2 = 1/4`.

### 3.1 Displacement fields

    v_x(x, y) = x/2 + sum A[i,j] B_i(x) B'_j(y),       v_y(x, y) = sum B[i,j] B_i(x) B'_j(y).

Conditions (exact equalities of stored floats): `A[0, :] = 0`, `A[n_x−1, :] = 0`, `B[:, 0] = 0`. With open knot
vectors this gives `v_x = 0` on x = 0, `v_x = 1/2 = δ/2` on x = 1, `v_y = 0` on y = 0: every such `v` is
kinematically admissible for the quarter model.

Strain and energy: `eps = (dv_x/dx, dv_y/dy, dv_x/dy + dv_y/dx)`, `a_q(v, v) = ∫_{Omega_q} eps.C.eps dA` with `C` of
definitions.py (`C = [[100/91, 30/91, 0], [30/91, 100/91, 0], [0, 0, 5/13]]`). Stress of `v`: `C eps(v)`.

### 3.2 Stress fields

    tau = g(y) e_x⊗e_x + Airy(rho^2 Psi),
    Psi(x, y) = sum C[i,j] B_i(x) B'_j(y)                (stress space, degree ps),
    g(y)      = sum_k gamma_k (y − 1/2)^{e_k} for y >= 1/2,   g(y) = 0 for y < 1/2,
    rho(x, y) = 2x^2 − x^4 + y^4 + (3/2)y^2 − 7/16 = (x^2 + y^2 − 1/4)(7/4 + y^2 − x^2)   (RHO in definitions.py),
    Airy(phi): tau_xx = phi_yy, tau_yy = phi_xx, tau_xy = −phi_xy.

The exponents `e_k` (non-negative integers) are `stress_space.g_exponents`, the base point 1/2 is
`stress_space.g_base`; `gamma_k = NAME_par[k]`. So on each cell `tau` is a polynomial; `g` adds a polynomial in y on
cells with y >= 1/2 (`1/2` is always a knot of `knots_y` of the stress space, so `g` is polynomial on every cell).

Conditions on the coefficients (exact equalities of stored floats):

- `C[0, :] = C[1, :]` and `C[n_x−1, :] = C[n_x−2, :]` (so `Psi_x = 0` on x = 0 and on x = 1);
- `C[:, 0] = C[:, 1]` (`Psi_y = 0` on y = 0);
- `C[:, n_y−1] = C[:, n_y−2] = 0` (`Psi = Psi_y = 0` on y = 1).

Identities that make `tau` statically admissible (all are polynomial identities, to be checked exactly):

1. `div tau = 0` on every cell (true for any Airy field and for `g(y) e_x⊗e_x`).
2. Traction continuity across every interior cell edge: across x = const, `(tau_xx, tau_xy)` continuous; across
   y = const, `(tau_xy, tau_yy)` continuous (follows from `rho^2 Psi` being C^1: knot multiplicities are at most
   `ps − 1`; `g` jumps only across y = 1/2, where it is tangential).
3. Hole: on every cell meeting the circle, both components of `tau ∇rho` are divisible by `rho` (as polynomials). For
   the Airy part this holds because `rho^2 Psi` and its gradient vanish on `rho = 0`; the `g` part vanishes on all
   cells with y < 1/2, and the circle has no points with y > 1/2. Since `rho = 0` on the circle exactly, the
   traction on the hole is zero.
4. Free edge y = 1: `tau_xy = tau_yy = 0`. Mirror lines and grip: `tau_xy = 0` on x = 0, on y = 0 and on x = 1.
5. Force: `F_q(tau) = ∫_0^1 tau_xx(1, y) dy = sum_k gamma_k (1/2)^{e_k+1}/(e_k+1)` (the Airy part contributes
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
    problem, variant    "P1" "frictionless" | "P1" "clamped" | "P2" "standard"
    model               "quarter" (P1) | "full" (P2)
    displacement_space  {"degree": pu, "knots_x": [...], "knots_y": [...]}
    stress_space        {"degree": ps, "knots_x": [...], "knots_y": [...], "g_exponents": [...], "g_base": "1/2"} (P1)
    fields              {"Q1": {"v": NAME, "tau": NAME},
                         "Q2": {"kappa": "p/q", "v_plus": NAME, "v_minus": NAME, "tau_plus": NAME, "tau_minus": NAME}}
    claims              {"Q1": {"lower": "...", "upper": "..."}, "Q2": {"lower": "...", "upper": "..."}}
    computed            the builder's rigorous enclosures [lo, hi] (floats) of a_q, Q_q, F_q, l for every field and of
                        the m bounds; diagnostic only
    label               "CERTIFIED (checked by the agent only)" or another label of the ground rules

A certificate may contain only `Q1` (then `fields.Q2` and `claims.Q2` are absent).

## 6. Claims

The builder computes rigorous bounds `lo <= Q <= hi` (outward-rounded interval arithmetic, exact rationals and
ball arithmetic for the transcendental constants of P1; see section 7). The claim is weakened before it is stored:
with `s = max(|lo|, |hi|)` and `w = max(1e-9 s, 1e-3 (hi − lo))`, the claimed lower number is `lo − w` rounded down
and the claimed upper number `hi + w` rounded up, both to 10 significant decimal digits. A checker confirms a claim
if its own rigorous bounds satisfy `lower_claim <= lo_checker` and `hi_checker <= upper_claim`.

## 7. How the builder evaluates the integrals (informational)

- Cells entirely in `Omega` (P1: `min x^2+y^2 over the cell >= 1/4`, decided exactly): tensor Gauss–Legendre rules
  exact for the degree (n points per direction, `2n − 1 >=` the degree of the integrand in each variable), nodes and
  weights enclosed in rational intervals, values in outward-rounded float64 interval arithmetic, sums rounded up.
- P1 cells cut by the circle: exact polynomial integrand (rational coefficients) integrated over the cell minus the
  disc by Green's formula; the arc part reduces to `∫ cos^n t dt` between endpoints whose cosines and sines are
  square roots of rationals; the only transcendental numbers are the endpoint angles, enclosed with Arb balls
  (python-flint).
- P1 patch integrals: divergence theorem on `P_h` (`∫ dv_x/dx = ∮ v_x n_x`, `∫ dv_y/dy = ∮ v_y n_y`,
  `∫ tau_xx = ∮ (phi_y) n_y` with `phi` the full Airy function of `tau`, which is C^1 and has `phi_y = 0` on the
  hole); the boundary pieces are the two arcs, the ray at angle `pi/2 − 1/10` and the segment on x = 0, split at the
  cell lines; trigonometric constants (`sin(1/10)`, `cos(1/10)`, crossing angles) as Arb balls.
- P2: all cells are rectangles inside `Omega`; Gauss rules as above.

## 8. What a checker has to do

1. Read the data; rebuild every field as one polynomial with rational coefficients per cell.
2. Check the coefficient conditions (3.1, 3.2, 4.1, 4.2) exactly, and the identities 1–5 of 3.2 (P1) or of 4.2 (P2).
3. Compute rigorous enclosures of `a`, `Q`, `F`, `l` for every field over the exact domain (or patch).
4. Form the bounds of 3.4 / 4.4 and compare with the claims (section 6).
