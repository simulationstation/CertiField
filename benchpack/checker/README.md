# Independent checker (P1, P1c, P2, P3)

Re-checks a certificate of `benchpack/certificates/` (format `benchpack-certificate-1`, described in
`benchpack/certificates/FORMAT.md`). It shares no code with the program that made the certificate; the only shared
file is the problem definition `benchpack/problems/definitions.py`.

## Use

    cd benchpack/checker
    python check.py ../certificates/STEM [--procs 2] [--log logs/STEM.txt]        # P1, P1c, P2
    python p3check.py ../certificates/STEM [--procs 2] [--log logs/STEM.txt]      # P3 (FORMAT.md section 9)
    python tests.py

`check.py` prints every check with `[pass]` or `[FAIL]`, the recomputed integrals and bounds, and ends with
`RESULT: ALL CHECKS PASSED` (exit code 0) or the list of failures (exit code 1). Needs Python 3 with NumPy and mpmath.
It runs 2 processes by default (2 threads each). Per-cell results are saved in `cache/` (one folder per certificate,
keyed by the bytes of the certificate and of the checker's per-cell code); an interrupted run resumes where it stopped.
The P1 certificate `P1_frictionless_N16_p9` takes about a minute on one process.

`p3check.py` checks the new P3 certificates (FORMAT.md section 9; octant, degree-4 Beltrami potential) slab by slab
in two processes with a checkpoint after every slab (resumable). The older curved3d certificates of P3 have their own
exact checker, `Reference/Curved3D/exact_check.py`; for large grids `p3_exact.py` runs the same checks (a)-(h) with
that program's own functions, slab by slab, resumable; it gives the same result as `exact_check.py` on the grid-8
certificate:

    python p3_exact.py clamped ../../Reference/Curved3D/certificates/psi_clamped_N24.npz --v ../../Reference/Curved3D/certificates/v_clamped_N24.npz

## What is checked

For every displacement field `v` and stress field `tau` of the certificate, on every cell of its spline grid, the
field is rebuilt as a polynomial with rational coefficients (the stored float64 numbers are taken as the exact binary
rationals they are). Then, in exact rational arithmetic:

- equilibrium, `div tau = 0`, in every cell of the domain;
- continuity of the normal traction across every interior cell edge;
- P1: on every cell cut by the hole, both components of `tau grad(rho)` are divisible by `x^2 + y^2 - 1/4` and by
  `7/4 + y^2 - x^2` (so by `rho`), hence the traction on the hole is zero; zero traction on y = 1; zero shear on the
  mirror lines x = 0, y = 0 and (frictionless) on the grip x = 1;
- P2: zero traction on the four free edges;
- the grip force, integrated from the field and compared with the closed form of FORMAT.md (P1: the `g` part;
  P2: `-beta`);
- the displacement is continuous across every interior edge and has the grip, clamp and mirror values exactly;
- the coefficient conditions of FORMAT.md, the knot vectors (open, spanning the box, multiplicities), array shapes;
- the cells tile the box; each is classified as inside, cut or outside exactly; the areas of the covered domain and
  of the patch enclose `1 - pi/16` (P1 quarter), 3 (P2), `21/8000` and `1/100`.

Bounds (FORMAT.md sections 3.4 and 4.4) are recomputed from scratch:

- energies `a(v,v)`, `Q(tau)`, the force `F(tau)` and the patch averages `l(v)`, `l(tau)` on cells inside the domain:
  exact rationals (products of polynomials integrated with `int s^a t^b = 1/((a+1)(b+1))`, every cell, no sampling);
- on cells cut by the hole and on the patch `P_h`: the region is sliced in x; the integrals have closed forms
  (`int x^m sqrt(R^2 - x^2) dx` by an exact recursion, verified symbolically on each use); only the end-point values
  (square roots, `asin`, `sin(1/10)`, `cos(1/10)`) are evaluated, with `mpmath.iv` (outward rounding) at 1200 bits;
- the bounds `K >= F^2/Q`, `K <= 4 a_q` (P1) or `a` (P2), and the two Q2 quotients, in interval arithmetic.

P3 (`p3check.py`, with `exact3.py` and `ball.py`): the cells are the products of the knot spans (uniform `k/N`, or
graded towards x, y, z = 1 with the JSON key `grade`; the edges are checked against FORMAT.md section 9); the stress
`g(y,z) e_x e_x + inc(rho^2 Psi)` and the displacement are rebuilt per cell as trivariate rational polynomials; exact checks of equilibrium, traction continuity across all
faces meeting the domain, divisibility of `tau grad(rho)` by `rho` on every cell cut by the sphere, the free faces, the
mirror planes and the grip, the grip force, the displacement's continuity and boundary values, the coefficient
conditions. Energies: exact rationals on whole cells; on cut cells either over the octree boxes of the certificate
(superset mode, exact rationals) or over the cell minus the ball (exact mode): integrated in z, then the polynomial
parts over rectangles cut by circles (closed forms as above) and the part with `sqrt(R^2 - x^2 - y^2)` exactly in y
and, by parts, in x (square roots, asin; verified closed forms); end points in `mpmath.iv`.

A claim is confirmed if `lower_claim <= recomputed lower` and `recomputed upper <= upper_claim`.

## Tests (`tests.py`)

1. B-splines: exact reproduction of `x`; smoothness `C^(p-m)` across a knot of multiplicity `m`.
2. Curved integrals against closed forms (`1 - pi/16`, `1/3 - pi/256`, the area and first moments of `P_h`), on grids
   with grid points exactly on the circles; against `mpmath.quad` (diagnostic agreement, 1e-40).
3. Known answers: the plate without the hole, frictionless grips: `K = 1` exactly and the patch average `1/2`; the
   Q2 bracket of these fields is `[1/2 - kappa c/4, 1/2 + kappa c/4]` (`c = C11/|P_h|`) as derived by hand.
4. Green identity `int tau:eps(v) = (work factor) F(tau)` for random admissible fields: exact for P2, to 1e-100 for
   P1 frictionless and clamped. It needs equilibrium, all boundary conditions and the cut-cell integrals to be right.
5. Falsification: random admissible fields never give a lower bound above an upper bound (K and `m(k)`).
6. Negative tests: shear on a mirror line, a wrong grip value, a C^0 stress function, a discontinuous displacement,
   a stress field without the weight `rho^2`, claims tighter than the bounds, traction on a free edge of P2 — each
   is rejected.
7. P3: the solid cube with frictionless grips (K = 2 exactly); the Green identity on the solid cube for random fields
   (exact); octree box integrals against direct integration; random fields with the cavity; negative tests.
8. P3 exact mode: the closed forms of `ball.py` (exact identities); the moments of the octant of the ball against
   the Gamma-function values; the lift alone (`a_o = (lambda + 2 mu)(1 - pi/48)/4`); the Green identity over the
   octant minus the ball for random admissible fields, frictionless and clamped.

## Files

- `check.py` the checker (P1, P2); `p3check.py`, `exact3.py`, `ball.py` the P3 checker; `p3_exact.py` the resumable
  run of curved3d's exact checker; `fields.py` certificate data to per-cell polynomials and the per-cell work;
  `exact.py` exact polynomials, B-splines and integrals; `region.py` integrals over the curved regions;
  `tests.py`, `synth.py` (random test certificates), `ritz.py` (near-optimal test fields from a small float Ritz
  solver; only their coefficients come from floating point, the bounds are the checker's); `logs/` run logs.
- `logs/p3_exact_*_N24.txt`: P3, the separate exact checker of `Reference/Curved3D` (`exact_check.py`) run on both
  grid-24 certificates.

## Limits

- The argument that the bounds hold (complementary energy, the duality for Q2) is in `benchpack/NOTES.md`; this
  program checks the identities and recomputes the numbers, not the argument.
- The certificates are checked by the agent only.
