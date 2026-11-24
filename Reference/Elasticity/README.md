# Linear elasticity: reference implementation (Python)

Certified two-sided bounds for

    E[u] = ∫_Ω ½ e(u)ᵀ D e(u) − b·u dx,     Ω a rectangle, plane strain, e(u) = (∂₁u₁, ∂₂u₂, ∂₁u₂ + ∂₂u₁),
    u = g on the Dirichlet faces, no traction on the others,

with a constant symmetric positive definite stiffness D (rational entries) and a constant body force b. Theory and
proofs: `docs/ELASTICITY.md`.

This is the reference implementation. `Kernel/Elasticity.wl` is its Wolfram Language port (see `PORTING.md`); the
port has its own tests (`Tests/Elasticity.wlt`), reproduces the test vectors stored here and recomputes every
certificate in `results/`.

Requirements: Python 3.10 or later, NumPy, SciPy. Nothing else. Developed and run on Windows 11 (x86-64) with
Python 3.11.9, NumPy 2.4.6, SciPy 1.17.1. Not yet run on Linux or macOS. Run the scripts from this directory.

## Files

| File | What it is |
|---|---|
| `elasticity.py` | the implementation: problem, exact local matrices, exact nodal values on the traction-free faces, candidate and certified upper bound, stress field and certified lower bound, certificates on disk |
| `ia.py` | outward-rounded interval arithmetic on NumPy arrays (the file of `Reference/ScalarND`) |
| `run_examples.py` | the examples of `docs/ELASTICITY.md`; writes `results/NAME.json` + `.npz` (seconds) |
| `verify_all.py` | recomputes both bounds of every stored certificate from its exact data |
| `check_exact.py` | rebuilds both fields as polynomials in rational arithmetic: static admissibility exactly, both exact values against the enclosures, the hypercircle identity |
| `check_theorem.py` | minima known as rational numbers, and falsification with random admissible fields (about 15 s) |
| `make_testvectors.py` | small inputs and expected outputs for the Wolfram Language port (`testvectors/`) |
| `results/` | stored certificates (`NAME.json`: the problem in exact form and the claimed bounds; `NAME.npz`: boundary data, candidate, stress field) |
| `testvectors/` | inputs and expected outputs for `Tests/Elasticity.wlt` (70 kB) |
| `PORTING.md` | the Wolfram Language port: state, function map, conventions, differences, log, observations |

## Use

```python
from fractions import Fraction as Fr
import elasticity as E

# a block clamped on the left face, pulled on the right face, no traction on the other two
P = E.Problem(box=[(0, 1), (0, 1)], N=(32, 32), lam=Fr(3, 2), mu=Fr(1), degree=2,
              bc=(('D', 'D'), ('N', 'N')), g=lambda x, y: (0.1 * x, 0 * y))
out, U, st = E.certify(P)                    # out['lower'] <= inf E <= out['upper']

E.upper_bound(P, U)                          # (lo, hi): enclosure of the energy of any nodal array U
E.lower_bound(P, st)                         # the bound of any stress field st (E.Stress)
```

`bc[i] = (type at the lower face, type at the upper face)` in direction i, `'D'` (displacement prescribed) or `'N'`
(no traction). `D=` takes a 3×3 stiffness matrix instead of `lam`, `mu`; `b=` a constant body force. Any nodal
array with the Dirichlet values imposed is an admissible candidate, and any nodal values and parameters give a
statically admissible stress field: both evaluators are rigorous whatever they are given. The two linear solves
that choose the fields (`solve_candidate`, `solve_stress`) are ordinary floating point and are never trusted.

**The certified datum.** The Dirichlet datum of the certified problem is the trace of the continuous piecewise
polynomial of degree `degree` (1 or 2) in each variable through the nodal values of `g` on the boundary grid. Each
stored float is taken as the exact binary number it is.

Run the checks (`python3` instead of `python` where that is the name):

    python verify_all.py        # ALL CERTIFICATES VERIFIED
    python check_exact.py       # ALL EXACT CHECKS PASSED
    python check_theorem.py     # ALL TESTS PASSED

## Stored certificates

| Name | Problem | Grid, degree | Claim |
|---|---|---|---|
| `block_p1_16` | block clamped on the left, pulled on the right, free sides | 16×16, 1 | 0.01501720485 ≤ inf E ≤ 0.01506421420 |
| `block_p1_32` | the same | 32×32, 1 | 0.01502346570 ≤ inf E ≤ 0.01504167619 |
| `block_p2_16` | the same | 16×16, 2 | 0.01501720485 ≤ inf E ≤ 0.01503403927 |
| `block_p2_32` | the same | 32×32, 2 | 0.01502346570 ≤ inf E ≤ 0.01502973919 |
| `bending` | pure bending with free sides; the exact minimum is 1/32 | 4×2, 2 | 0.03124999996 ≤ inf E ≤ 0.03125000004 |
| `gravity` | column clamped at the bottom under its own weight, three free faces | 16×32, 2 | −0.4514566136 ≤ inf E ≤ −0.4511230036 |

The claims are the computed bounds weakened by one part in 10⁹ of the larger modulus and rounded outward to ten
significant digits (`claim_strings`).
