# ScalarND reference implementation (Python)

Certified two-sided bounds for

    E[u] = ∫_Ω k |∇u|² + W(u) dx,     Ω a box in 2D or 3D (1D also works),
    u = g on the Dirichlet faces, free on the natural faces,

with a polynomial potential W that need not be convex. Theory and proofs: `docs/SCALAR_ND.md`.

This is the reference implementation. `Kernel/ScalarFieldND.wl` is its Wolfram Language port (see `PORTING.md`); the
port has its own tests (`Tests/ScalarFieldND.wlt`), reproduces the test vectors stored here and recomputes every
certificate in `results/`.

Requirements: Python 3.10 or later, NumPy, SciPy. Nothing else. Developed and tested with Python 3.13, NumPy 2.5,
SciPy 1.18 on Linux (x86-64); the three check scripts were also run on Windows 11 (x86-64) with Python 3.11.9,
NumPy 2.4.6, SciPy 1.17.1. Not yet run on macOS. Run the scripts from this directory.

## Files

| File | What it is |
|---|---|
| `scalarnd.py` | the implementation: problem, candidate, certified upper bound, fields, certified lower bound, the two ways of choosing a field, certificates on disk |
| `ia.py` | outward-rounded interval arithmetic on NumPy arrays |
| `run_examples.py` | the examples of `docs/SCALAR_ND.md`; writes `results/NAME.json` + `.npz` |
| `verify_all.py` | recomputes both bounds of every stored certificate from its exact data |
| `check_exact.py` | independent evaluator in exact rational arithmetic, compared with the interval one |
| `check_theorem.py` | falsification tests against problems with known minima (8 to 15 minutes) |
| `make_testvectors.py` | small inputs and expected outputs for the Wolfram Language port (`testvectors/`) |
| `make_table.py` | the results table of `docs/SCALAR_ND.md` from the stored certificates |
| `ref1d.py`, `refseries.py` | floating-point reference values (1D first integral, sine series); not rigorous |
| `results/` | stored certificates (`NAME.json`, `NAME.npz`), level-by-level logs of the calibration runs (`NAME_levels.txt`), `table.md` |
| `testvectors/` | inputs and expected outputs for `Tests/ScalarFieldND.wlt` (180 kB) |
| `logs/` | output of the check scripts as run for this version (Linux and Windows) |
| `PORTING.md` | the Wolfram Language port: function map, conventions, differences, port log, observations |

## Use

```python
from fractions import Fraction as Fr
import scalarnd as S

# convex energy: seconds
P = S.Problem(box=[(0, 1)] * 3, N=(32, 32, 32), k=Fr(1, 2),
              W=[0, -1, Fr(1, 2), 0, Fr(1, 4)],          # W(u) = u^4/4 + u^2/2 - u, ascending coefficients
              bc=[('D', 'D')] * 3, g=lambda x, y, z: 0 * x)
out, U, cal = S.certify_convex(P)                        # out['lower'] <= inf E <= out['upper']

# non-convex energy: minutes; coarse-to-fine list of (cells, knot intervals)
mk = lambda N: S.Problem(box=[(-1, 1)] * 2, N=N, k=Fr(1, 2), W=[4, 0, -8, 0, 4],
                         bc=[('D', 'D')] * 2, g=lambda x, y: x + 0 * y)
rows, U, cal, P = S.certify_calibration(mk, [((8, 8), 8), ((16, 16), 8), ((16, 16), 16), ((32, 32), 16)],
                                        trange=(-1, 1), U0=lambda x, y: x + 0 * y, scale=4)
```

`bc[i] = (type at the lower face, type at the upper face)` in direction i, `'D'` or `'N'`. Any nodal array can be
passed as the candidate (`S.upper_bound(P, U)`), and any coefficients as the field (`S.lower_bound(P, cal)`): both
functions are rigorous whatever they are given.

Run the checks (`python3` instead of `python` where that is the name):

```
python verify_all.py          # stored certificates, about 20 s
python check_exact.py         # about 10 s
python check_theorem.py       # 8 to 15 minutes
python run_examples.py fast   # the convex examples again, about 20 s (overwrites their certificates)
```

Each of the first three ends with a line in capitals saying that everything passed, and exits with status 0.
`check_theorem.py` prints SciPy integration warnings from the floating-point reference values; they are harmless.

## What is rigorous and what is not

Rigorous (outward-rounded arithmetic on exact inputs): `upper_bound`, `lower_bound`, `flux_lower`, `conj_upper`,
`conj_cells`, `tail_inf`, `verify_certificate`. Heuristic (plain floating point, never trusted): `solve_candidate`,
`fv_solve`, `affine_from_candidate`, `optimise_calibration`, `refine_t`, `prolong_x`, and the reference values in
`ref1d.py` and `refseries.py`.

The Dirichlet datum of the certified problem is the piecewise multilinear interpolant of the nodal values of `g` on
the boundary grid. It equals `g` when `g` is multilinear on each boundary cell face.

## Stored certificates

`results/NAME.json` holds the problem in exact form (rationals as strings), the t-range, K, and the two claims
`lower` and `upper` as exact decimal strings. `results/NAME.npz` holds the boundary data, the candidate's nodal
values and the field coefficients as doubles; they are taken as the exact binary numbers they are.

`verify_certificate` recomputes both bounds and accepts if they are at least as strong as the claims. A claim is the
computed bound weakened by one part in 10⁹ and rounded outward to ten significant digits (`claim_strings`). The
margin is there because one heuristic ingredient is not bitwise reproducible across platforms: the floating-point
root finder that places the fine part of the conjugate tables (`conj_upper`). Any placement gives a valid bound;
only the last digits can move. Measured here: perturbing the roots by a relative 10⁻⁶ moved the lower bounds of the
eleven examples by at most 5×10⁻¹³ relative.

Sizes: the eleven certificates take 14 MB, of which `ramp3d_32.npz` is 8.4 MB. The four largest data files
(`ramp2d.npz`, `saddle2d.npz`, `screened3d_48.npz`, `ramp3d_32.npz`; 12.8 MB together) are not in the repository.
Their claims (`NAME.json`) are, and `verify_all.py` reports a certificate whose `.npz` file is absent as skipped.
`python run_examples.py NAME` computes such a certificate again; it overwrites the stored claim, and the new bounds
differ from the old ones in their last digits, because the optimiser takes a different path on another platform.
