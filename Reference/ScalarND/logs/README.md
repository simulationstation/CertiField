# Logs and development history of the reference

## The logs

The first three were produced on 23 November 2026 with Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1 on Linux (x86-64),
one core. The `_windows` files were produced on the same day on Windows 11 (x86-64) with Python 3.11.9, NumPy 2.4.6,
SciPy 1.17.1, from the files of the handover archive, all eleven certificates present.

| File | Command | Result |
|---|---|---|
| `log_verify_all.txt` | `python verify_all.py` | all 11 stored certificates verified |
| `log_check_exact.txt` | `python check_exact.py` | all exact checks passed |
| `log_check_theorem.txt` | `python check_theorem.py` | all tests passed; largest (lower bound − reference) = −1.087×10⁻⁴ |
| `log_verify_all_windows.txt` | `python verify_all.py` | all 11 stored certificates verified |
| `log_check_exact_windows.txt` | `python check_exact.py` | all exact checks passed |
| `log_check_theorem_windows.txt` | `python check_theorem.py` | all tests passed; largest (lower bound − reference) = −1.087×10⁻⁴ |

The two platforms agree to the last printed digit in every rigorous number that does not depend on an optimiser
(the convex rows of `check_exact.py` and of `check_theorem.py`, including the closest approach of −1.087×10⁻⁴), and
`verify_all.py` accepts the same eleven claims. Rows that start from a field found by L-BFGS (the calibration rows
of `check_exact.py`, most of `check_theorem.py`) differ from the fifth or sixth digit on, because the optimiser takes
a different path with another SciPy; they pass on both.

The Linux logs were produced with the files as shipped, the first two also from a fresh unpacking of the handover
archive under an ASCII-only locale. The integration warnings in `log_check_theorem.txt` come from SciPy's quadrature
in the floating-point reference values and are harmless (the directory name in their file paths has been removed;
nothing else in the logs was edited).

## What was tried, what failed, what was fixed

Kept because a port can repeat the same mistakes. Where a number here differs from `docs/SCALAR_ND.md` or
`results/`, those are the current ones.

1. **Clamped splines only (first version).** Enough for the domain wall. A convex test showed two things: the
   tails outside the t-range bind unless the range covers the minimum of W, and with few samples in t the optimiser
   hides peaks of Φ between the samples. The second is why the curvature allowance of Lemma 3 is part of the
   optimiser's objective and not only of the verifier.
2. **Affine part added** (t σ₀, the classical convex dual; tails by conjugate tables). Optimising the affine and
   spline parts together converges very slowly, because the two overlap inside the t-range. Freezing the affine
   part at the finite-volume solution is bad for non-convex W: the tails bind and the bound collapses. Hence two
   modes: convex (affine part alone) and calibration (spline part alone).
3. **Convex mode.** Fluxes built from cell means of the candidate gave first-order accuracy only; solving the
   finite-volume problem gives second order. Integrating |σ₀|² exactly over each cell for K = 0 (Proposition 6)
   instead of using vertex values made the deficit four times smaller.
4. **Two defects in the rigorous evaluation, both found by tests and fixed.** Both made the bound weaker, not
   wrong.
   * Boundary flux: the Lipschitz penalty used the range of the datum over the whole face instead of sub-faces
     (a loss of 0.19 on the ramp example). Now 8 sub-faces per direction.
   * `_spline_eval` at t = t_hi: after rounding, the last knot interval was hulled with a successor that does not
     exist, which brought in the value at its other end. The condition `kap < K - 1` fixes it. The Wolfram Language
     `splineEval` needs the same care.
5. **Rounding.** `math.fsum` is widened by two units in the last place; `_spline_eval` asserts that its argument
   intervals are thin (less than half a knot interval), so that at most two pieces are visited; `tail_inf` is
   cached.
6. **Stored claims** carry a margin of one part in 10⁹ (see `../README.md`).
