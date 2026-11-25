# Changelog

## v0.1 (November 2026)

Certificates (all CERTIFIED, re-checked by the independent checker; checked by the agent only):

- P1, frictionless grips: `P1_frictionless_N16_p9`.
- P1c, clamped grips: `P1c_N32_p7`.
- P2: `P2_h8_L8_p6`, and the finer `P2_h32_L7_p6` (Q2 to relative width 8.3e-4); both intervals of the finer
  certificate lie inside those of the first.
- P3: the grid-24 certificates of `Reference/Curved3D` (clamped and frictionless), unchanged; their exact check was
  rerun and passed.

Comparison harness: scikit-fem, CalculiX and Wolfram NDSolve on all problems (`results/VERDICTS.md`).

Errors found:

- None in a stored result.
- In the builder before any certificate was stored: a cut-cell rule that lost accuracy to cancellation (rewritten),
  and a solver rule that could drop the degrees of freedom of the particular field (now protected). Both could only
  have weakened a bound, never made it false.
