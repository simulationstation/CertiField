# Inventory: what the benchmark pack reuses

| Item | Where | What it gives the pack | State |
|---|---|---|---|
| curved3d package (P3) | `Reference/Curved3D/` (`curved3d.py`, `ia.py`, `verify.py`, `exact_check.py`, `tests.py`, `smoke_test.py`, `certificates/`, `logs/`, `NOTES.md`) | P3 as it stands: block with a spherical cavity, clamped and frictionless grips; CERTIFIED grid-24 intervals (clamped 1.81732246 <= K <= 1.83348613, frictionless 1.74510549 <= K <= 1.75827651); the method of 3.1 (step-column particular field, Beltrami of rho^2 Psi, octree superset integration, exact Gauss rules, interval evaluation) and a separate exact-rational checker | Python; base check on this machine in `ENVIRONMENT.md` |
| 2D linear elasticity reference | `Reference/Elasticity/` (`elasticity.py`, `check_exact.py`, `check_theorem.py`, certificates, test vectors) | Airy stress functions in a C^1 Hermite space with exact traction-free faces, Q_p displacements, exact local matrices, claim strings, known-answer and falsification tests; rectangles only (no curved boundary) | Python, tested |
| 3D linear elasticity reference | `Reference/Elasticity/*3d*` | Beltrami stress functions on boxes, exact traction-free states | Python, tested |
| LinearElasticity (paclet) | `Kernel/Elasticity.wl`, `Kernel/Elasticity3D.wl` | 2D and 3D boxes in the Wolfram Language, with VerifyCertificate | not needed by the Python pack |
| curved3d port (paclet) | branch `curved3d`, `Kernel/Curved3D.wl`, model `ElasticCavity3D` | P3 in the Wolfram Language | archived, not merged (`archive/STATE_curved3d_wl_port.md`) |
| NDSolve comparison | `Reference/Curved3D/ndsolve_compare.wls` | P3 against Wolfram's FEM | written, not run |

Interval and exact arithmetic available here:

- `ia.py` (three copies: `Reference/ScalarND`, `Reference/Elasticity`, `Reference/Curved3D`): NumPy float64 intervals widened by one ulp per operation with `nextafter`; exact `Fraction` enclosures. Used by every Python certificate so far.
- `fractions.Fraction` (Python): exact rational arithmetic for the checkers.
- `mpmath.iv` 1.3.0: arbitrary-precision intervals (scalar).
- python-flint (Arb balls): not installed globally; the audit installs it in its own environment.
- Wolfram Engine 15.0: `Interval`, exact rationals, and the paclet's `RigorousArithmetic.wl` (two kernel seats on this machine).

Not present: scikit-fem, gmsh, meshio (pip), DOLFINx, CalculiX, conda, R, Julia.
