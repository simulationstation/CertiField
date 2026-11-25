# Environment

Machine used for the benchmark pack and the audit (recorded 2026-11-24).

- OS: Microsoft Windows [Version 10.0.26200.9974]; 32 logical processors.
- Python 3.11.9 with NumPy 2.4.6, SciPy 1.17.1, mpmath 1.3.0, SymPy 1.14.0, matplotlib 3.11.1, PyTorch 2.13.0+cpu (global environment; the audit uses its own environment, versions in its summary).
- Wolfram Engine 15.0.0 for Microsoft Windows (64-bit), Windows-x86-64; two kernels can run at a time.
- Not installed at the start: scikit-fem, gmsh, meshio, python-flint, DOLFINx, CalculiX, conda, R, Julia.

## Base check of P3 (Reference/Curved3D)

- `python tests.py`: ALL TESTS PASSED (159 s).
- `python verify.py 24 clamped certificates/psi_clamped_N24.npz certificates/v_clamped_N24.npz 3 5`: VERIFIED 1.817059 <= K <= 1.834215 (336 s), as expected.
- `python verify.py 24 frictionless certificates/psi_frictionless_N24.npz certificates/v_frictionless_N24.npz 3 5`: VERIFIED 1.744878 <= K <= 1.758911 (332 s), as expected.
