# Comparison harness

Runs open-source finite element codes on P1, P1c, P2 and P3 (definitions in `../problems/`) and sets their values of
Q1 = K and Q2 against the certified intervals. Results: `../results/`; verdicts: `../results/VERDICTS.md`.
All finite element values are FLOAT; checked by the agent only.

## Install

    pip install scikit-fem gmsh meshio pypardiso matplotlib scipy     # pypardiso optional (faster direct solver)
    conda create -n ccx -c conda-forge calculix                         # CalculiX 2.24 (any conda; Miniforge used here)

`run_calculix.py` finds `ccx` on PATH or in `~/miniforge3/envs/ccx` (set `CCX_PREFIX` for another env).
NDSolve needs a Wolfram kernel (`wolframscript`; set `WOLFRAMSCRIPT` to a wrapper command if needed).

## Run (from the repository root)

    python benchpack/harness/run_skfem.py               # all cases, or e.g.  P1 P2 P3c
    python benchpack/harness/run_calculix.py
    python benchpack/harness/run_ndsolve.py             # NDSolve on the same gmsh meshes
    wolframscript -file benchpack/harness/ndsolve_P3.wls   # NDSolve with its own mesher (copy of Reference/Curved3D/ndsolve_compare.wls)
    python benchpack/harness/geometry_check.py          # where the meshed hole / cavity lies relative to the exact one
    python benchpack/harness/summarize.py               # plots results/<case>.png and results/verdicts_table.md
    wolframscript -file benchpack/harness/repro_ndsolve_gradient.wls   # minimal NIntegrate reproducer

Cases (`cases.py`): `P1solid` (P1 without the hole: K = 1, Q2 = 1/2 exactly), `P1` (frictionless grips), `P1c`
(clamped grips), `P2`, `P3c` / `P3f` (P3 with clamped / frictionless grips). P1 and P3 use the quarter / octant model
with the mirror conditions of the problem files; P2 the full model.

## Files

- `meshes.py` gmsh geometry and meshes (the Q2 patch is a separate surface, so meshes conform to it; graded towards the
  re-entrant corner and the clamped-grip corners of P2, P1c and P3c);
- `cases.py` grip conditions, scale factors from the model to the full body, mesh levels;
- `run_skfem.py`, `run_calculix.py` + `ccx_io.py`, `run_ndsolve.py` + `ndsolve_fem.wls`, `ndsolve_P3.wls`: the runs;
- `geometry_check.py`, `summarize.py`, `repro_ndsolve_gradient.wls`: checks, plots, verdict table, reproducer.
