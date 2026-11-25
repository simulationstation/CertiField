# Decisions

Decisions taken under prompt_benchpack_and_audit.md, with the reason.

- 2026-11-24. Day 0, stop and archive: no worker was running (the previous session had ended); the unpushed work of the stopped tracks was committed and pushed to their own branches, nothing merged or deleted, and archive/STATE_*.md written. No Jakobson-Naud or fractal-uncertainty track exists in this repository or this session, so there was nothing to stop or archive for it.
- 2026-11-24. `outreach/` and `audit/drafts/` are git-ignored: the repository is public, and the prompt forbids outside actions; drafts stay on this machine for the owner.
- 2026-11-24. The pack is built in Python (reusing Reference/Curved3D and Reference/Elasticity) in new folders only; the paclet and its frozen Skyrme code are not touched.
- 2026-11-24. Harness, second code: CalculiX 2.24 (conda-forge `calculix`, win-64, in a user-level Miniforge at ~/miniforge3, env `ccx`). It installed in about five minutes and ran at once; DOLFINx would in addition need its JIT compiler chain on Windows (MSVC through cffi), an extra point of failure, so it was not tried. CalculiX expands plane-stress elements into 3D elements of the section thickness, so the harness uses a thin section (1e-3) and divides the force by it; with thickness 1 the model is a 3D slab, not plane stress.
- 2026-11-24. Harness, NDSolve for P3: a copy of Reference/Curved3D/ndsolve_compare.wls runs as benchpack/harness/ndsolve_P3.wls (the original lies outside benchpack/ and stays untouched); body, material, grips and the stiffness formula are unchanged, only CSV output and the grid-24 ranges were added.
