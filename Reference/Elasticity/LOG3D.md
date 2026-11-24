# Log: linear elasticity in three dimensions (reference implementation)

Branch `elasticity-3d`. Note: `docs/ELASTICITY_3D.md`. Code: `Reference/Elasticity/elasticity3d.py`.
Entries are dated; the newest state is at the top.

## State when the work was interrupted (23 November 2026, evening)

**Done and committed**

* Note `docs/ELASTICITY_3D.md`: statements and proofs (written first), states, solver, checks, results, limits,
  references and what was searched. Section 3.5 is a sketch about tightness and says so.
* `elasticity3d.py`: complete for the Python reference (problem, exact local matrices, states, rigorous
  `upper_bound` and `lower_bound`, heuristic solvers, certificates, raw export).
* `check_exact3d.py`: 24 exact checks, all passed with the final evaluator (`logs3d/check_exact3d.txt`, 4 min).
* `check_theorem3d.py`: parts A (ten known minima) and B (falsification) passed with the final evaluator
  (`logs3d/check_theorem3d_partial.txt`).
* `run_examples3d.py` and seven stored certificates in `results3d/`, each re-verified when written
  (`logs3d/run_examples3d.txt`).
* Acceptance: E-T4 gap 0.859 % at 5×5×5 in about 10 s and 0.667 % at 6×6×6 in about 14 s. The stop rule does
  not apply.

**Half-done: exact state**

1. `check_theorem3d.py`, part C (number of states for all 63 arrangements of faces, with and without body force)
   was stopped while running. It had crashed once before on a missing attribute, which is fixed (`_all_states` is
   now set when there is no traction-free face), and it has not completed since. Run the script again (about 8
   minutes alone) and replace `logs3d/check_theorem3d_partial.txt`; if it passes, restore the statement about all
   63 arrangements in Sections 3.4 and 7 of the note.
2. `verify_all3d.py` and `make_testvectors3d.py` are written and have never been run. `testvectors3d/` does not
   exist yet. Run `python make_testvectors3d.py` (under a minute), check the output, add the directory.
3. The table in Section 8 of the note was assembled from runs made while other jobs were running; the times are
   upper estimates. An 8×8×8 run of E-T4 was planned (about 5000 iterations, 3 minutes) and not made.

**Next, in order**

1. Items 1 to 3 above.
2. `python verify_all3d.py`, and a CI step for the three scripts (the exact check takes 4 minutes; `--quick` runs
   four of the twelve problems).
3. Faster exact local matrix (integer arithmetic) and a better preconditioner or nested iteration for grids
   beyond 8×8×8; optional.
4. The Wolfram Language port. The states are reproducible with `RowReduce` (the echelon form is unique; columns
   numbered 64f + 16γ₁ + 4γ₂ + γ₃) and the exact resultant rule; `testvectors3d/*/states.json` will hold them.

## Earlier state (23 November 2026)


**Done**

* The note with the statements and proofs (Theorem 1, Lemma 2, Lemma 3, Corollary 4): written before the code.
* `elasticity3d.py`: problem class (box, cells, 6×6 rational stiffness, body force, Dirichlet or traction-free
  faces, Q_p datum), the Beltrami table generated from σ_ij = ε_ikl ε_jmn ∂_k∂_m Φ_ln, exact local matrices,
  the polynomial states for traction-free faces (exact sparse row reduction), rigorous `upper_bound` and
  `lower_bound`, heuristic solvers (dense eigen-decomposition for small grids, Jacobi-preconditioned conjugate
  gradients otherwise), certificates and raw export.
* First runs (floating-point solvers, rigorous evaluation):
  * known minima enclosed with gaps of 10⁻¹²: homogeneous strain with all faces Dirichlet (Beltrami, Maxwell and
    Morera functions; p = 1, 2); uniaxial tension and pure bending of a bar with four traction-free faces;
  * E-T4 (unit cube, λ = 3/2, μ = 1, clamped at x₁ = 0, displaced by (0.1, 0, 0) at x₁ = 1, four free faces,
    p = 2, six stress functions): gap 3.17 % (2³), 1.75 % (3³), 1.17 % (4³), 0.86 % (5³, about 20 s).
    So the acceptance criterion (below 1 % in under ten minutes) is met at 5³.
  * Maxwell's three functions give the same bound as all six for the symmetric pull, and a weaker one for
    transverse displacement and for twist of the end face (twist, 3³: 0.0007175 against 0.0007215).

**In progress**

* Which polynomial states to keep. With the rule "independent of the face data of the polynomial null sextuples"
  a ring of four free faces gets 18 states and two opposite free faces get 12, where the continuum argument
  expects 6 and 0. The extra ones make the dual system badly conditioned. New default rule: for a ring, six states
  with independent resultants (force and moment through the cross-section, exact); none otherwise. The comparison
  of the two rules on small grids is running.

**Next**

1. `check_exact3d.py`: exact rational cross-check on tiny grids (equilibrium, traction continuity, free faces,
   both exact values inside the enclosures, the hypercircle identity).
2. `check_theorem3d.py`: known minima with dyadic data (E-T1, E-T2) and falsification (E-T3).
3. `run_examples3d.py`, stored certificates in `results3d/`, `verify_all3d.py`, test vectors in `testvectors3d/`.
4. Sections 7 to 9 of the note (checks, results, what is not covered), references, what was searched.
5. Faster dual solve for grids beyond 6³ (nested iteration or a block preconditioner); optional.

## Observations so far

* The system for the stress field is singular (fields ½(∇w + ∇wᵀ) give no stress). Conjugate gradients converge
  when the regularisation is a multiple of the diagonal; with a multiple of the identity they stall, because the
  diagonal entries differ by ten orders of magnitude between the types of nodal values.
* Costs at 5³ with six functions: exact local matrix about 4 s, dual solve 14 s (10 700 iterations), rigorous
  evaluation 4 s.

2026-11-24: ported to the paclet (`Kernel/Elasticity3D.wl`, model `LinearElasticity` with a three-dimensional domain); test vectors agree to 10⁻¹², E-T4 cube 5³ certified with gap 0.86 % in 9 s; items 1 and 3 of the list above dropped.
