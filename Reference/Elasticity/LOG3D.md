# Log: linear elasticity in three dimensions (reference implementation)

Branch `elasticity-3d`. Note: `docs/ELASTICITY_3D.md`. Code: `Reference/Elasticity/elasticity3d.py`.
Entries are dated; the newest state is at the top.

## State (23 November 2026)

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
