# Certified benchmark pack for linear elasticity

Verification problems in solid mechanics whose benchmark quantities come with a **guaranteed two-sided interval**:
the exact value of the stated quantity, for the stated body, is proved to lie inside. Curved surfaces are
represented exactly (not meshed), every rounding error is enclosed, and every certificate is re-checked by a separate
program that shares nothing with the code that produced it except the problem definition.

It is meant for people who develop simulation codes and for engineers who have to show that a simulation can be
trusted: run your code on a problem, and the interval tells you whether your number is right, and by how much it can
still be wrong.

## Problems

| | Body | Quantities | Certified interval |
|---|---|---|---|
| [P1](problems/P1.md) | square plate with a circular hole, pulled by rigid grips (frictionless; clamped variant P1c) | stiffness K; mean of σ<sub>xx</sub> over a small patch at the edge of the hole | — |
| [P2](problems/P2.md) | L-bracket, clamped on one edge, the opposite grip moved | stiffness K; mean of σ<sub>xx</sub> over a square at the re-entrant corner | — |
| [P3](problems/P3.md) | cube with a spherical cavity, pulled by rigid grips | stiffness K | clamped: 1.81732246 ≤ K ≤ 1.83348613; frictionless: 1.74510549 ≤ K ≤ 1.75827651 |

The definitions in `problems/` are exact: geometry, material, grips, and the quantities to compute.

## How the intervals are obtained

- **Upper bound:** the energy of a displacement field that satisfies the grip conditions exactly. By the principle of
  minimum potential energy no displacement has lower energy than the true one, so its energy bounds K from above.
- **Lower bound:** the complementary energy of a stress field that is exactly in equilibrium and exactly free of
  traction on the free surfaces, including the curved ones (principle of Prager and Synge). The stress field is a
  simple equilibrated field plus stress-function terms multiplied by ρ², where ρ is a polynomial that vanishes on the
  curved surface, so the traction there vanishes as a polynomial identity.
- **Arithmetic:** Gauss rules that are exact for the polynomial degree, with nodes and weights enclosed in rational
  intervals; outward-rounded interval arithmetic for every operation; the body is covered by boxes decided in exact
  integer arithmetic.
- **Independent check:** a second program verifies the identities (equilibrium, zero traction, grip conditions, total
  force, coverage of the body) in exact rational arithmetic.

## Check P3 yourself

```
cd Reference/Curved3D
python tests.py                                                                                     # about 3 minutes
python verify.py 24 frictionless certificates/psi_frictionless_N24.npz certificates/v_frictionless_N24.npz 3 5
python exact_check.py frictionless certificates/psi_frictionless_N24.npz --v certificates/v_frictionless_N24.npz
```

Python 3 with NumPy and SciPy; nothing else.

## What is guaranteed and what is not

Guaranteed: the stated quantity of linear elasticity (small strains, the given material) for the exact body and grip
conditions lies in the interval. Not claimed: anything about real specimens (material laws, large deformation,
contact), or about quantities other than the ones defined. The arguments have been checked by their author only.

Code: MIT licence; data: CC BY 4.0 (see `LICENSE` in this folder).
