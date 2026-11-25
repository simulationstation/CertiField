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
| [P1](problems/P1.md) | square plate with a circular hole, pulled by frictionless grips | stiffness K; mean of σ<sub>xx</sub> over a small patch at the edge of the hole | 0.5912009869 ≤ K ≤ 0.5912009882; 1.032737550 ≤ Q2 ≤ 1.032768058 |
| [P1c](problems/P1.md) | the same plate, clamped grips | the same | 0.5957196143 ≤ K ≤ 0.5957222828; 1.045047061 ≤ Q2 ≤ 1.046096446 |
| [P2](problems/P2.md) | L-bracket, clamped on one edge, the opposite grip moved | stiffness K; mean of σ<sub>xx</sub> over a square at the re-entrant corner | 0.1226044392 ≤ K ≤ 0.1226073846; −0.3095832144 ≤ Q2 ≤ −0.3093253171 |
| [P3](problems/P3.md) | cube with a spherical cavity, pulled by rigid grips | stiffness K | clamped: 1.832048486 ≤ K ≤ 1.832438844; frictionless: 1.757411214 ≤ K ≤ 1.757622167 |

All intervals are CERTIFIED and re-checked by the independent checker; the arguments have been checked by their
author only. Certificates are in `certificates/` (format in `certificates/FORMAT.md`), the method and derivations in
`NOTES.md`.

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

## Compare your own result

```
python benchpack/compare/compare.py P1 --Q1 0.59121 --Q2 1.03275
python benchpack/compare/compare.py --list
```

reports whether your values lie inside the certified intervals, and what a value outside usually means.
`results/VERDICTS.md` has scikit-fem, CalculiX and Wolfram NDSolve on every problem: with quadratic elements on curved
boundary edges all three converge into every interval; straight-edged meshes land above, by an amount that is almost
entirely the polygonal approximation of the hole (`harness/README.md` has the commands).

## Check a certificate yourself

```
cd benchpack/checker
python check.py ../certificates/P1_frictionless_N16_p9        # about a minute; P1c, P2 likewise
python tests.py
python p3check.py ../certificates/P3_frictionless_N12_p4_exact   # P3 (see checker/README.md)

cd Reference/Curved3D                                         # P3, from the repository root
python verify.py 24 frictionless certificates/psi_frictionless_N24.npz certificates/v_frictionless_N24.npz 3 5
python exact_check.py frictionless certificates/psi_frictionless_N24.npz --v certificates/v_frictionless_N24.npz
```

Python 3 with NumPy, SciPy and mpmath; nothing else.

## What is guaranteed and what is not

Guaranteed: the stated quantity of linear elasticity (small strains, the given material) for the exact body and grip
conditions lies in the interval. Not claimed: anything about real specimens (material laws, large deformation,
contact), or about quantities other than the ones defined. The arguments have been checked by their author only.

Code: MIT licence; data: CC BY 4.0 (see `LICENSE` in this folder).
