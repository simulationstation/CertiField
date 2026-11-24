# Skyrme fields of degree B: heuristic data for the certified programme

The certified lower bound for Skyrme fields of degree ±B is computed in the Wolfram Language
(`Kernel/SkyrmeLevelSetDegree.wl`); the statement and its proof are in `docs/SKYRME_DEGREE_TWO.md`. The two
programmes here are ordinary floating point. Nothing they compute enters a certified number except as a
heuristic choice.

| File | What it does |
|---|---|
| `relaxed.py` | Integrates the extremals of the relaxed level-set problem in degree B and finds the one that satisfies the end condition. Prints the value of the problem along it, in units of 12π². |
| `make_grids.py` | Writes the grid of the certified programme for a degree, a pion mass and a resolution into `Data/SkyrmeDegree`: levels, the volume of the extremal at each level, node offsets, and the cost of the extremal beyond each level. |

```
python relaxed.py                  # degree one (1.19210, 1.41439) and degree two (2.20844, 2.61111), mu = 0 and 1
python make_grids.py 2 0 quick     # Data/SkyrmeDegree/b2_mu0_quick
```

Requirements: Python 3 with NumPy and SciPy.

How the data are used.

* The nodes of level i are 0 and `vref[i] * eug[offs[i]:]`. The certified programme builds its piecewise affine
  lower bounds on these cells. Any increasing levels in (0, π) and any positive increasing nodes give a valid
  bound; these follow the extremal so that the bound is tight.
* `ctg.f64` holds the cost of the extremal beyond each level. The extremal is an admissible path of the relaxed
  problem, so this cost is an upper bound of the value function at the extremal's volume. The path test
  (`tools/skyrme_degree_claims.wls --write`, `Tests/SkyrmeDegree.wlt`) checks at every level that the certified
  function does not exceed it. This compares a certified number with a floating-point one and is a test, not
  part of a proof.

For degree one and μ = 0 the generator reproduces the grid `Data/Skyrme/mu0_quick` that ships with the paclet:
the same number of levels (1012) and node offsets, levels to 10⁻⁹ and reference volumes to 5×10⁻⁵.

For B ≥ 2 the extremal starts where the area constraint of the relaxed problem is not active: the level surface
has more area than a sphere of the same volume (B separate cores). `relaxed.py` therefore starts the integration
with a free multiplier instead of the central slope used for B = 1.
