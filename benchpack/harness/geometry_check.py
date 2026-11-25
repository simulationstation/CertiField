"""Where does the meshed boundary of the hole (P1) or cavity (P3) lie relative to the exact circle / sphere?

    python benchpack/harness/geometry_check.py

For every mesh level and variant it samples the boundary edges (P1) or faces (P3) of the mesh that lie on the hole
or cavity and reports min and max of r - 1/2 (r = distance to the centre), and the area / volume of the meshed model
minus the exact value.  If r <= 1/2 everywhere on the meshed hole boundary, the meshed body contains the exact body;
then, for these displacement-loaded problems (grips with prescribed displacement, every other boundary free), the
minimal energy over the meshed body is at least that over the exact body, and a conforming finite element value of K
on the meshed body is an upper bound for the exact K (FLOAT statement about the computed numbers; the inequality
itself is the standard monotonicity of the minimum energy under enlarging a body whose added part is traction-free).
Output: benchpack/results/geometry_check.csv.  Checked by the agent only.
"""
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases  # noqa: E402

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
EXACT = {"P1": 1 - math.pi / 16, "P3f": 1 - math.pi / 48}


def tri_edges_on_circle(P, C, order):
    """boundary edges with both end vertices on r = 1/2: list of (a, b, mid or None)"""
    edges = {}
    loc = [(0, 1, 3), (1, 2, 4), (2, 0, 5)]
    for c in C:
        for i, j, k in loc:
            key = tuple(sorted((c[i], c[j])))
            edges.setdefault(key, []).append(c[k] if order == 2 else None)
    out = []
    for (a, b), mids in edges.items():
        if len(mids) == 1 and abs(np.hypot(*P[a]) - 0.5) < 1e-9 and abs(np.hypot(*P[b]) - 0.5) < 1e-9:
            out.append((a, b, mids[0]))
    return out


def check_2d(m):
    P, C, order = m["points"], m["cells"], m["order"]
    s = np.linspace(0, 1, 41)
    rmin, rmax = np.inf, -np.inf
    for a, b, k in tri_edges_on_circle(P, C, order):
        if order == 1:
            X = np.outer(1 - s, P[a]) + np.outer(s, P[b])
        else:
            N = np.array([(1 - s) * (1 - 2 * s), 4 * s * (1 - s), s * (2 * s - 1)])
            X = np.outer(N[0], P[a]) + np.outer(N[1], P[k]) + np.outer(N[2], P[b])
        r = np.hypot(X[:, 0], X[:, 1]) - 0.5
        rmin, rmax = min(rmin, r.min()), max(rmax, r.max())
    return rmin, rmax


def check_3d(m):
    P, C, order = m["points"], m["cells"], m["order"]
    faces = {}
    # local faces of a tet (vertices) and their edge-node indices in VTK order (edges 01 12 20 03 13 23 -> 4..9)
    fl = [((0, 1, 2), (4, 5, 6)), ((0, 1, 3), (4, 8, 7)), ((1, 2, 3), (5, 9, 8)), ((0, 2, 3), (6, 9, 7))]
    for c in C:
        for vs, es in fl:
            key = tuple(sorted(c[list(vs)]))
            faces.setdefault(key, []).append((c[list(vs)], c[list(es)] if order == 2 else None))
    # sample points in a triangle: barycentric grid
    n = 8
    bary = [(i / n, j / n, 1 - i / n - j / n) for i in range(n + 1) for j in range(n + 1 - i)]
    rmin, rmax = np.inf, -np.inf
    for key, lst in faces.items():
        if len(lst) != 1:
            continue
        vs, es = lst[0]
        if np.any(np.abs(np.linalg.norm(P[vs], axis=1) - 0.5) > 1e-9):
            continue
        for l1, l2, l3 in bary:
            if order == 1:
                X = l1 * P[vs[0]] + l2 * P[vs[1]] + l3 * P[vs[2]]
            else:
                # quadratic triangle: vertices v0 v1 v2, edges (0,1) (1,2) (2,0) = es[0] es[1] es[2]
                X = (l1 * (2 * l1 - 1) * P[vs[0]] + l2 * (2 * l2 - 1) * P[vs[1]] + l3 * (2 * l3 - 1) * P[vs[2]]
                     + 4 * l1 * l2 * P[es[0]] + 4 * l2 * l3 * P[es[1]] + 4 * l3 * l1 * P[es[2]])
            r = np.linalg.norm(X) - 0.5
            rmin, rmax = min(rmin, r), max(rmax, r)
    return rmin, rmax


def measure(case, h, variant):
    """area / volume of the meshed model domain (scikit-fem quadrature of 1, isoparametric for curved meshes)"""
    import run_skfem
    from skfem import CellBasis, Functional
    mesh, e = run_skfem.build(case, h, variant)
    b = CellBasis(mesh, e, intorder=6)

    @Functional
    def one(w):
        return 1.0 + 0.0 * w.x[0]
    return one.assemble(b)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    rows = []
    for case in ("P1", "P3f"):
        for variant, order, curved in [("P1", 1, False), ("P2-straight", 2, False), ("P2-curved", 2, True)]:
            levels = cases.LEVELS[case] if case == "P1" else cases.LEVELS[case][:4]
            for h in levels:
                m = cases.make_mesh(case, h, order, curved)
                rmin, rmax = (check_2d if case == "P1" else check_3d)(m)
                meas = measure(case, h, variant)
                rows.append(dict(case=case, element=variant, h=h, r_minus_half_min="%.3e" % rmin,
                                 r_minus_half_max="%.3e" % rmax, measure_minus_exact="%.3e" % (meas - EXACT[case])))
                print(rows[-1], flush=True)
    with open(os.path.join(RESULTS, "geometry_check.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
