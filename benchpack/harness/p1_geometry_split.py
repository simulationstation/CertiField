"""P1 on straight-edged meshes: how much of K - K_exact is geometry (polygonal hole) and how much discretisation?

    python benchpack/harness/p1_geometry_split.py

For a straight-edged mesh of level h the domain is a polygon Omega_h (vertices on the circle, so Omega_h contains
Omega).  Write
    K_h - K(Omega) = [K_h - K(Omega_h)] + [K(Omega_h) - K(Omega)]  =  discretisation + geometry.
K(Omega_h), the exact stiffness of the polygonal plate, is estimated by refining the same straight mesh uniformly
(new nodes on the chords, so the domain stays Omega_h) and solving with quadratic elements; the last two refinements
give the estimate and its spread.  K(Omega) is taken as the midpoint of the certified interval.  The same split is
made for Q2 (average over the polygonal patch P_h of the mesh).
Output: benchpack/results/p1_geometry_split.csv.  All values FLOAT; checked by the agent only.
"""
import csv
import os
import sys
import time

import numpy as np
from skfem import CellBasis, ElementTriP1, ElementTriP2, ElementVector, Functional, MeshTri1, condense
from skfem.helpers import sym_grad, trace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases  # noqa: E402
import run_skfem  # noqa: E402

PACK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PACK, "compare"))
import intervals as IV  # noqa: E402


def solve(mesh, elem, patch):
    """K and Q2 for P1 (frictionless) on the given straight mesh; patch = element indices of the patch."""
    c = cases.CASES["P1"]
    e = ElementVector(elem)
    basis = CellBasis(mesh, e)
    A = run_skfem.stiffness(basis, 2)
    x = np.zeros(basis.N)
    D = []
    for (axis, val), comp, g in c["bc"]:
        f = mesh.facets_satisfying(lambda X, axis=axis, val=val: np.abs(X[axis] - val) < 1e-9, boundaries_only=True)
        d = basis.get_dofs(f).all("u^%d" % (comp + 1))
        x[d] = g
        D.append(d)
    D = np.unique(np.concatenate(D))
    AII, bI, xx, I = condense(A, x=x, D=D)
    xx[I] = run_skfem.spsolve(AII, bI)
    K = c["energy_scale"] * float(xx @ (A @ xx))
    lam, mu = run_skfem.lame(2)
    pb = CellBasis(mesh, e, elements=patch, intorder=6)

    @Functional
    def sxx(w):
        eps = sym_grad(w["u"])
        return 2 * mu * eps[0, 0] + lam * trace(eps)

    @Functional
    def one(w):
        return 1.0 + 0.0 * w.x[0]
    return K, sxx.assemble(pb, u=pb.interpolate(xx)) / one.assemble(pb), basis.N


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    iv = IV.load()["P1"]
    Kref = float(iv["Q1"][0] + iv["Q1"][1]) / 2
    Qref = float(iv["Q2"][0] + iv["Q2"][1]) / 2
    rows = []
    plan = {0.05: 3, 0.025: 2, 0.0125: 2, 0.00625: 1}
    for h, nref in plan.items():
        m = cases.make_mesh("P1", h, 1, False)
        mesh0 = MeshTri1(m["points"].T.copy(), m["cells"].T.copy())
        patch0 = np.nonzero(cases.in_p1_patch(mesh0.p[:, mesh0.t].mean(axis=1)))[0]
        tag0 = np.zeros(mesh0.nelements, bool)
        tag0[patch0] = True
        res = {}
        for variant, elem in (("P1", ElementTriP1()), ("P2-straight", ElementTriP2())):
            t0 = time.perf_counter()
            res[variant] = solve(mesh0, elem, patch0)
            print("h=%g %-12s K=%.12f Q2=%.10f ndof=%d (%.1fs)" % (h, variant, *res[variant],
                                                                   time.perf_counter() - t0), flush=True)
        # the polygonal domain Omega_h: refine the straight mesh, quadratic elements
        seq = []
        mesh = mesh0
        for r in range(1, nref + 1):
            mesh = mesh.refined()
            finder = mesh0.element_finder()
            cen = mesh.p[:, mesh.t].mean(axis=1)
            parent = finder(cen[0], cen[1])
            patch = np.nonzero(tag0[parent])[0]
            t0 = time.perf_counter()
            Kr, Qr, n = solve(mesh, ElementTriP2(), patch)
            seq.append((Kr, Qr, n))
            print("h=%g refined %d: K(Omega_h) ~ %.12f  Q2 ~ %.10f  ndof=%d (%.1fs)" % (
                h, r, Kr, Qr, n, time.perf_counter() - t0), flush=True)
        KO, QO = seq[-1][0], seq[-1][1]
        prev = res["P2-straight"] if len(seq) == 1 else seq[-2]
        dK, dQ = abs(KO - prev[0]), abs(QO - prev[1])
        for variant in ("P1", "P2-straight"):
            Kh, Qh, n = res[variant]
            rows.append(dict(h=h, element=variant, ndof=n, K_h=Kh, K_polygon=KO, K_polygon_spread=dK,
                             K_discretisation=Kh - KO, K_geometry=KO - Kref, K_total=Kh - Kref,
                             Q2_h=Qh, Q2_polygon=QO, Q2_polygon_spread=dQ, Q2_discretisation=Qh - QO,
                             Q2_geometry=QO - Qref, Q2_total=Qh - Qref))
            print("   %-12s K: total %+.3e = discretisation %+.3e + geometry %+.3e (spread %.1e)" % (
                variant, Kh - Kref, Kh - KO, KO - Kref, dK), flush=True)
    path = os.path.join(PACK, "results", "p1_geometry_split.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: ("%.6e" % v if isinstance(v, float) and k not in ("h",) else v) for k, v in r.items()})
    print("reference: K(Omega) = %.12f, Q2 = %.10f (midpoints of the certified intervals)" % (Kref, Qref))


if __name__ == "__main__":
    main()
