"""P1 on straight-edged meshes: how much of K - K_exact is geometry (polygonal hole) and how much discretisation?

    python benchpack/harness/p1_geometry_split.py

For a straight-edged mesh of level h the domain is a polygon Omega_h (vertices on the circle, so Omega_h contains
Omega).  Write
    K_h - K(Omega) = [K_h - K(Omega_h)] + [K(Omega_h) - K(Omega)]  =  discretisation + geometry.
K(Omega_h), the exact stiffness of the polygonal plate, is estimated by refining the same straight mesh uniformly
(new nodes on the chords, so the domain stays Omega_h) and solving with quadratic elements; the last two refinements
give the estimate (Richardson extrapolation when there are two or more) and its spread.  K(Omega) is taken as the
midpoint of the certified interval.  The same split is made for Q2 (average over the polygonal patch P_h of the mesh).
Refinements are made while the problem stays under the memory cap of runlog.py; for the finest level (no refinement
fits) the geometry part is estimated as c * (area of Omega_h - area of Omega) with c fitted on the coarser levels
(marked "estimated" in the note column).
Output: benchpack/results/p1_geometry_split.csv (resumable: a level with a row is not recomputed).
All values FLOAT; checked by the agent only.
"""
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runlog  # noqa: E402  (caps the threads before numpy is loaded)

import numpy as np  # noqa: E402
from skfem import CellBasis, ElementTriP1, ElementTriP2, ElementVector, Functional, MeshTri1, condense  # noqa: E402
from skfem.helpers import sym_grad, trace  # noqa: E402

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



FIELDS = ["h", "element", "ndof", "K_h", "K_polygon", "K_polygon_spread", "K_discretisation", "K_geometry", "K_total",
          "Q2_h", "Q2_polygon", "Q2_polygon_spread", "Q2_discretisation", "Q2_geometry", "Q2_total",
          "area_excess", "refinements", "note"]
LEVELS = [0.05, 0.025, 0.0125, 0.00625]


def area_excess(h):
    m = cases.make_mesh("P1", h, 1, False)
    P, C = m["points"], m["cells"]
    a, b, c = P[C[:, 0]], P[C[:, 1]], P[C[:, 2]]
    area = np.sum(0.5 * np.abs((b - a)[:, 0] * (c - a)[:, 1] - (b - a)[:, 1] * (c - a)[:, 0]))
    return float(area - (1 - math.pi / 16))


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    iv = IV.load()["P1"]
    Kref = float(iv["Q1"][0] + iv["Q1"][1]) / 2
    Qref = float(iv["Q2"][0] + iv["Q2"][1]) / 2
    log = runlog.CsvLog(os.path.join(PACK, "results", "p1_geometry_split.csv"), FIELDS, key=("element", "h"))
    for h in LEVELS:
        if log.has(element="P2-straight", h="%.12g" % h):
            continue
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
        # the polygonal domain Omega_h: refine the straight mesh (quadratic elements) while under the memory cap
        seq = [res["P2-straight"]]
        mesh = mesh0
        finder = mesh0.element_finder()
        while run_skfem.memory_gb(2, 4 * seq[-1][2], True) <= runlog.MAX_GB:
            mesh = mesh.refined()
            cen = mesh.p[:, mesh.t].mean(axis=1)
            patch = np.nonzero(tag0[finder(cen[0], cen[1])])[0]
            t0 = time.perf_counter()
            seq.append(solve(mesh, ElementTriP2(), patch))
            print("h=%g refined %d: K(Omega_h) ~ %.12f  Q2 ~ %.10f  ndof=%d (%.1fs)" % (
                h, len(seq) - 1, seq[-1][0], seq[-1][1], seq[-1][2], time.perf_counter() - t0), flush=True)
        ae = area_excess(h)
        note = ""
        if len(seq) >= 3:
            # Richardson from the last three values (geometric convergence)
            (k0, q0, _), (k1, q1, _), (k2, q2, _) = seq[-3:]
            r = (k2 - k1) / (k1 - k0) if k1 != k0 else 0.0
            KO = k2 + (k2 - k1) * r / (1 - r) if 0 < r < 1 else k2
            rq = (q2 - q1) / (q1 - q0) if q1 != q0 else 0.0
            QO = q2 + (q2 - q1) * rq / (1 - rq) if 0 < rq < 1 else q2
            dK, dQ = abs(KO - k2), abs(QO - q2)
        elif len(seq) == 2:
            KO, QO = seq[-1][0], seq[-1][1]
            dK, dQ = abs(seq[-1][0] - seq[-2][0]), abs(seq[-1][1] - seq[-2][1])
        else:
            # no refinement fits under the cap: geometry from the fitted factor of the coarser levels
            fits = [(float(r["K_geometry"]) / float(r["area_excess"]), float(r["Q2_geometry"]) / float(r["area_excess"]))
                    for r in log.rows if r["element"] == "P2-straight" and not r.get("note")]
            cK = float(np.mean([f[0] for f in fits]))
            cQ = float(np.mean([f[1] for f in fits]))
            KO, QO = Kref + cK * ae, Qref + cQ * ae
            dK = abs(np.std([f[0] for f in fits]) * ae)
            dQ = abs(np.std([f[1] for f in fits]) * ae)
            note = "estimated: geometry = c * area excess, c fitted on the coarser levels (no refinement fits the cap)"
        for variant in ("P1", "P2-straight"):
            Kh, Qh, n = res[variant]
            log.append(dict(h="%.12g" % h, element=variant, ndof=n, K_h="%.12e" % Kh, K_polygon="%.12e" % KO,
                            K_polygon_spread="%.2e" % dK, K_discretisation="%.4e" % (Kh - KO),
                            K_geometry="%.4e" % (KO - Kref), K_total="%.4e" % (Kh - Kref),
                            Q2_h="%.12e" % Qh, Q2_polygon="%.12e" % QO, Q2_polygon_spread="%.2e" % dQ,
                            Q2_discretisation="%.4e" % (Qh - QO), Q2_geometry="%.4e" % (QO - Qref),
                            Q2_total="%.4e" % (Qh - Qref), area_excess="%.4e" % ae, refinements=len(seq) - 1,
                            note=note))
            print("   %-12s K: total %+.3e = discretisation %+.3e + geometry %+.3e (spread %.1e)  %s" % (
                variant, Kh - Kref, Kh - KO, KO - Kref, dK, note), flush=True)
    print("reference: K(Omega) = %.12f, Q2 = %.10f (midpoints of the certified intervals)" % (Kref, Qref))


if __name__ == "__main__":
    main()
