"""Run CalculiX (ccx) on the benchmark cases and write benchpack/results/calculix_<case>.csv.

    python benchpack/harness/run_calculix.py P1 [P2 P3c ...]   (all cases if none given)

Element variants (same gmsh meshes as run_skfem.py):
    2D: CPS3; CPS6-straight (mid-edge nodes at the chord midpoints); CPS6-curved (mid-edge nodes on the circle)
    3D: C3D4; C3D10-straight; C3D10-curved
CalculiX expands plane-stress elements into 3D elements of the section thickness (CPS3 -> C3D6, CPS6 -> C3D15).
With thickness 1 that is a 3D slab, not plane stress (sigma_zz != 0 near the hole; K differs by about 0.7 %).
The harness therefore uses a thin section, THICKNESS = 1e-3, and divides the grip force by it
(K changed by less than 1e-6 relative between 1e-3 and 1e-4).

Q1: total reaction on the grip node set (*NODE PRINT, TOTALS=ONLY), 7 significant digits in the .dat file.
    K_energy: total strain energy (*EL PRINT ELSE, TOTALS=ONLY) as a cross-check.
Q2: the integration-point stresses of the patch elements (*EL PRINT S, COORD).  The in-plane rule of each element is
    recovered from the printed integration-point coordinates (inverse isoparametric map), the weights times the
    Jacobian are computed from the element nodes, and the through-thickness points are averaged with their Gauss
    weights.  So Q2 is the code's own quadrature of the integral of sigma_xx over the patch, divided by the patch
    area under the same quadrature.
All values are FLOAT; checked by the agent only.
"""
import csv
import os
import shutil
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases   # noqa: E402
import ccx_io  # noqa: E402

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
THICKNESS = 1e-3
VARIANTS_2D = [("CPS3", 1, False), ("CPS6-straight", 2, False), ("CPS6-curved", 2, True)]
VARIANTS_3D = [("C3D4", 1, False), ("C3D10-straight", 2, False), ("C3D10-curved", 2, True)]
NTHREADS = int(os.environ.get("CCX_THREADS", "4"))
# ccx (SPOOLES solver, plane-stress expansion to 3D) stopped without a message at about 1e6 unknowns on the
# harness machine; larger runs are skipped
MAX_DOF = int(os.environ.get("CCX_MAX_DOF", "700000"))

# triangle rules: reference points (xi, eta) and weights (reference area 1/2)
TRI_RULES = {
    1: ([(1 / 3, 1 / 3)], [0.5]),
    3: ([(1 / 6, 1 / 6), (2 / 3, 1 / 6), (1 / 6, 2 / 3)], [1 / 6] * 3),
}
LINE_GAUSS = {1: [1.0], 2: [0.5, 0.5], 3: [5 / 18, 8 / 18, 5 / 18]}


def tri_shape(xi, eta, order):
    l1, l2, l3 = 1 - xi - eta, xi, eta
    if order == 1:
        N = np.array([l1, l2, l3])
        dN = np.array([[-1, -1], [1, 0], [0, 1]], dtype=float)
        return N, dN
    N = np.array([l1 * (2 * l1 - 1), l2 * (2 * l2 - 1), l3 * (2 * l3 - 1), 4 * l1 * l2, 4 * l2 * l3, 4 * l3 * l1])
    dl = np.array([[-1, -1], [1, 0], [0, 1]], dtype=float)
    dN = np.array([(4 * l1 - 1) * dl[0], (4 * l2 - 1) * dl[1], (4 * l3 - 1) * dl[2],
                   4 * (l1 * dl[1] + l2 * dl[0]), 4 * (l2 * dl[2] + l3 * dl[1]), 4 * (l3 * dl[0] + l1 * dl[2])])
    return N, dN


def invert(X, xy, order):
    """Reference coordinates of the point xy in the triangle with nodes X (Newton)."""
    xi = np.array([1 / 3, 1 / 3])
    for _ in range(30):
        N, dN = tri_shape(xi[0], xi[1], order)
        r = N @ X - xy
        J = X.T @ dN
        step = np.linalg.solve(J, r)
        xi = xi - step
        if np.linalg.norm(step) < 1e-14:
            break
    return xi


def patch_average(dat, points, cells, patch, order):
    S, CO = dat["S"], dat["COORD"]
    num = den = 0.0
    for e in patch:
        s = np.array(S[e + 1])
        co = np.array(CO[e + 1])
        X = points[cells[e]]
        nip = len(s)
        nz = {2: 2, 9: 3, 6: 2, 3: 1, 1: 1}[nip]          # C3D6: 1x2, C3D15: 3x3
        npl = nip // nz
        zs = co[:, 3]
        # ccx order (checked below): point k lies in layer k // npl at in-plane point j = k % npl of the rule
        pts, wts = TRI_RULES[npl]
        size = np.ptp(X[:, 0]) + np.ptp(X[:, 1])
        for j in range(npl):
            ks = [k for k in range(nip) if k % npl == j]
            if np.ptp(co[ks, 1]) + np.ptp(co[ks, 2]) > 1e-5 * max(size, 1e-3) + 1e-7 or \
                    np.any(np.diff(zs[ks]) <= 0):
                raise RuntimeError("unexpected integration-point layout in element %d" % (e + 1))
            sxx = sum(w * s[k, 1] for w, k in zip(LINE_GAUSS[len(ks)], ks))
            if size > 1e-3:       # the printed coordinates have 7 digits: check the rule on elements not too small
                xi = invert(X, co[ks, 1:3].mean(axis=0), order)
                if np.hypot(xi[0] - pts[j][0], xi[1] - pts[j][1]) > 1e-3:
                    raise RuntimeError("integration point does not match the %d-point rule: xi=%s" % (npl, xi))
            _, dN = tri_shape(pts[j][0], pts[j][1], order)
            detJ = abs(np.linalg.det(X.T @ dN))
            num += wts[j] * detJ * sxx
            den += wts[j] * detJ
    return num / den, den


def solve(case, h, variant, order, curved, workdir):
    c = cases.CASES[case]
    dim = c["dim"]
    t0 = time.perf_counter()
    m = cases.make_mesh(case, h, order, curved)
    t_mesh = time.perf_counter() - t0
    P, C = m["points"], m["cells"]
    if dim * len(P) > MAX_DOF:
        return None
    on = lambda ax, v: np.nonzero(np.abs(P[:, ax] - v) < 1e-9)[0]
    bcs = [(on(ax, val), comp + 1, g) for (ax, val), comp, g in c["bc"]]
    (gax, gval), gcomp = c["grip"]
    patch = None
    if c["q2"]:
        cent = P[C[:, :3]].mean(axis=1).T
        sel = cases.in_p1_patch(cent) if c["q2"] == "p1" else cases.in_p2_patch(cent)
        patch = np.nonzero(sel)[0]
    eltype = variant.split("-")[0]
    job = "%s_%s_%g" % (case, variant, h)
    inp = os.path.join(workdir, job + ".inp")
    ccx_io.write_inp(inp, P, C, eltype, bcs, on(gax, gval), patch, thickness=THICKNESS, el_vars=("S", "COORD"),
                     energy=True)
    t0 = time.perf_counter()
    ccx_io.run(workdir, job, nthreads=NTHREADS)
    t_solve = time.perf_counter() - t0
    dat = ccx_io.parse_dat(os.path.join(workdir, job + ".dat"))
    thick = THICKNESS if dim == 2 else 1.0
    F = dat["RF"][gcomp]
    K = c["kscale"] * F / thick * np.sign(c["d_model"])
    Ken = c["energy_scale"] * 2 * dat["ELSE"] / thick
    q2, area = float("nan"), float("nan")
    if patch is not None:
        q2, area = patch_average(dat, P, C, patch, order)
        area /= 1.0
    ndof = dim * len(P)
    for ext in (".inp", ".dat", ".frd", ".sta", ".cvg", ".log", ".12d"):
        try:
            os.remove(os.path.join(workdir, job + ext))
        except OSError:
            pass
    return dict(code="CalculiX", case=case, element=variant, h=h, ndof=ndof, nelem=len(C),
                time_s=round(t_solve, 3), mesh_time_s=round(t_mesh, 3), K=K, K_energy=Ken, Q2=q2, patch_area=area)


FIELDS = ["code", "case", "element", "h", "ndof", "nelem", "time_s", "mesh_time_s", "K", "K_energy", "Q2",
          "patch_area"]


def main(argv):
    sys.stdout.reconfigure(encoding="utf-8")
    todo = argv or list(cases.CASES)
    os.makedirs(RESULTS, exist_ok=True)
    workdir = tempfile.mkdtemp(prefix="ccx_bench_")
    for case in todo:
        c = cases.CASES[case]
        variants = VARIANTS_2D if c["dim"] == 2 else VARIANTS_3D
        rows = []
        path = os.path.join(RESULTS, "calculix_%s.csv" % case)
        for variant, order, curved in variants:
            if curved and c["mesh"] == "lbracket":
                continue                   # no curved boundary: identical to the straight variant
            levels = cases.LEVELS[case] if order == 1 else cases.LEVELS_QUAD.get(case, cases.LEVELS[case])
            for h in levels:
                row = solve(case, h, variant, order, curved, workdir)
                if row is None:
                    print("%-8s %-15s h=%-8g skipped (more than %d unknowns)" % (case, variant, h, MAX_DOF))
                    continue
                rows.append(row)
                print("%-8s %-15s h=%-8g ndof=%-8d t=%7.2fs  K=%.8f  (energy %.8f)  Q2=%.7f" % (
                    case, variant, h, row["ndof"], row["time_s"], row["K"], row["K_energy"], row["Q2"]), flush=True)
                with open(path, "w", newline="") as f:
                    wr = csv.DictWriter(f, FIELDS)
                    wr.writeheader()
                    for rr in rows:
                        wr.writerow({k: (("%.12g" % v) if isinstance(v, float) else v) for k, v in rr.items()})
    shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main(sys.argv[1:])
