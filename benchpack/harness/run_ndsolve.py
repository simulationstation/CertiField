"""Run Wolfram NDSolve (finite elements) on the harness meshes and write benchpack/results/ndsolve_<case>.csv.

    python benchpack/harness/run_ndsolve.py P1 [P2 P3c ...]

The kernel is started as  $WOLFRAMSCRIPT -file benchpack/harness/ndsolve_fem.wls <jobs.wl>  (default: wolframscript
on PATH; set WOLFRAMSCRIPT to a wrapper command, e.g. '"C:/Program Files/Git/bin/bash.exe" /path/wl.sh', where the
kernel is shared).  If the kernel stops, it is started again for the jobs not yet done (at most 3 times).
Meshes are the gmsh meshes of run_skfem.py / run_calculix.py; element variants:
    linear (order-1 mesh), quadratic-straight (order 2, mid-edge nodes at the chord midpoints),
    quadratic-curved (order 2, mid-edge nodes on the circle / sphere).
Columns: K = K_system of ndsolve_fem.wls (a(u_h, u_h) from NDSolve's assembled system);
K_nintegrate, Q2_nintegrate = NIntegrate of the InterpolatingFunction (the usual NDSolve post-processing);
Q2 and K_check = the quantities of NDSolve's nodal displacement field, evaluated here with the element shape
functions (scikit-fem basis on the same mesh, same element space).
All values are FLOAT; checked by the agent only.
"""
import csv
import os
import shlex
import subprocess
import sys
import tempfile

import numpy as np
from scipy.spatial import cKDTree
from skfem import CellBasis

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases  # noqa: E402
import run_skfem  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")
VARIANTS = [("linear", 1, False, "P1"), ("quadratic-straight", 2, False, "P2-straight"),
            ("quadratic-curved", 2, True, "P2-curved")]
# NDSolve runs: the first levels of cases.LEVELS (NIntegrate over large meshes is slow)
NLEVELS = {"P1solid": 3, "P1": 5, "P1c": 5, "P2": 5, "P3c": 4, "P3f": 4}
NLEVELS_QUAD_3D = 3


def wpath(p):
    return p.replace("\\", "/")


def job_entry(name, case, h, order, curved, workdir):
    c = cases.CASES[case]
    m = cases.make_mesh(case, h, order, curved)
    P, C = m["points"], m["cells"]
    stem = os.path.join(workdir, name)
    np.savetxt(stem + "_coords.csv", P, delimiter=",", fmt="%.17g")
    np.savetxt(stem + "_cells.csv", C + 1, delimiter=",", fmt="%d")
    pc = pcells = "None"
    if c["q2"]:
        cent = P[C[:, :3]].mean(axis=1).T
        sel = cases.in_p1_patch(cent) if c["q2"] == "p1" else cases.in_p2_patch(cent)
        sub = C[sel]
        used = np.unique(sub)
        remap = np.full(len(P), -1)
        remap[used] = np.arange(len(used))
        np.savetxt(stem + "_pcoords.csv", P[used], delimiter=",", fmt="%.17g")
        np.savetxt(stem + "_pcells.csv", remap[sub] + 1, delimiter=",", fmt="%d")
        pc = '"%s"' % wpath(stem + "_pcoords.csv")
        pcells = '"%s"' % wpath(stem + "_pcells.csv")
    bc = "{" + ", ".join("{%d, %r, %d, %r}" % (ax, float(val), comp, float(g)) for (ax, val), comp, g in c["bc"]) + "}"
    bc = bc.replace("e-", "*^-").replace("e+", "*^")
    out, nodal = wpath(stem + "_out.csv"), wpath(stem + "_nodal.csv")
    entry = ('<|"name" -> "%s", "dim" -> %d, "coords" -> "%s", "cells" -> "%s", "patchcoords" -> %s, '
             '"patchcells" -> %s, "bc" -> %s, "energyScale" -> %r, "out" -> "%s", "nodal" -> "%s"|>') % (
        name, c["dim"], wpath(stem + "_coords.csv"), wpath(stem + "_cells.csv"), pc, pcells, bc,
        float(c["energy_scale"]), out, nodal)
    return entry, out, nodal, P


def evaluate_nodal(case, h, skvariant, P, nodal_csv):
    """K and Q2 of NDSolve's nodal displacements, with the scikit-fem basis of the same mesh and element space."""
    c = cases.CASES[case]
    U = np.loadtxt(nodal_csv, delimiter=",", ndmin=2)
    mesh, e = run_skfem.build(case, h, skvariant)
    basis = CellBasis(mesh, e)
    tree = cKDTree(P)
    u = np.zeros(basis.N)
    for k, dofs in enumerate(basis.split_indices()):
        dist, idx = tree.query(basis.doflocs[:, dofs].T)
        if dist.max() > 1e-9:
            raise RuntimeError("degree of freedom without a mesh node (%g)" % dist.max())
        u[dofs] = U[idx, k]
    A = run_skfem.stiffness(basis, c["dim"])
    kc = c["energy_scale"] * float(u @ (A @ u))
    q2 = float("nan")
    if c["q2"]:
        q2, _ = run_skfem.patch_average(case, mesh, e, u)
    return kc, q2


def main(argv):
    sys.stdout.reconfigure(encoding="utf-8")
    todo = argv or list(cases.CASES)
    workdir = tempfile.mkdtemp(prefix="nd_bench_")
    jobs, meta = [], []
    for case in todo:
        c = cases.CASES[case]
        for variant, order, curved, skv in VARIANTS:
            if curved and c["mesh"] == "lbracket":
                continue
            n = NLEVELS[case] if not (c["dim"] == 3 and order == 2) else NLEVELS_QUAD_3D
            n = min(n, int(os.environ.get("ND_MAXLEVELS", "99")))
            for h in cases.LEVELS[case][:n]:
                name = "%s_%s_%g" % (case, variant, h)
                entry, out, nodal, P = job_entry(name, case, h, order, curved, workdir)
                jobs.append(entry)
                meta.append((case, variant, skv, h, out, nodal, P))
    jobfile = os.path.join(workdir, "jobs.wl")
    with open(jobfile, "w") as f:
        f.write("{\n" + ",\n".join(jobs) + "\n}\n")
    cmd = shlex.split(os.environ.get("WOLFRAMSCRIPT", "wolframscript"), posix=True)
    cmd += ["-file", wpath(os.path.join(HERE, "ndsolve_fem.wls")), wpath(jobfile)]
    for attempt in range(3):
        print("running (attempt %d):" % (attempt + 1), " ".join(cmd), flush=True)
        p = subprocess.run(cmd, capture_output=True, text=True)
        print(p.stdout[-8000:], p.stderr[-2000:], flush=True)
        if all(os.path.exists(m[4]) for m in meta):
            break
        # a job that stopped the kernel is skipped on the next attempt
        for m in meta:
            if not os.path.exists(m[4]):
                print("not done:", os.path.basename(m[4]), flush=True)
                if attempt == 1:
                    with open(m[4], "w") as f:
                        f.write("failed\nkernel stopped\n")
                break
    by_case = {}
    for case, variant, skv, h, out, nodal, P in meta:
        row = dict(code="NDSolve", case=case, element=variant, h=h)
        with open(out) as f:
            r = list(csv.reader(f))
        if r[0][0] == "failed":
            row["K"] = "failed (kernel stopped)"
        else:
            row.update(dict(zip(r[0], r[1])))
            row["K"] = row.pop("K_system")
            kc, q2 = evaluate_nodal(case, h, skv, P, nodal)
            row["K_check"] = "%.12g" % kc
            row["Q2"] = "%.12g" % q2
        by_case.setdefault(case, []).append(row)
        print(case, variant, h, row.get("K"), row.get("K_check"), row.get("K_nintegrate"), row.get("Q2"),
              row.get("Q2_nintegrate"), flush=True)
    fields = ["code", "case", "element", "h", "ndof", "nelem", "time_s", "ndsolvevalue_time_s", "integrate_time_s",
              "K", "K_check", "K_nintegrate", "Q2", "Q2_nintegrate", "patch_area_nintegrate", "mesh_measure"]
    for case, rows in by_case.items():
        with open(os.path.join(RESULTS, "ndsolve_%s.csv" % case), "w", newline="") as f:
            w = csv.DictWriter(f, fields, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
    print("work directory:", workdir)


if __name__ == "__main__":
    main(sys.argv[1:])
