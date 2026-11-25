"""CalculiX input writer, runner and .dat parser used by run_calculix.py."""
import os
import re
import shutil
import subprocess

import numpy as np

CCX_ENV = os.environ.get("CCX_PREFIX", os.path.expanduser(r"~\miniforge3\envs\ccx"))


def ccx_command():
    exe = shutil.which("ccx") or shutil.which("ccx.exe")
    if exe:
        return exe, None
    lib = os.path.join(CCX_ENV, "Library")
    exe = os.path.join(lib, "bin", "ccx.exe")
    if not os.path.exists(exe):
        raise FileNotFoundError("ccx not found; install with: conda create -n ccx -c conda-forge calculix")
    path = os.pathsep.join([os.path.join(lib, "bin"), os.path.join(lib, "mingw-w64", "bin"),
                            os.path.join(lib, "usr", "bin"), os.environ.get("PATH", "")])
    return exe, path


def _lines(ids, per=8):
    ids = list(ids)
    return "\n".join(", ".join(str(i) for i in ids[k:k + per]) for k in range(0, len(ids), per))


def write_inp(path, points, cells, eltype, bcs, grip_nodes, patch_elems, thickness=1.0, E=1.0, nu=0.3,
              el_vars=("S", "COORD"), energy=True):
    """points (n, d), cells (m, k) in Abaqus order (0-based); bcs = list of (node index array, dof 1..3, value).
    Numbers are written with 13 significant digits: ccx reads at most 20 characters per field."""
    n, d = points.shape
    pts = np.where(np.abs(points) < 1e-14, 0.0, points)
    with open(path, "w") as f:
        f.write("*HEADING\nbenchpack harness\n*NODE, NSET=NALL\n")
        for i, p in enumerate(pts):
            if d == 2:
                f.write("%d, %.13g, %.13g\n" % (i + 1, p[0], p[1]))
            else:
                f.write("%d, %.13g, %.13g, %.13g\n" % (i + 1, p[0], p[1], p[2]))
        f.write("*ELEMENT, TYPE=%s, ELSET=EALL\n" % eltype)
        for j, c in enumerate(cells):
            f.write("%d, %s\n" % (j + 1, ", ".join(str(int(k) + 1) for k in c)))
        for k, (nodes, dof, val) in enumerate(bcs):
            f.write("*NSET, NSET=BC%d\n%s\n" % (k, _lines(np.asarray(nodes) + 1)))
        f.write("*NSET, NSET=GRIP\n%s\n" % _lines(np.asarray(grip_nodes) + 1))
        if patch_elems is not None and len(patch_elems):
            f.write("*ELSET, ELSET=PATCH\n%s\n" % _lines(np.asarray(patch_elems) + 1))
        f.write("*MATERIAL, NAME=M\n*ELASTIC\n%.13g, %.13g\n" % (E, nu))
        if d == 2:
            f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=M\n%.13g\n" % thickness)
        else:
            f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=M\n")
        f.write("*STEP\n*STATIC, SOLVER=PARDISO\n" if os.environ.get("CCX_PARDISO") else "*STEP\n*STATIC\n")
        f.write("*BOUNDARY\n")
        for k, (nodes, dof, val) in enumerate(bcs):
            f.write("BC%d, %d, %d, %.13g\n" % (k, dof, dof, val))
        f.write("*NODE PRINT, NSET=GRIP, TOTALS=ONLY\nRF\n")
        if patch_elems is not None and len(patch_elems):
            f.write("*EL PRINT, ELSET=PATCH\n%s\n" % ", ".join(el_vars))
        if energy:
            f.write("*EL PRINT, ELSET=EALL, TOTALS=ONLY\nELSE\n")
        f.write("*END STEP\n")


def run(jobdir, job, nthreads=1, timeout=7200):
    exe, path = ccx_command()
    env = dict(os.environ)
    if path:
        env["PATH"] = path
    env["OMP_NUM_THREADS"] = str(nthreads)
    env["CCX_NPROC_EQUATION_SOLVER"] = str(nthreads)
    p = subprocess.run([exe, "-i", job], cwd=jobdir, env=env, capture_output=True, text=True, timeout=timeout)
    log = p.stdout + p.stderr
    with open(os.path.join(jobdir, job + ".log"), "w") as f:
        f.write(log)
    if p.returncode != 0 or "*ERROR" in log:
        raise RuntimeError("ccx failed:\n" + log[-3000:])
    return log


_num = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def parse_dat(path):
    """Return {'RF': (fx, fy, fz) total over GRIP, 'S': {elem: [[ip, s11..s23], ...]}, 'COORD': {...}}."""
    out = {"S": {}, "COORD": {}, "EVOL": {}}
    with open(path) as f:
        lines = f.read().splitlines()
    mode = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if not s:
            continue
        low = s.lower()
        if low.startswith("total force"):
            mode = "RF"
            continue
        if low.startswith("total internal energy"):
            mode = "ELSE"
            continue
        if low.startswith("stresses (elem"):
            mode = "S"
            continue
        if low.startswith("global coordinates"):
            mode = "COORD"
            continue
        if low.startswith("volume (element"):
            mode = "EVOL"
            continue
        if low.startswith("displacements") or low.startswith("forces") or low.startswith("step") \
                or low.startswith("increment") or low.startswith("s t e p") or low.startswith("i n c"):
            mode = None
            continue
        nums = re.findall(_num, s)
        if mode == "RF" and len(nums) == 3 and "RF" not in out:
            out["RF"] = tuple(float(x) for x in nums)
        elif mode == "S" and len(nums) == 8:
            e = int(nums[0])
            out["S"].setdefault(e, []).append([float(x) for x in nums[1:]])
        elif mode == "COORD" and len(nums) == 5:
            e = int(nums[0])
            out["COORD"].setdefault(e, []).append([float(x) for x in nums[1:]])
        elif mode == "ELSE" and len(nums) == 1 and "ELSE" not in out:
            out["ELSE"] = float(nums[0])
        elif mode == "EVOL" and len(nums) == 2:
            out["EVOL"][int(nums[0])] = float(nums[1])
    return out
