"""Plots of value against DOF with the certified interval shaded, and the verdict table.

    python benchpack/harness/summarize.py

Reads benchpack/results/*_<case>.csv and the stored intervals (benchpack/compare/intervals.py); writes
benchpack/results/<case>.png and benchpack/results/verdicts_table.md.

Verdict per code, case, quantity and element variant, from the sequence of values on the refined meshes:
  v = value on the finest mesh, d = change from the previous mesh, r = ratio of the last two changes,
  e = estimated remaining change = d r / (1 - r) if 0 < r < 0.8, else 10 d (slow or erratic convergence);
  e is at least 2e-7 |v| for CalculiX (7 significant digits in its output) and 1e-10 |v| otherwise.
  - "converges into the interval": v and v +- e lie inside the interval;
  - "lands outside by x": v lies outside and its distance x to the interval exceeds 2 e;
  - "not converged": otherwise.
These are FLOAT diagnostics of the codes, not certificates; checked by the agent only.
"""
import csv
import glob
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
RESULTS = os.path.join(PACK, "results")
sys.path.insert(0, os.path.join(PACK, "compare"))
import intervals as IV  # noqa: E402

CASES = ["P1solid", "P1", "P1c", "P2", "P3c", "P3f"]
INTERVAL_KEY = {"P1solid": "P1-nohole", "P1": "P1", "P1c": "P1c", "P2": "P2", "P3c": "P3c", "P3f": "P3f"}
TITLE = {"P1solid": "P1 without the hole (exact K = 1, Q2 = 1/2)", "P1": "P1, frictionless grips",
         "P1c": "P1c, clamped grips", "P2": "P2, L-bracket", "P3c": "P3, clamped grips",
         "P3f": "P3, frictionless grips"}
VARIANT = {  # element name in the CSV files -> (class, label)
    "P1": "lin", "CPS3": "lin", "C3D4": "lin", "linear": "lin",
    "P2-straight": "qs", "CPS6-straight": "qs", "C3D10-straight": "qs", "quadratic-straight": "qs",
    "P2-curved": "qc", "CPS6-curved": "qc", "C3D10-curved": "qc", "quadratic-curved": "qc",
}
VLABEL = {"lin": "linear", "qs": "quadratic, straight edges", "qc": "quadratic, curved edges",
          "sys": "quadratic, K = a(u_h, u_h)", "nint": "quadratic, K by NIntegrate (original script)"}
# categorical slots 1-3 of the reference palette (validated all-pairs), marker shape as second encoding
STYLE = {"lin": ("#257de0", "o"), "qs": ("#f47230", "s"), "qc": ("#00b894", "^"),
         "sys": ("#257de0", "D"), "nint": ("#f47230", "v")}
INK, INK2, GRID, AXIS = "#191919", "#666666", "#e7e7e7", "#c2c2c2"
CODES = ["scikit-fem", "CalculiX", "NDSolve", "NDSolve, own mesh"]


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def load(case):
    """{code: {variant: [(ndof, K, Q2), ...] sorted by ndof}}"""
    out = {}
    for path in sorted(glob.glob(os.path.join(RESULTS, "*_%s.csv" % case))):
        stem = os.path.basename(path)[:-len("_%s.csv" % case)]
        code = {"skfem": "scikit-fem", "calculix": "CalculiX", "ndsolve": "NDSolve"}.get(stem)
        if code is None:
            continue
        with open(path) as f:
            for r in csv.DictReader(f):
                v = VARIANT.get(r["element"])
                n = fnum(r.get("ndof"))
                if v is None or n is None:
                    continue
                out.setdefault(code, {}).setdefault(v, []).append((n, fnum(r.get("K")), fnum(r.get("Q2"))))
    if case in ("P3c", "P3f"):
        path = os.path.join(RESULTS, "ndsolve-toelementmesh_P3.csv")
        if os.path.exists(path):
            with open(path) as f:
                for r in csv.DictReader(f):
                    if r["case"] != case:
                        continue
                    n = fnum(r["ndof"])
                    d = out.setdefault("NDSolve, own mesh", {})
                    if fnum(r.get("K_system")) is not None:
                        d.setdefault("sys", []).append((n, fnum(r["K_system"]), None))
                    d.setdefault("nint", []).append((n, fnum(r.get("K_nintegrate", r.get("K"))), None))
    for code in out:
        for v in out[code]:
            out[code][v].sort()
    return out


def verdict(seq, interval, code):
    """seq: values on refined meshes (coarse to fine)."""
    seq = [s for s in seq if s is not None]
    if len(seq) < 3:
        return None
    v = seq[-1]
    d1, d0 = abs(seq[-1] - seq[-2]), abs(seq[-2] - seq[-3])
    r = d1 / d0 if d0 > 0 else 0.0
    e = d1 * r / (1 - r) if 0 < r < 0.8 else (10 * d1 if r >= 0.8 else d1)
    floor = (2e-7 if code == "CalculiX" else 1e-10) * abs(v)
    e = max(e, floor)
    res = dict(value=v, est=e)
    if interval is None:
        res["text"] = "no certified interval stored yet"
        return res
    lo, hi = float(interval[0]), float(interval[1])
    dist = lo - v if v < lo else (v - hi if v > hi else 0.0)
    side = "below" if v < lo else "above"
    if code == "CalculiX":
        resolution = 0.5 * 10 ** (math.floor(math.log10(abs(v))) - 6)   # ccx prints 7 significant digits
        if dist <= resolution and (hi - lo) < 2 * resolution:
            res["text"] = ("agrees to the 7 digits ccx prints (the interval is narrower than that); "
                           "remaining change about %.2g" % e)
            return res
    if lo <= v - e and v + e <= hi:
        res["text"] = "converges into the interval"
    elif v < lo or v > hi:
        if dist > 2 * e:
            res["text"] = "lands outside: %s by %.3g (%.2g relative)" % (side, dist, dist / abs(v))
        else:
            res["text"] = "not converged (%s by %.3g, remaining change about %.2g)" % (side, dist, e)
    else:
        res["text"] = "not converged (inside, remaining change about %.2g)" % e
    return res


def window(values, interval):
    vals = [v for v in values if v is not None]
    if interval is not None:
        lo, hi = float(interval[0]), float(interval[1])
        w = hi - lo
        if w == 0:
            w = max(1e-9, 1e-4 * abs(lo))
        return lo - 2.5 * w, hi + 2.5 * w
    if not vals:
        return None
    # no interval: centre on the finest values
    tail = sorted(vals)[len(vals) // 4: 3 * len(vals) // 4 + 1] or vals
    a, b = min(tail), max(tail)
    pad = max(b - a, 1e-6 * abs(b))
    return a - pad, b + pad


def plot_case(case, data, ivs, path):
    qs = ["Q1", "Q2"] if any(t[2] is not None for c in data.values() for s in c.values() for t in s) else ["Q1"]
    codes = [c for c in CODES if c in data]
    fig, axes = plt.subplots(len(qs), len(codes), figsize=(3.6 * len(codes), 3.0 * len(qs) + 0.6), squeeze=False)
    for qi, q in enumerate(qs):
        interval = ivs.get(q)
        allvals = [t[1 + qi] for c in codes for s in data[c].values() for t in s]
        win = window(allvals, interval)
        for ci, code in enumerate(codes):
            ax = axes[qi][ci]
            ax.set_facecolor("#ffffff")
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                ax.spines[sp].set_color(AXIS)
            ax.tick_params(colors=INK2, labelsize=8)
            ax.grid(True, color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)
            if interval is not None:
                lo, hi = float(interval[0]), float(interval[1])
                ax.axhspan(lo, hi, color="#d9d9d9", alpha=0.7, linewidth=0,
                           label="certified interval" if (qi == 0 and ci == 0) else None)
                ax.axhline(lo, color=INK2, linewidth=0.6)
                ax.axhline(hi, color=INK2, linewidth=0.6)
            for v, pts in data[code].items():
                pts = [(n, t[qi]) for n, *t in [(p[0], p[1], p[2]) for p in pts] if t[qi] is not None]
                if not pts:
                    continue
                col, mk = STYLE[v]
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                if win is not None:
                    yc = [min(max(y, win[0]), win[1]) for y in ys]
                else:
                    yc = ys
                ax.plot(xs, yc, color=col, linewidth=1.5, label=VLABEL[v], zorder=3)
                inside = [(x, y) for x, y, y0 in zip(xs, yc, ys) if y == y0]
                outside = [(x, y) for x, y, y0 in zip(xs, yc, ys) if y != y0]
                if inside:
                    ax.plot(*zip(*inside), linestyle="none", marker=mk, markersize=5.5, color=col,
                            markeredgecolor="#ffffff", markeredgewidth=0.8, zorder=4)
                if outside:
                    ax.plot(*zip(*outside), linestyle="none", marker=mk, markersize=5.5, markerfacecolor="none",
                            markeredgecolor=col, markeredgewidth=1.2, zorder=4)
            ax.set_xscale("log")
            if win is not None:
                ax.set_ylim(*win)
            ax.ticklabel_format(axis="y", useOffset=False)
            if qi == 0:
                ax.set_title(code, color=INK, fontsize=10)
            if ci == 0:
                ax.set_ylabel("Q1 = K" if q == "Q1" else "Q2 = mean of sigma_xx on the patch", color=INK2, fontsize=9)
            if qi == len(qs) - 1:
                ax.set_xlabel("degrees of freedom", color=INK2, fontsize=9)
            ax.legend(fontsize=7, frameon=False, labelcolor=INK2, loc="best")
    fig.suptitle("%s: finite element values (FLOAT) against the certified interval%s" % (
        TITLE[case], "" if ivs else " (no interval stored yet)"), color=INK, fontsize=11)
    fig.text(0.01, 0.005, "Open markers: value outside the plotted window, drawn at its edge.  "
             "Grey band: certified interval.  Data: benchpack/results/*_%s.csv" % case, color=INK2, fontsize=7)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(path, dpi=130, facecolor="#ffffff")
    plt.close(fig)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ivs_all = IV.load()
    lines = ["| case | quantity | code | elements | finest DOF | value (finest) | certified interval | verdict |",
             "|---|---|---|---|---|---|---|---|"]
    for case in CASES:
        data = load(case)
        if not data:
            continue
        ivs = ivs_all.get(INTERVAL_KEY[case], {})
        ivs = {q: ivs[q] for q in ("Q1", "Q2") if q in ivs}
        plot_case(case, data, ivs, os.path.join(RESULTS, "%s.png" % case))
        for qi, q in enumerate(["Q1", "Q2"]):
            for code in CODES:
                for v in ["lin", "qs", "qc", "sys", "nint"]:
                    pts = data.get(code, {}).get(v)
                    if not pts:
                        continue
                    seq = [p[1 + qi] for p in pts]
                    if all(s is None for s in seq):
                        continue
                    res = verdict(seq, ivs.get(q), code)
                    if res is None:
                        continue
                    iv = ivs.get(q)
                    ivs_txt = "[%.10g, %.10g]" % (float(iv[0]), float(iv[1])) if iv else "-"
                    lines.append("| %s | %s | %s | %s | %d | %.10g | %s | %s |" % (
                        case, q, code, VLABEL[v], int(pts[-1][0]), res["value"], ivs_txt, res["text"]))
    with open(os.path.join(RESULTS, "verdicts_table.md"), "w", encoding="utf-8") as f:
        f.write("# Verdict table (generated by benchpack/harness/summarize.py)\n\n")
        f.write("All values FLOAT; intervals CERTIFIED (P1-nohole: exact); checked by the agent only.\n\n")
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
