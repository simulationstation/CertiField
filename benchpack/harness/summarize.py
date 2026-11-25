"""Plots of value against DOF with the certified interval shaded, and the verdict table.

    python benchpack/harness/summarize.py

Reads benchpack/results/*_<case>.csv and the stored intervals (benchpack/compare/intervals.py); writes
benchpack/results/<case>.png and benchpack/results/verdicts_table.md.

Verdict per code, case, quantity and element variant, from the values on the three finest meshes (N unknowns):
  v = value on the finest mesh; the change is modelled as C N^(-p), p fitted to the three values (monotone sequences
  only), and e = estimated remaining change |v - v_inf|; without a fit (non-monotone or p outside [0.05, 5])
  e = 10 x the last change.  e is at least half a unit in the 7th significant digit for CalculiX (its printed precision) and
  1e-10 |v| otherwise.
  - "converges into the interval": v and v +- e lie inside the interval;
  - "lands outside by x": v lies outside and its distance x to the interval exceeds 2 e;
  - "not converged": otherwise;
  - exact reference values (zero width): "reproduces the exact value" if |v - exact| <= 1e-11 |v| (rounding).
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
CODES = ["scikit-fem", "CalculiX", "NDSolve", "NDSolve, NIntegrate", "NDSolve, own mesh"]
# NDSolve:            a(u_h, u_h) of NDSolve's solution (assembled system) and Q2 of its nodal field
# NDSolve, NIntegrate: the same solution post-processed with NIntegrate of the InterpolatingFunction
# NDSolve, own mesh:   P3 meshed by ToElementMesh (Reference/Curved3D/ndsolve_compare.wls)


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
                if code == "NDSolve" and fnum(r.get("K_nintegrate")) is not None:
                    out.setdefault("NDSolve, NIntegrate", {}).setdefault(v, []).append(
                        (n, fnum(r.get("K_nintegrate")), fnum(r.get("Q2_nintegrate"))))
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


def remaining(ns, vs):
    """estimated |v_finest - v_limit| from the three finest (N, value) pairs, model v = v_inf + C N^(-p)"""
    (n1, v1), (n2, v2), (n3, v3) = zip(ns[-3:], vs[-3:])
    a, b = v2 - v1, v3 - v2
    if a == 0 or b == 0 or (a > 0) != (b > 0) or not (n1 < n2 < n3):
        return 10 * abs(b)

    def f(p):
        return (n1 ** -p - n2 ** -p) / (n2 ** -p - n3 ** -p) - a / b
    lo, hi = 0.05, 5.0
    if f(lo) * f(hi) > 0:
        return 10 * abs(b)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    p = 0.5 * (lo + hi)
    return abs(b) * n3 ** -p / (n2 ** -p - n3 ** -p)


def verdict(pts, interval, code):
    """pts: (N, value) on refined meshes (coarse to fine); values None are skipped."""
    pts = [(n, s) for n, s in pts if s is not None]
    if len(pts) < 3:
        return None
    ns, seq = [p[0] for p in pts], [p[1] for p in pts]
    v = seq[-1]
    e = remaining(ns, seq)
    if code == "CalculiX":   # half a unit in the 7th significant digit (ccx's printed precision)
        floor = 0.5 * 10 ** (math.floor(math.log10(abs(v))) - 6) if v != 0 else 0.0
    else:
        floor = 1e-10 * abs(v)
    e = max(e, floor)
    res = dict(value=v, est=e, ndof=ns[-1])
    if interval is None:
        res["text"] = "no certified interval stored yet"
        return res
    lo, hi = float(interval[0]), float(interval[1])
    dist = lo - v if v < lo else (v - hi if v > hi else 0.0)
    side = "below" if v < lo else "above"
    if hi - lo <= 1e-14 * max(1.0, abs(lo)) and code != "CalculiX":
        res["text"] = ("reproduces the exact value (difference %.1e)" % dist if dist <= 1e-11 * abs(v)
                       else "differs from the exact value by %.3g" % dist)
        return res
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
        # with an interval: relative deviation from its midpoint on a symmetric-log axis, linear inside the interval
        rel = None
        if interval is not None:
            lo, hi = float(interval[0]), float(interval[1])
            mid, hw = (lo + hi) / 2, (hi - lo) / 2
            hwr = hw / abs(mid)
            rel = (mid, max(hwr, 1e-14))
            win = None
        else:
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
            if rel is not None:
                ax.set_yscale("symlog", linthresh=rel[1], linscale=0.6)
                ax.axhspan(-hwr, hwr, color="#d9d9d9", alpha=0.8, linewidth=0,
                           label="certified interval" if (qi == 0 and ci == 0) else None)
                ax.axhline(-hwr, color=INK2, linewidth=0.6)
                ax.axhline(hwr, color=INK2, linewidth=0.6)
            for v, pts in data[code].items():
                pts = [(n, t[qi]) for n, *t in [(p[0], p[1], p[2]) for p in pts] if t[qi] is not None]
                if not pts:
                    continue
                col, mk = STYLE[v]
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                if rel is not None:
                    ys = [(y - rel[0]) / abs(rel[0]) for y in ys]
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
            if rel is not None:
                # inside the band only the midpoint gets a tick (the band is labelled by its edges' lines)
                lim = ax.get_ylim()
                ticks = [t for t in ax.get_yticks() if t == 0 or abs(t) >= 1.5 * rel[1]]
                ax.set_yticks(ticks)
                ax.set_ylim(*lim)
            if qi == 0:
                ax.set_title(code, color=INK, fontsize=10)
            if ci == 0:
                name = "Q1 = K" if q == "Q1" else "Q2"
                if rel is not None:
                    ax.set_ylabel("%s: (value - m) / m\nm = %.10g" % (name, rel[0]), color=INK2, fontsize=8)
                else:
                    ax.set_ylabel(name if q == "Q1" else "Q2 = mean of sigma_xx on the patch", color=INK2, fontsize=9)
            if qi == len(qs) - 1:
                ax.set_xlabel("degrees of freedom", color=INK2, fontsize=9)
            ax.legend(fontsize=7, frameon=False, labelcolor=INK2,
                      loc="center right" if code == "NDSolve, own mesh" else "upper right")
    fig.suptitle("%s: finite element values (FLOAT) against the certified interval%s" % (
        TITLE[case], "" if ivs else " (no interval stored yet)"), color=INK, fontsize=11)
    if ivs:
        note = ("Grey band: certified interval; m = its midpoint; y axis symmetric-log, linear inside the band.  "
                "Data: benchpack/results/*_%s.csv" % case)
    else:
        note = ("Open markers: value outside the plotted window, drawn at its edge.  "
                "Data: benchpack/results/*_%s.csv" % case)
    fig.text(0.01, 0.005, note, color=INK2, fontsize=7)
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
                    seq = [(p[0], p[1 + qi]) for p in pts]
                    if all(s[1] is None for s in seq):
                        continue
                    res = verdict(seq, ivs.get(q), code)
                    if res is None:
                        continue
                    iv = ivs.get(q)
                    ivs_txt = "[%.10g, %.10g]" % (float(iv[0]), float(iv[1])) if iv else "-"
                    lines.append("| %s | %s | %s | %s | %d | %.10g | %s | %s |" % (
                        case, q, code, VLABEL[v], int(res["ndof"]), res["value"], ivs_txt, res["text"]))
    with open(os.path.join(RESULTS, "verdicts_table.md"), "w", encoding="utf-8") as f:
        f.write("# Verdict table (generated by benchpack/harness/summarize.py)\n\n")
        f.write("All values FLOAT; intervals CERTIFIED (P1-nohole: exact); checked by the agent only.\n\n")
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
