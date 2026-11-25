"""Set your own values of Q1 and Q2 against the certified intervals of the benchmark pack.

Usage (from the repository root or anywhere; Python 3, no other packages):

    python benchpack/compare/compare.py P1 --Q1 0.59121 --Q2 1.03275
    python benchpack/compare/compare.py P3c --Q1 1.8329
    python benchpack/compare/compare.py --list

Problems: P1 (frictionless grips), P1c (clamped grips), P2, P3c / P3f (P3 with clamped / frictionless grips),
P1-nohole (the known answers K = 1, Q2 = 1/2).  Quantities and units are those of benchpack/problems/*.md
(E = 1, delta = 1; Q1 = K, Q2 = patch average of sigma_xx).  P3 has only Q1.

Exit status: 0 if every value given lies inside its interval, 1 if one lies outside, 2 on a usage error.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intervals  # noqa: E402

NOTE_Q1 = ("Note on Q1: the grips are displacement-loaded, so a conforming displacement finite-element model on a mesh\n"
           "that contains the body (e.g. straight edges with vertices on a hole) gives a value ABOVE the exact K.\n"
           "A value above the upper end should fall towards the interval as the mesh is refined (discretisation\n"
           "error, or a polygonal hole: see benchpack/results/VERDICTS.md). A value below the lower end points to an\n"
           "error: geometry (a mesh hole larger than the true hole), grip conditions, the quantity (e.g. the force of\n"
           "half the grip), non-conforming or under-integrated elements, inaccurate post-processing, or an\n"
           "unconverged solver.")


def fmt(x):
    return "%.10g" % float(x)


def report(key, values, data):
    ok = True
    lines = ["Problem %s: %s" % (key, intervals.DESCRIPTION[key])]
    for q, v in values.items():
        if v is None:
            continue
        if q not in data:
            lines.append("  %s = %s : no certified interval stored for this quantity" % (q, fmt(v)))
            continue
        lo, hi = data[q]
        mid = (float(lo) + float(hi)) / 2
        width = float(hi - lo)
        where, dist, rel = intervals.judge(v, data[q])
        lines.append("  %s = %s" % (q, fmt(v)))
        lines.append("     certified: %s <= %s <= %s   (width %.3g, relative %.3g)" % (
            fmt(lo), q, fmt(hi), width, width / abs(mid) if mid else float("nan")))
        if where == "inside":
            lines.append("     INSIDE the certified interval")
        else:
            ok = False
            lines.append("     OUTSIDE: %s the interval by %.3g (%.3g relative, %.3g interval widths)" % (
                where.upper(), dist, rel, dist / width if width else float("inf")))
    lines.append("  intervals from: " + "; ".join(data.get("sources", [])))
    if values.get("Q1") is not None:
        lines.append(NOTE_Q1)
    return ok, "\n".join(lines)


def main(argv=None):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except AttributeError:
        pass
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("problem", nargs="?", help="P1, P1c, P2, P3c, P3f or P1-nohole")
    ap.add_argument("--Q1", type=float, help="your value of Q1 = K")
    ap.add_argument("--Q2", type=float, help="your value of Q2 (patch average of sigma_xx)")
    ap.add_argument("--list", action="store_true", help="list the stored intervals")
    a = ap.parse_args(argv)
    data = intervals.load()
    if a.list or not a.problem:
        for key in ["P1", "P1c", "P2", "P3c", "P3f", "P1-nohole"]:
            d = data.get(key, {})
            qs = ", ".join("%s in [%s, %s]" % (q, fmt(d[q][0]), fmt(d[q][1])) for q in ("Q1", "Q2") if q in d)
            print("%-10s %s" % (key, qs if qs else "(no certificate stored yet)"))
        return 0
    try:
        key = intervals.resolve(a.problem)
    except KeyError as e:
        print(e)
        return 2
    if a.Q1 is None and a.Q2 is None:
        print("give --Q1 and/or --Q2")
        return 2
    if key not in data or not any(q in data[key] for q in ("Q1", "Q2")):
        print("No certified interval is stored yet for %s." % key)
        return 2
    ok, text = report(key, {"Q1": a.Q1, "Q2": a.Q2}, data[key])
    print(text)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
