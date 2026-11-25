"""The stored certified intervals of the benchmark quantities.

Sources:
- P3 (block with a spherical cavity): the grid-24 certificates of Reference/Curved3D (benchpack/problems/P3.md).
- P1, P1c, P2: every certificate benchpack/certificates/*.json in the format of benchpack/certificates/FORMAT.md
  ("benchpack-certificate-1"); its "claims" are read.  When several certificates cover the same quantity, their
  intersection is used (each one is a guaranteed enclosure, so the intersection is one too).
- P1 without the hole (sanity case): K = 1 and Q2 = 1/2 exactly (P1.md, known answers).

Keys: P1 (frictionless grips), P1c (clamped grips), P2, P3c (clamped), P3f (frictionless), P1-nohole.
"""
import glob
import json
import os
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CERT_DIR = os.path.join(os.path.dirname(HERE), "certificates")

ALIASES = {"p1": "P1", "p1f": "P1", "p1-frictionless": "P1", "p1c": "P1c", "p1-clamped": "P1c", "p2": "P2",
           "p3c": "P3c", "p3-clamped": "P3c", "p3f": "P3f", "p3-frictionless": "P3f", "p3": "P3f",
           "p1-nohole": "P1-nohole", "p1solid": "P1-nohole"}

DESCRIPTION = {
    "P1": "gripped plate with a circular hole, frictionless grips (P1.md)",
    "P1c": "gripped plate with a circular hole, clamped grips (P1.md, variant P1c)",
    "P2": "L-bracket (P2.md)",
    "P3c": "block with a spherical cavity, clamped grips (P3.md)",
    "P3f": "block with a spherical cavity, frictionless grips (P3.md)",
    "P1-nohole": "P1 without the hole, frictionless grips: exact known answers (P1.md)",
}

BUILTIN = {
    "P3c": {"Q1": (Fraction("1.81732246"), Fraction("1.83348613"))},
    "P3f": {"Q1": (Fraction("1.74510549"), Fraction("1.75827651"))},
    "P1-nohole": {"Q1": (Fraction(1), Fraction(1)), "Q2": (Fraction(1, 2), Fraction(1, 2))},
}
BUILTIN_SOURCE = {
    "P3c": "Reference/Curved3D grid-24 certificate (P3.md); CERTIFIED, checked by the agent only",
    "P3f": "Reference/Curved3D grid-24 certificate (P3.md); CERTIFIED, checked by the agent only",
    "P1-nohole": "exact (uniform uniaxial stress 1/2); PROVED in P1.md",
}


def _key(problem, variant):
    if problem == "P1":
        return "P1c" if str(variant).lower().startswith("clamp") else "P1"
    if problem == "P2":
        return "P2"
    if problem == "P3":
        return "P3c" if str(variant).lower().startswith("clamp") else "P3f"
    return None


def load(cert_dir=CERT_DIR):
    """{key: {"Q1": (lo, hi), "Q2": (lo, hi), "sources": [...]}} with Fraction endpoints."""
    out = {k: dict(v, sources=[BUILTIN_SOURCE[k]]) for k, v in BUILTIN.items()}
    for path in sorted(glob.glob(os.path.join(cert_dir, "*.json"))):
        try:
            with open(path, encoding="utf-8") as f:
                d = json.load(f)
        except (OSError, ValueError):
            continue
        if not isinstance(d, dict) or d.get("format") != "benchpack-certificate-1" or "claims" not in d:
            continue
        key = _key(d.get("problem"), d.get("variant"))
        if key is None:
            continue
        entry = out.setdefault(key, {"sources": []})
        used = False
        for q in ("Q1", "Q2"):
            c = d["claims"].get(q)
            if not c:
                continue
            lo, hi = Fraction(c["lower"]), Fraction(c["upper"])
            if q in entry:
                lo, hi = max(lo, entry[q][0]), min(hi, entry[q][1])
                if lo > hi:
                    raise ValueError("certificates disagree for %s %s (empty intersection): %s" % (key, q, path))
            entry[q] = (lo, hi)
            used = True
        if used:
            entry["sources"].append("%s (%s)" % (os.path.relpath(path, os.path.dirname(os.path.dirname(HERE))),
                                                 d.get("label", "no label")))
    return out


def resolve(name):
    k = ALIASES.get(str(name).lower())
    if k is None:
        raise KeyError("unknown problem %r; use one of %s" % (name, ", ".join(sorted(DESCRIPTION))))
    return k


def judge(value, interval):
    """('inside' | 'below' | 'above', distance to the interval, distance relative to the interval's end)."""
    lo, hi = float(interval[0]), float(interval[1])
    if value < lo:
        return "below", lo - value, (lo - value) / abs(lo) if lo else float("inf")
    if value > hi:
        return "above", value - hi, (value - hi) / abs(hi) if hi else float("inf")
    return "inside", 0.0, 0.0
