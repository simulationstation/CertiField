"""Resource limits and resumable CSV logs shared by the harness runners.

Import this module first (before numpy / scipy / pypardiso): it caps the threads of the numerical libraries at
HARNESS_THREADS (default 2) unless they are set already.

CsvLog: every finished run is appended to its CSV file at once, and a rerun skips the runs whose row exists, so a stop
loses at most the run in progress.  A run that is not done because it would exceed the memory cap gets a row with
the reason in the column "note".
"""
import csv
import os

THREADS = os.environ.get("HARNESS_THREADS", "2")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, THREADS)

# memory cap per process (GB) and the measured cost of the direct solves (GB per million unknowns, with margin);
# runs above the cap are skipped and recorded as such
MAX_GB = float(os.environ.get("HARNESS_MAX_GB", "4"))


def peak_gb():
    try:
        import psutil
        p = psutil.Process()
        mi = p.memory_info()
        return getattr(mi, "peak_wset", mi.rss) / 2 ** 30
    except Exception:  # pragma: no cover
        return float("nan")


class CsvLog:
    def __init__(self, path, fields, key=("element", "h")):
        self.path, self.fields, self.key = path, list(fields), key
        if "note" not in self.fields:
            self.fields.append("note")
        self.rows = []
        if os.path.exists(path):
            with open(path, newline="") as f:
                self.rows = list(csv.DictReader(f))

    def _k(self, row):
        return tuple(str(row.get(k)) for k in self.key)

    def has(self, **kw):
        k = tuple(str(kw[k]) for k in self.key)
        return any(self._k(r) == k for r in self.rows)

    def append(self, row):
        row = {k: (("%.12g" % v) if isinstance(v, float) else v) for k, v in row.items()}
        new = not os.path.exists(self.path)
        with open(self.path, "a", newline="") as f:
            w = csv.DictWriter(f, self.fields, extrasaction="ignore")
            if new:
                w.writeheader()
            w.writerow(row)
        self.rows.append(row)
