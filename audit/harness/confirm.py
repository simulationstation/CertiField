"""Confirmation of screening candidates (spec 4.2).

Candidates of each target are grouped into clusters (flag type and decade of every
argument); the worst points of each cluster (by excess score and by relative error) are
recomputed with mpmath at 50 and 100 digits and cross-checked with Arb.  Output:
audit/results/confirm/<target>.jsonl.

    python confirm.py --targets scipy.special.stdtrit [--per-target 30] [--workers 12]
"""
import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), 'results')
OUT = os.path.join(RESULTS, 'confirm')
DBL_MIN = sys.float_info.min
FLT_MIN = float(np.finfo(np.float32).tiny)


def fname(name):
    return name.replace('[', '_').replace(']', '')


def load(name):
    fn = os.path.join(RESULTS, fname(name) + '.jsonl')
    if not os.path.exists(fn):
        return []
    with open(fn, encoding='utf-8') as f:
        return [json.loads(l) for l in f]


def bucket(t, x):
    key = []
    for k, v in zip(t.kinds, x):
        if v == 0 or not math.isfinite(v):
            key.append(('z', v))
        elif k == 'unit' and v > 0.5:
            key.append(('1-', int(math.floor(math.log10(max(1 - v, 1e-300))))))
        elif k == 'fixed':
            key.append(('f', v))
        else:
            key.append((v > 0, int(math.floor(math.log10(abs(v))))))
    return tuple(key)


def kind_of(t, r):
    f = r.get('flag')
    if f == 'zero':
        try:
            rv = abs(float(r['r']))
        except (TypeError, ValueError):
            rv = 0.0
        tiny = DBL_MIN if t.dtype == 'float64' else FLT_MIN
        return 'subnormal_flush' if rv < tiny else 'zero'
    return f or 'accuracy'


def select(t, recs, per_target):
    cands = [r for r in recs if r.get('cand')]
    clusters = {}
    for r in cands:
        key = (kind_of(t, r), bucket(t, r['x']))
        clusters.setdefault(key, []).append(r)
    # rank clusters by their best excess score, then relative error
    def best(c):
        return max(((r.get('score') or 0), (r.get('relerr') or 0)) for r in c)
    ranked = sorted(clusters.items(), key=lambda kv: best(kv[1]), reverse=True)
    picks = []
    for key, c in ranked[:per_target]:
        a = max(c, key=lambda r: (r.get('score') or 0))
        b = max(c, key=lambda r: (r.get('relerr') if isinstance(r.get('relerr'), float) and math.isfinite(r['relerr']) else 0))
        for r in {id(a): a, id(b): b}.values():
            picks.append((key, len(c), r))
    return picks, len(clusters), len(cands)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    ap = argparse.ArgumentParser()
    ap.add_argument('--targets', default='')
    ap.add_argument('--per-target', type=int, default=30)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--timeout', type=float, default=600)
    args = ap.parse_args()
    import targets
    from pool import Pool
    pre = [p for p in args.targets.split(',') if p]
    sel = [n for n in targets.names() if not pre or any(n.startswith(p) for p in pre)]
    os.makedirs(OUT, exist_ok=True)
    pool = Pool(args.workers, args.timeout)
    try:
        for nm in sel:
            t = targets.get(nm)
            recs = load(nm)
            if not recs:
                continue
            picks, ncl, nc = select(t, recs, args.per_target)
            if not picks:
                print('%-40s no candidates (%d points)' % (nm, len(recs)), flush=True)
                with open(os.path.join(OUT, fname(nm) + '.jsonl'), 'w', encoding='utf-8') as f:
                    pass
                continue
            res = pool.map([('confirm', nm, tuple(r['x']), r['y']) for _, _, r in picks])
            with open(os.path.join(OUT, fname(nm) + '.jsonl'), 'w', encoding='utf-8') as f:
                for (key, size, r), c in zip(picks, res):
                    c = dict(c or {})
                    c.update(target=nm, cluster=repr(key), cluster_kind=key[0], cluster_size=size,
                             screen_relerr=r.get('relerr'), screen_flag=r.get('flag'), x=r['x'])
                    f.write(json.dumps(c, default=float) + '\n')
            st = {}
            for c in res:
                s = (c or {}).get('status', (c or {}).get('flag', '?'))
                st[s] = st.get(s, 0) + 1
            print('%-40s points=%d candidates=%d clusters=%d confirmed-runs=%d status=%s' % (
                nm, len(recs), nc, ncl, len(res), st), flush=True)
    finally:
        pool.close()


if __name__ == '__main__':
    main()
