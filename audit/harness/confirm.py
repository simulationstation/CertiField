"""Confirmation of screening candidates (spec 4.2), in one process, resumable.

Candidates of each target are grouped into clusters (kind of failure and decade of every
argument); the worst points of each cluster (by excess score and by relative error) are
recomputed with mpmath at 50 and 100 digits and cross-checked with Arb.  Output, appended
point by point: audit/results/confirm/<target>.jsonl (<target>.done when finished).

    bash run.sh confirm.py --targets scipy.special.stdtrit [--per-target 12]
"""
import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import runner  # noqa: E402
import numpy as np  # noqa: E402

RESULTS = os.environ.get('AUDIT_RESULTS') or os.path.join(os.path.dirname(HERE), 'results')
OUT = os.path.join(RESULTS, 'confirm')
DBL_MIN = sys.float_info.min
FLT_MIN = float(np.finfo(np.float32).tiny)


def fname(name):
    return name.replace('[', '_').replace(']', '')


def load(path):
    if not os.path.exists(path):
        return []
    out = []
    with open(path, encoding='utf-8') as f:
        for line in f:
            try:
                out.append(json.loads(line))
            except ValueError:
                pass
    return out


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

    def best(c):
        return max(((r.get('score') or 0),
                    (r['relerr'] if isinstance(r.get('relerr'), float) and math.isfinite(r['relerr']) else 0))
                   for r in c)
    ranked = sorted(clusters.items(), key=lambda kv: best(kv[1]), reverse=True)
    special = [kv for kv in ranked if kv[0][0] not in ('accuracy', 'subnormal_flush')]
    acc = [kv for kv in ranked if kv[0][0] == 'accuracy']
    sub = [kv for kv in ranked if kv[0][0] == 'subnormal_flush']
    picks = []
    for key, c in special[:per_target] + acc[:per_target] + sub[:2]:
        a = max(c, key=lambda r: (r.get('score') or 0))
        b = max(c, key=lambda r: (r.get('relerr') if isinstance(r.get('relerr'), float) and math.isfinite(r['relerr']) else 0))
        for r in ([a] if a is b else [a, b]):
            picks.append((key, len(c), r))
    return picks, len(clusters), len(cands)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    ap = argparse.ArgumentParser()
    ap.add_argument('--targets', default='')
    ap.add_argument('--per-target', type=int, default=12)
    ap.add_argument('--timeout', type=float, default=600)
    ap.add_argument('--redo-unverified', action='store_true',
                    help='drop unverified/timed-out confirmations of the selected targets and redo them')
    args = ap.parse_args()
    import targets
    import tasks
    pre = [p for p in args.targets.split(',') if p]
    sel = [n for n in targets.names() if not pre or any(n.startswith(p) for p in pre)]
    os.makedirs(OUT, exist_ok=True)
    R = runner.Runner('confirm', args.timeout)
    stuck = R.take_timed_out()
    for nm in sel:
        t = targets.get(nm)
        src = os.path.join(RESULTS, fname(nm))
        out = os.path.join(OUT, fname(nm))
        if args.redo_unverified and os.path.exists(out + '.done'):
            keep = [c for c in load(out + '.jsonl') if c.get('status') not in ('unverified', 'timeout')]
            with open(out + '.jsonl', 'w', encoding='utf-8') as f:
                for c in keep:
                    f.write(json.dumps(c, default=float) + '\n')
            os.remove(out + '.done')
        if os.path.exists(out + '.done') or not os.path.exists(src + '.done'):
            continue
        picks, ncl, nc = select(t, load(src + '.jsonl'), args.per_target)
        have = {r.get('k') for r in load(out + '.jsonl')}
        with open(out + '.jsonl', 'a', encoding='utf-8') as f:
            for k, (key, size, r) in enumerate(picks):
                if k in have:
                    continue
                if stuck and stuck.get('t') == nm and stuck.get('k') == k:
                    c = {'status': 'timeout'}
                else:
                    c = R.run({'t': nm, 'k': k}, lambda: tasks.confirm(nm, tuple(r['x']), r['y'], r.get('r')))
                c = dict(c or {})
                c.update(k=k, target=nm, cluster=repr(key), cluster_kind=key[0], cluster_size=size,
                         screen_relerr=r.get('relerr'), screen_flag=r.get('flag'), x=r['x'], y=r['y'])
                f.write(json.dumps(c, default=float) + '\n')
                f.flush()
        st = {}
        for c in load(out + '.jsonl'):
            st[c.get('status', '?')] = st.get(c.get('status', '?'), 0) + 1
        with open(out + '.done', 'w') as f:
            f.write('done\n')
        print('%-40s candidates=%d clusters=%d confirmed-runs=%d status=%s' % (nm, nc, ncl, len(picks), st), flush=True)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except BaseException:
        import traceback
        traceback.print_exc()
        sys.exit(2)
