"""Screening sweep: grid + random sampling + adversarial search, per target, in ONE process.

    bash run.sh sweep.py --targets scipy.special,torch.special [--n 2000]

A target is selected when its name starts with one of the given prefixes.  Records are
appended to audit/results/<target>.jsonl in batches as they are produced; the adversarial
search keeps a checkpoint (<target>.adv.json); a finished target gets <target>.done.  A
rerun skips what is done.  A point that runs longer than --timeout ends the process
(runner.py); run.sh restarts it and the point is recorded as a timeout.
"""
import argparse
import json
import math
import os
import sys
import time
import warnings
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import runner  # noqa: E402  (sets the thread limits before numpy is imported)
import numpy as np  # noqa: E402

RESULTS = os.environ.get('AUDIT_RESULTS') or os.path.join(os.path.dirname(HERE), 'results')
LOGS = os.path.join(os.path.dirname(RESULTS), 'logs')  # follows AUDIT_RESULTS (tests)


def fbase(name):
    return os.path.join(RESULTS, name.replace('[', '_').replace(']', ''))


def cast(t, pts):
    """Make the points exactly representable in the target dtype; drop invalid ones."""
    pts = np.asarray(pts, dtype=float).reshape(-1, t.nargs)
    if t.dtype == 'float32':
        with np.errstate(all='ignore'):
            p32 = pts.astype(np.float32)
        ok = np.all(np.isfinite(p32), axis=1) & np.all((p32 != 0) | (pts == 0), axis=1)
        pts = p32[ok].astype(np.float64)
    keep = np.array([bool(t.valid(*p)) for p in pts], dtype=bool) if len(pts) else np.zeros(0, bool)
    return pts[keep]


def call(t, pts):
    """Library values; returns (values, errors). errors[i] is None or a short reason when
    the library call itself failed (exception, or an external R/Julia call that did not
    return): the value is then NaN and the record says why."""
    out = np.empty(len(pts))
    err = [None] * len(pts)
    external = t.lib in ('R', 'julia')
    with warnings.catch_warnings(), np.errstate(all='ignore'):
        warnings.simplefilter('ignore')
        step = 500 if external else 4096
        for i in range(0, len(pts), step):
            _call_block(t, pts[i:i + step], out, err, i, external)
    return out, err


def _call_block(t, c, out, err, off, external):
    import targets
    if external:
        targets.EXT_TIMEOUT[0] = 30 + 0.05 * len(c)
    try:
        out[off:off + len(c)] = t.call(*[c[:, j] for j in range(t.nargs)])
        return
    except Exception as e:
        reason = 'timeout' if 'Timeout' in type(e).__name__ else type(e).__name__
    if len(c) == 1:
        out[off] = np.nan
        err[off] = 'library call failed: ' + reason
        return
    if not external:  # cheap in-process calls: point by point
        for k in range(len(c)):
            _call_block(t, c[k:k + 1], out, err, off + k, external)
        return
    h = len(c) // 2  # external: bisect, so one hanging point costs O(log n) calls
    _call_block(t, c[:h], out, err, off, external)
    _call_block(t, c[h:], out, err, off + h, external)


def perturb(rng, t, x, sigma):
    x = list(x)
    idx = [i for i, k in enumerate(t.kinds) if k != 'fixed']
    which = idx if rng.random() < 0.5 else [idx[rng.integers(len(idx))]]
    for i in which:
        k, v = t.kinds[i], x[i]
        g = rng.normal() * sigma
        if k == 'pos':
            v = v * math.exp(g) if v > 0 else 1e-300
        elif k == 'real':
            if v == 0:
                v = rng.normal() * 1e-300
            elif rng.random() < 0.05:
                v = v + g * (abs(v) + 1)
            else:
                v = v * math.exp(g)
        elif k == 'unit':
            v = v * math.exp(g) if v < 0.5 else 1 - (1 - v) * math.exp(g)
        elif k == 'p2':
            v = v * math.exp(g) if v < 1 else 2 - (2 - v) * math.exp(g)
        elif k == 'int':
            if rng.random() < 0.3:
                v = v + rng.choice([-1.0, 1.0])
            else:
                v = float(round(v * math.exp(g)))
        x[i] = v
    return tuple(x)


class Sink:
    """Appends records to the results file in batches."""

    def __init__(self, path):
        self.path, self.buf, self.t0 = path, [], time.time()

    def add(self, rec):
        self.buf.append(rec)
        if len(self.buf) >= 200 or time.time() - self.t0 > 60:
            self.flush()

    def flush(self):
        if self.buf:
            with open(self.path, 'a', encoding='utf-8') as f:
                for r in self.buf:
                    f.write(json.dumps(r, default=float) + '\n')
            self.buf = []
        self.t0 = time.time()


def load(path):
    if not os.path.exists(path):
        return []
    out = []
    with open(path, encoding='utf-8') as f:
        for line in f:
            try:
                out.append(json.loads(line))
            except ValueError:
                pass  # a line cut by a hard stop
    return out


def run_target(R, t, args, stuck):
    import tasks
    base = fbase(t.name)
    if os.path.exists(base + '.done'):
        return
    # two sweeps may run (one from each end of the list): a target belongs to whoever
    # started it first
    lock = base + '.lock'
    if os.path.exists(lock):
        with open(lock) as f:
            if f.read().strip() != args.owner:
                return
    else:
        with open(lock, 'w') as f:
            f.write(args.owner)
    t0 = time.time()
    rng = np.random.default_rng(args.seed * 1000003 + zlib.crc32(t.name.encode()))
    grid = cast(t, t.grid) if t.grid else np.zeros((0, t.nargs))
    with np.errstate(all='ignore'):
        rnd = cast(t, t.sample(rng, args.n))
    pts = np.concatenate([grid, rnd])
    phase = ['grid'] * len(grid) + ['random'] * len(rnd)
    path = base + '.jsonl'
    recs = load(path)
    done_i = {r['i'] for r in recs if 'i' in r}
    sink = Sink(path)
    if stuck and stuck.get('t') == t.name and stuck.get('ph') == 's' and stuck['i'] not in done_i:
        i = stuck['i']
        sink.add({'flag': 'timeout', 'x': [float(v) for v in pts[i]], 'phase': phase[i], 'i': i})
        sink.flush()  # must survive a further stop, or the point would hang again
        done_i.add(i)
    todo = [i for i in range(len(pts)) if i not in done_i]
    if todo:
        y, yerr = call(t, pts[todo])
        for i, yy, ye in zip(todo, y, yerr):
            p = tuple(float(v) for v in pts[i])
            rec = R.run({'t': t.name, 'ph': 's', 'i': i}, lambda: tasks.screen(t.name, p, float(yy)))
            rec = dict(rec or {})
            rec.update(x=list(p), y=float(yy), phase=phase[i], i=i)
            if ye:
                rec.update(lib_error=ye, cand=True, flag='lib_error')
            sink.add(rec)
        sink.flush()
    if args.adv_steps > 0:
        adversarial(R, t, args, stuck, path, base + '.adv.json')
    with open(base + '.done', 'w') as f:
        f.write('%s\n' % time.strftime('%Y-%m-%d %H:%M'))
    summarize(t, load(path), time.time() - t0)


def adversarial(R, t, args, stuck, path, ckpt):
    """(1+1) evolution strategy from the worst screening points (by score: the error in
    excess of what conditioning explains) and from random points."""
    import tasks
    if os.path.exists(ckpt):
        with open(ckpt, encoding='utf-8') as f:
            st = json.load(f)
    else:
        recs = [r for r in load(path) if r.get('phase') in ('grid', 'random')]
        if not recs:
            return
        arng = np.random.default_rng(args.seed * 7919 + zlib.crc32(t.name.encode()))
        order = sorted(range(len(recs)), key=lambda i: -(recs[i].get('score') or 0))
        starts, seen = [], set()
        for i in order:
            if recs[i].get('flag') in ('noref', 'timeout'):
                continue
            key = tuple(np.round(np.log10(np.abs(np.array(recs[i]['x'])) + 1e-300), 1))
            if key in seen:
                continue
            seen.add(key)
            starts.append(i)
            if len(starts) >= args.adv_starts * 2 // 3:
                break
        pool = [i for i in range(len(recs)) if recs[i].get('flag') not in ('noref', 'timeout')]
        extra = arng.choice(pool, min(len(pool), args.adv_starts - len(starts)), replace=False) if pool else []
        starts += [int(i) for i in extra]
        st = {'step': 0, 'cur': [[recs[i]['x'], recs[i].get('score') or 0.0] for i in starts],
              'sig': [0.3] * len(starts), 'rng': arng.bit_generator.state, 'skip': []}
    if stuck and stuck.get('t') == t.name and stuck.get('ph') == 'a':
        st['skip'].append([stuck['step'], stuck['j']])
        st['elapsed'] = st.get('elapsed', 0.0) + args.timeout  # the lost task's time
        with open(ckpt, 'w', encoding='utf-8') as f:
            json.dump(st, f)
    floor = 2.0 ** -22 if t.dtype == 'float32' else 1e-14
    skip = {tuple(s) for s in st['skip']}
    while st['step'] < args.adv_steps:
        if st.get('elapsed', 0.0) > args.adv_time:
            break  # wall-time cap for the adversarial search of one target
        t_step = time.time()
        s = st['step']
        rng = np.random.default_rng()
        rng.bit_generator.state = st['rng']
        props = []
        for (x, sc), sg in zip(st['cur'], st['sig']):
            for _ in range(30):
                q = perturb(rng, t, x, sg)
                qa = cast(t, [q])
                if len(qa):
                    props.append(tuple(float(v) for v in qa[0]))
                    break
            else:
                props.append(tuple(x))
        ya, yaerr = call(t, np.array(props, dtype=float))
        out = []
        for j, (p, v) in enumerate(zip(props, ya)):
            if (s, j) in skip:
                r = {'flag': 'timeout'}
            else:
                r = R.run({'t': t.name, 'ph': 'a', 'step': s, 'j': j},
                          lambda: tasks.screen(t.name, p, float(v)))
            r = dict(r or {})
            r.update(x=list(p), y=float(v), phase='adv', step=s)
            if yaerr[j]:
                r.update(lib_error=yaerr[j], cand=True, flag='lib_error', score=0.0)
            out.append(r)
            sc = r.get('score') or 0.0
            if sc > st['cur'][j][1]:
                st['cur'][j] = [list(p), sc]
                st['sig'][j] = min(st['sig'][j] * 1.6, 3.0)
            else:
                st['sig'][j] = max(st['sig'][j] * 0.8, floor)
                if st['sig'][j] <= floor * 1.0001:
                    st['sig'][j] = 0.3
        with open(path, 'a', encoding='utf-8') as f:
            for r in out:
                f.write(json.dumps(r, default=float) + '\n')
        st['step'] = s + 1
        st['rng'] = rng.bit_generator.state
        st['elapsed'] = st.get('elapsed', 0.0) + (time.time() - t_step)
        with open(ckpt, 'w', encoding='utf-8') as f:
            json.dump(st, f)


def summarize(t, recs, dt):
    n = len(recs)
    flags = {}
    for r in recs:
        if r.get('flag'):
            flags[r['flag']] = flags.get(r['flag'], 0) + 1
    good = [r for r in recs if r.get('flag') is None and r.get('relerr') is not None]
    mx = max(good, key=lambda r: r['relerr'], default=None)
    mxs = max(good, key=lambda r: r.get('score') or 0, default=None)
    ncand = sum(1 for r in recs if r.get('cand'))
    line = '%-40s n=%6d cand=%5d  max_rel=%.2e at %s  max_excess=%.2e at %s  flags=%s  %.0fs' % (
        t.name, n, ncand, mx['relerr'] if mx else float('nan'), mx['x'] if mx else None,
        (mxs.get('score') or 0) if mxs else float('nan'), mxs['x'] if mxs else None, flags, dt)
    if t.dtype == 'float32':
        mu = max(good, key=lambda r: r.get('ulps') or 0, default=None)
        line += '  max_ulps=%s at %s' % (mu.get('ulps') if mu else None, mu['x'] if mu else None)
    print(line, flush=True)
    os.makedirs(LOGS, exist_ok=True)
    with open(os.path.join(LOGS, 'sweep_summary.txt'), 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    ap = argparse.ArgumentParser()
    ap.add_argument('--targets', default='')
    ap.add_argument('--exclude', default='')
    ap.add_argument('--n', type=int, default=2000)
    ap.add_argument('--adv-starts', type=int, default=16)
    ap.add_argument('--adv-steps', type=int, default=40)
    ap.add_argument('--timeout', type=float, default=20)
    ap.add_argument('--adv-time', type=float, default=900, help='wall-time cap (s) of the adversarial search per target')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--owner', default='A', help='lock tag of this sweep (two sweeps may run)')
    ap.add_argument('--reverse', action='store_true', help='targets from the end of the list')
    args = ap.parse_args()
    import targets
    names = targets.names()
    pre = [p for p in args.targets.split(',') if p]
    exc = [p for p in args.exclude.split(',') if p]
    sel = [n for n in names if (not pre or any(n.startswith(p) for p in pre)) and not any(n.startswith(p) for p in exc)]
    if args.list:
        print('\n'.join(sel))
        return 0
    os.makedirs(RESULTS, exist_ok=True)
    if args.reverse:
        sel = sel[::-1]
    R = runner.Runner(os.environ.get('AUDIT_RUNNER_NAME', 'sweep'), args.timeout)
    stuck = R.take_timed_out()
    if stuck:
        print('resuming after a timeout at %s' % stuck, flush=True)
    for nm in sel:
        run_target(R, targets.get(nm), args, stuck)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except BaseException:
        import traceback
        traceback.print_exc()
        sys.exit(2)  # a real error: run.sh does not restart
