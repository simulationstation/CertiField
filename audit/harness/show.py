"""Print the worst records of a results file.

    python show.py scipy.stats.t.isf [--by score|relerr|ulps] [--flag inf] [-n 20]
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), 'results')


def load(name):
    fn = os.path.join(RESULTS, name.replace('[', '_').replace(']', '') + '.jsonl')
    with open(fn, encoding='utf-8') as f:
        return [json.loads(l) for l in f]


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--by', default='score')
    ap.add_argument('--flag', default=None)
    ap.add_argument('--cand', action='store_true')
    ap.add_argument('-n', type=int, default=20)
    a = ap.parse_args()
    recs = load(a.name)
    if a.flag:
        recs = [r for r in recs if r.get('flag') == a.flag]
    if a.cand:
        recs = [r for r in recs if r.get('cand')]
    key = lambda r: (r.get(a.by) if isinstance(r.get(a.by), (int, float)) else -1)
    recs.sort(key=key, reverse=True)
    for r in recs[:a.n]:
        print('x=%s y=%r r=%s rel=%s ulps=%s kappa=%s score=%s flag=%s src=%s %s' % (
            r.get('x'), r.get('y'), r.get('r'), r.get('relerr'), r.get('ulps'), r.get('kappa'),
            r.get('score'), r.get('flag'), r.get('src'), r.get('err', '')))


if __name__ == '__main__':
    main()
