"""P2: rigorous (exact rational) evaluation of stored fields.

Every cell of the L-bracket is a rectangle on which each field is one polynomial with rational coefficients
(python-flint fmpq_mpoly; stored floats are exact binary rationals); its energy density is integrated exactly.
So a(v,v), Q(tau), l(v), l(tau) and F(tau) are exact rationals.  The cells are split over worker processes.
"""
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as Fr

import numpy as np
from flint import fmpq

from exact import X, Y, ZERO, q, uni, local_to_x, rect_integral
import p2float as pf

C11, C12, C66 = [q(v) for v in (pf.D.C[0][0], pf.D.C[0][1], pf.D.C[2][2])]
NU = q(pf.D.NU)
CHI = uni(pf.CHI, Y)


class _Cache:
    def __init__(self, S, var):
        self.S, self.var, self.p = S, var, {}

    def __call__(self, e):
        if e not in self.p:
            S = self.S
            self.p[e] = [local_to_x(S.pieces[e][m], S.a[e], S.h[e], self.var) for m in range(S.p + 1)]
        return self.p[e]


def _spline_cell(cx, cy, Sx, Sy, coef, ex, ey):
    fx, fy = int(Sx.first[ex]), int(Sy.first[ey])
    px, py = cx(ex), cy(ey)
    tot = ZERO
    for m in range(len(px)):
        inner = ZERO
        for n in range(len(py)):
            c = float(coef[fx + m, fy + n])
            if c != 0.0:
                inner += q(c) * py[n]
        if not inner.is_zero():
            tot += px[m] * inner
    return tot


def _fr(v):
    return Fr(int(v.p), int(v.q))


def _work(args):
    S, fields, cells, patch = args
    csx, csy, cux, cuy = _Cache(S.sx, X), _Cache(S.sy, Y), _Cache(S.ux, X), _Cache(S.uy, Y)
    out = {}
    for name, f in fields.items():
        tot, ptot = fmpq(0), fmpq(0)
        if 'A' in f:
            for (ex, ey) in cells['u']:
                vx = _spline_cell(cux, cuy, S.ux, S.uy, f['A'], ex, ey)
                vy = _spline_cell(cux, cuy, S.ux, S.uy, f['B'], ex, ey) - (X + 1) * q(Fr(1, 2))
                exx, eyy = vx.derivative('x'), vy.derivative('y')
                gxy = vx.derivative('y') + vy.derivative('x')
                dens = C11 * (exx * exx + eyy * eyy) + 2 * C12 * exx * eyy + C66 * gxy * gxy
                box = (S.ux.a[ex], S.ux.b[ex], S.uy.a[ey], S.uy.b[ey])
                tot += rect_integral(dens, *box)
                if (ex, ey) in patch['u']:
                    ptot += rect_integral(C11 * exx + C12 * eyy, *box)
        else:
            al, be, ga = [q(float(v)) for v in f['par']]
            php = CHI * (al + be * X + ga * Y)
            for (ex, ey) in cells['s']:
                phi = _spline_cell(csx, csy, S.sx, S.sy, f['C'], ex, ey)
                if S.sy.a[ey] >= 0:
                    phi = phi + php
                txx = phi.derivative('y').derivative('y')
                tyy = phi.derivative('x').derivative('x')
                txy = -phi.derivative('x').derivative('y')
                dens = txx * txx + tyy * tyy - 2 * NU * txx * tyy + 2 * (1 + NU) * txy * txy
                box = (S.sx.a[ex], S.sx.b[ex], S.sy.a[ey], S.sy.b[ey])
                tot += rect_integral(dens, *box)
                if (ex, ey) in patch['s']:
                    ptot += rect_integral(txx, *box)
        out[name] = (str(tot), str(ptot))
    return out


def evaluate(S, fields, nproc=None, checkpoint=None, nchunk=16):
    """exact a (or Q) and patch integrals of every field: dict name -> dict(a|Q=Fraction, l=Fraction[, F=Fraction]).
    checkpoint: a directory; the partial sums of finished chunks are stored there and reused when rerun."""
    nproc = nproc or int(os.environ.get('BENCHPACK_NPROC', '2'))
    cu, cs = S.cells(S.ux, S.uy), S.cells(S.sx, S.sy)
    patch = dict(u=set(S.patch_cells(S.ux, S.uy)), s=set(S.patch_cells(S.sx, S.sy)))
    jobs = [(S, fields, dict(u=cu[k::nchunk], s=cs[k::nchunk]), patch) for k in range(nchunk)]
    results = {}
    if checkpoint:
        os.makedirs(checkpoint, exist_ok=True)
        for k in range(nchunk):
            fn = os.path.join(checkpoint, f'chunk_{k:03d}.json')
            if os.path.exists(fn):
                with open(fn) as fh:
                    results[k] = json.load(fh)
    todo = [k for k in range(nchunk) if k not in results]

    def done(k, res):
        results[k] = res
        if checkpoint:
            fn = os.path.join(checkpoint, f'chunk_{k:03d}.json')
            with open(fn + '.tmp', 'w') as fh:
                json.dump(res, fh)
            os.replace(fn + '.tmp', fn)

    if todo:
        try:
            with ProcessPoolExecutor(max_workers=nproc) as ex:
                futs = {ex.submit(_work, jobs[k]): k for k in todo}
                for fu in as_completed(futs):
                    done(futs[fu], fu.result())
        except Exception as err:                 # e.g. a worker could not start (memory): evaluate serially
            print(f"    parallel evaluation failed ({type(err).__name__}); evaluating serially", flush=True)
            for k in todo:
                if k not in results:
                    done(k, _work(jobs[k]))
    tot = {name: [Fr(0), Fr(0)] for name in fields}
    for k in range(nchunk):
        for name, (a, p) in results[k].items():
            tot[name][0] += Fr(a)
            tot[name][1] += Fr(p)
    out = {}
    inv = 1 / pf.PATCH_AREA
    for name, f in fields.items():
        if 'A' in f:
            out[name] = dict(a=tot[name][0], l=tot[name][1] * inv)
        else:
            out[name] = dict(Q=tot[name][0], l=tot[name][1] * inv, F=-Fr(float(f['par'][1])))
    return out
