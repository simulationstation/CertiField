"""P1: rigorous evaluation of stored fields (exact rational polynomials per cell, Arb balls for the curved parts).

For a displacement field (A, B) and a stress field (gamma, C) as defined in certificates/FORMAT.md this module
returns ball enclosures of a_q(v,v), Q_q(tau), l(v), l(tau) and the exact F_q(tau).  Every polynomial integral
over a rectangle is an exact rational; cells cut by the circle use exact.disc_part_integral; the patch functional
uses the divergence theorem on the boundary of P_h (exact.arc_integral / exact.ray_integral).
"""
import math
from fractions import Fraction as Fr

import numpy as np
from flint import arb, fmpq

from exact import (CTX, X, Y, ZERO, q, uni, local_to_x, rect_integral, disc_part_integral, arc_integral,
                   ray_integral, APoint, eval_arb)
import p1float as pf

RHO = uni(pf.ETA, X) + uni(pf.ZETA, Y)
RHO2 = RHO * RHO
HALF = Fr(1, 2)
R1, R2 = pf.A_HOLE, pf.R2
TH0 = arb.pi() / 2 - arb(q(Fr(1, 10)))
C0, S0 = arb(q(Fr(1, 10))).sin(), arb(q(Fr(1, 10))).cos()          # cos(TH0), sin(TH0)
C11, C12, C66 = [q(v) for v in (pf.D.C[0][0], pf.D.C[0][1], pf.D.C[2][2])]
NU = q(pf.D.NU)


class CellPolys:
    """exact local basis polynomials of a Spline1D, in the physical variable"""

    def __init__(self, S, var):
        self.S = S
        self.p = [[local_to_x(S.pieces[e][m], S.a[e], S.h[e], var) for m in range(S.p + 1)] for e in range(S.ne)]


def spline_cell(cpx, cpy, coef, ex, ey):
    """exact polynomial of sum coef[i,j] B_i(x) B'_j(y) on cell (ex, ey)"""
    fx, fy = int(cpx.S.first[ex]), int(cpy.S.first[ey])
    px, py = cpx.p[ex], cpy.p[ey]
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


class Rig:
    """rigorous evaluator for one Setup"""

    def __init__(self, S):
        self.S = S
        self.csx, self.csy = CellPolys(S.sx, X), CellPolys(S.sy, Y)
        self.cux, self.cuy = CellPolys(S.ux, X), CellPolys(S.uy, Y)
        self.scells = S.cells(S.sx, S.sy)
        self.ucells = S.cells(S.ux, S.uy)

    # ------------------------------------------------------------------ stress
    def g_poly(self, gamma):
        return sum((q(float(g)) * (Y - q(self.S.gbase)) ** e for g, e in zip(gamma, self.S.gexp)), ZERO)

    def g_y_poly(self, gamma):
        """d Phi_g / dy = sum gamma_k (y - b)^(e_k + 1)/(e_k + 1) for y >= b = g_base"""
        return sum((q(float(g)) * (Y - q(self.S.gbase)) ** (e + 1) * fmpq(1, e + 1) for g, e in zip(gamma, self.S.gexp)), ZERO)

    def airy_phi(self, C, ex, ey):
        psi = spline_cell(self.csx, self.csy, C, ex, ey)
        return RHO2 * psi if self.S.hole else psi

    def stress_cell(self, C, gamma, ex, ey, gpoly):
        phi = self.airy_phi(C, ex, ey)
        txx = phi.derivative('y').derivative('y')
        tyy = phi.derivative('x').derivative('x')
        txy = -phi.derivative('x').derivative('y')
        if self.S.sy.a[ey] >= self.S.gbase:
            txx = txx + gpoly
        return txx, tyy, txy

    def Q(self, C, gamma):
        """ball enclosure of Q_q(tau) = int_{Omega_q} tau.S.tau"""
        gpoly = self.g_poly(gamma)
        exact_sum = fmpq(0)
        ball = arb(0)
        S = self.S
        for ex, ey, kind in self.scells:
            txx, tyy, txy = self.stress_cell(C, gamma, ex, ey, gpoly)
            dens = txx * txx + tyy * tyy - 2 * NU * txx * tyy + 2 * (1 + NU) * txy * txy
            x0, x1, y0, y1 = S.sx.a[ex], S.sx.b[ex], S.sy.a[ey], S.sy.b[ey]
            exact_sum += rect_integral(dens, x0, x1, y0, y1)
            if kind == 1:
                ball -= disc_part_integral(dens, x0, x1, y0, y1)
        return arb(exact_sum) + ball

    def F(self, gamma):
        """exact F_q = int_0^1 tau_xx(1, y) dy (Fraction)"""
        return sum((Fr(float(g)) * (1 - self.S.gbase) ** (e + 1) / (e + 1) for g, e in zip(gamma, self.S.gexp)), Fr(0))

    def l_tau(self, C, gamma):
        """ball enclosure of l(tau) = (1/|P_h|) int_{P_h} tau_xx = (1/|P_h|) oint Phi_y n_y ds"""
        gy = self.g_y_poly(gamma)
        tot = arb(0)
        pieces = ['outer', 'ray'] if self.S.hole else ['outer', 'ray', 'inner']
        for name in pieces:
            for (ex, ey, lo, hi) in self.boundary_pieces(self.S.sx, self.S.sy, name):
                phi = self.airy_phi(C, ex, ey)
                Fy = phi.derivative('y')
                if self.S.sy.a[ey] >= self.S.gbase:
                    Fy = Fy + gy
                tot += self._piece_integral(name, Fy, 0, 1, lo, hi)
        return tot / arb(q(pf.PH_AREA))

    # ------------------------------------------------------------------ displacement
    def disp_cell(self, Am, Bm, ex, ey):
        vx = spline_cell(self.cux, self.cuy, Am, ex, ey) + X * q(HALF)
        vy = spline_cell(self.cux, self.cuy, Bm, ex, ey)
        return vx, vy

    def a(self, Am, Bm):
        """ball enclosure of a_q(v,v) = int_{Omega_q} eps.C.eps"""
        S = self.S
        exact_sum = fmpq(0)
        ball = arb(0)
        for ex, ey, kind in self.ucells:
            vx, vy = self.disp_cell(Am, Bm, ex, ey)
            exx, eyy = vx.derivative('x'), vy.derivative('y')
            gxy = vx.derivative('y') + vy.derivative('x')
            dens = C11 * (exx * exx + eyy * eyy) + 2 * C12 * exx * eyy + C66 * gxy * gxy
            x0, x1, y0, y1 = S.ux.a[ex], S.ux.b[ex], S.uy.a[ey], S.uy.b[ey]
            exact_sum += rect_integral(dens, x0, x1, y0, y1)
            if kind == 1:
                ball -= disc_part_integral(dens, x0, x1, y0, y1)
        return arb(exact_sum) + ball

    def l_v(self, Am, Bm):
        """ball enclosure of l(v) = (1/|P_h|) oint (C11 v_x n_x + C12 v_y n_y) ds  (v_x = 0 on x = 0)"""
        tot = arb(0)
        for name in ('outer', 'inner', 'ray'):
            for (ex, ey, lo, hi) in self.boundary_pieces(self.S.ux, self.S.uy, name):
                vx, vy = self.disp_cell(Am, Bm, ex, ey)
                tot += self._piece_integral(name, vx, C11, 0, lo, hi)
                tot += self._piece_integral(name, vy, 0, C12, lo, hi)
        return tot / arb(q(pf.PH_AREA))

    # ------------------------------------------------------------------ boundary of P_h
    def _piece_integral(self, name, F, wx, wy, lo, hi):
        """int over a boundary piece of F (wx n_x + wy n_y) ds"""
        if name == 'outer':        # n = (c, s), ds = R2 dt
            return arc_integral(F, R2, q(wx) * q(R2), q(wy) * q(R2), lo, hi)
        if name == 'inner':        # n = -(c, s), ds = R1 dt
            return arc_integral(F, R1, -q(wx) * q(R1), -q(wy) * q(R1), lo, hi)
        # ray at angle TH0: n = (sin TH0, -cos TH0) = (S0, -C0), ds = dr
        return ray_integral(F, C0, S0, arb(q(wx)) * S0 - arb(q(wy)) * C0, lo, hi)

    def boundary_pieces(self, SX, SY, name):
        """split a boundary piece of P_h at the cell lines: list of (ex, ey, lo, hi) with lo/hi APoints (arcs) or
        balls r (ray); the order of the break points is decided with balls (an overlap raises)"""
        xs = sorted(set(SX.t))
        ys = sorted(set(SY.t))
        if name in ('outer', 'inner'):
            R = R2 if name == 'outer' else R1
            pts = [APoint(TH0, C0, S0), APoint(arb.pi() / 2, arb(0), arb(1))]
            for xv in xs:
                if 0 < xv < R:
                    p = APoint.from_x(R, xv)
                    if p.t > TH0:
                        pts.append(p)
                    elif not p.t < TH0:
                        raise ArithmeticError("break point not separated from the end of the arc")
            for yv in ys:
                if 0 < yv < R:
                    p = APoint.from_y(R, yv)
                    if p.t > TH0:
                        pts.append(p)
                    elif not p.t < TH0:
                        raise ArithmeticError("break point not separated from the end of the arc")
            key = [float(p.t.mid()) for p in pts]
            order = np.argsort(key)
            pts = [pts[i] for i in order]
            out = []
            for lo, hi in zip(pts[:-1], pts[1:]):
                if not lo.t < hi.t:
                    raise ArithmeticError("arc break points not separated")
                tm = (lo.t + hi.t) / 2
                ex, ey = self._locate(SX, SY, arb(q(R)) * tm.cos(), arb(q(R)) * tm.sin())
                out.append((ex, ey, lo, hi))
            return out
        r0, r1 = arb(q(R1)), arb(q(R2))
        rs = [r0, r1]
        for xv in xs:
            if xv > 0:
                r = arb(q(xv)) / C0
                if r > r0 and r < r1:
                    rs.append(r)
                elif r.overlaps(r0) or r.overlaps(r1):
                    raise ArithmeticError("ray break point not separated")
        for yv in ys:
            if yv > 0:
                r = arb(q(yv)) / S0
                if r > r0 and r < r1:
                    rs.append(r)
                elif r.overlaps(r0) or r.overlaps(r1):
                    raise ArithmeticError("ray break point not separated")
        rs = [rs[i] for i in np.argsort([float(r.mid()) for r in rs])]
        out = []
        for lo, hi in zip(rs[:-1], rs[1:]):
            if not lo < hi:
                raise ArithmeticError("ray break points not separated")
            rm = (lo + hi) / 2
            ex, ey = self._locate(SX, SY, rm * C0, rm * S0)
            out.append((ex, ey, lo, hi))
        return out

    @staticmethod
    def _locate(SX, SY, xm, ym):
        """the cell that contains the ball point (xm, ym) in its interior (decided with balls).  Used for the interior
        point of a boundary piece between two consecutive break points: such a piece crosses no cell line, so it
        lies in the closed cell that contains one of its interior points."""
        ex, ey = int(SX.cell_of(float(xm.mid()))), int(SY.cell_of(float(ym.mid())))
        if not (xm > arb(q(SX.a[ex])) and xm < arb(q(SX.b[ex])) and ym > arb(q(SY.a[ey])) and ym < arb(q(SY.b[ey]))):
            raise ArithmeticError("cannot place a boundary piece in a cell")
        return ex, ey
