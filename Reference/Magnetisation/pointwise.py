#!/usr/bin/env python3
"""Independent exact test of the pointwise inequality behind Theorem 4 of docs/MAGNETISATION.md.

At random rational points y of the unit sphere and random rational tangent vectors p, q, the quantity
    f(y, p, q) - N(y, p, q) + Gamma
is evaluated in exact rational arithmetic from the definitions in three-dimensional coordinates (the fields of a
certificate are differentiated by hand).  It shares no code with the checker of magnet.py except the evaluation of a
polynomial at a point.  For a certificate that the checker accepts it must be >= 0 everywhere."""
import random
from fractions import Fraction as Fr

import magnet as M

rnd = random.Random(20261123)


def sphere_point():
    """a random rational point of the unit sphere and two rational tangent vectors"""
    while True:
        a, b = Fr(rnd.randint(-40, 40), rnd.randint(1, 12)), Fr(rnd.randint(-40, 40), rnd.randint(1, 12))
        d = 1 + a * a + b * b
        y = (2 * a / d, 2 * b / d, (1 - a * a - b * b) / d)
        if abs(y[2]) != 1:
            break

    def tangent():
        v = [Fr(rnd.randint(-30, 30), rnd.randint(1, 9)) for _ in range(3)]
        dot = sum(v[i] * y[i] for i in range(3))
        return tuple(v[i] - dot * y[i] for i in range(3))

    return y, tangent(), tangent()


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return sum(x * z for x, z in zip(a, b))


def pointwise_defect(cert, y, p, q):
    """f(y, p, q) - N(y, p, q) + Gamma for the fields of the certificate, in exact arithmetic, from the definitions:
    f = A(|p|^2 + |q|^2) + kappa (y x e1).p + kappa (y x e2).q + V(y3),
    N = (u(y3) - 2A) y.(p x q) + grad S1 . p + grad S2 . q,   S = phi(y3) (y2, -y1)."""
    m = cert.model
    t = y[2]
    i = next(k for k in range(len(cert.u)) if cert.mesh[k] <= t <= cert.mesh[k + 1])
    phi, dphi, u = M.pval(cert.phi[i], t), M.pval(M.pder(cert.phi[i]), t), M.pval(cert.u[i], t)
    e1, e2 = (Fr(1), Fr(0), Fr(0)), (Fr(0), Fr(1), Fr(0))
    f = m.A * (dot(p, p) + dot(q, q)) + m.kappa * (dot(cross(y, e1), p) + dot(cross(y, e2), q)) + m.K * (1 - t * t) + m.h * (1 - t)
    # gradients in R^3 of S1 = phi(y3) y2 and S2 = -phi(y3) y1 (p, q are tangent, so the ambient gradient may be used)
    g1 = (Fr(0), phi, dphi * y[1])
    g2 = (-phi, Fr(0), -dphi * y[0])
    N = (u - 2 * m.A) * dot(y, cross(p, q)) + dot(g1, p) + dot(g2, q)
    return f - N + cert.Gamma


def pointwise_ok(cert, n=400):
    worst = None
    for _ in range(n):
        y, p, q = sphere_point()
        d = pointwise_defect(cert, y, p, q)
        worst = d if worst is None else min(worst, d)
    return worst >= 0, worst
