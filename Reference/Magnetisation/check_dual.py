#!/usr/bin/env python3
"""Known answers and falsification for the exact lower bound (Proposition 5 and Proposition 7 of docs/MAGNETISATION.md).

A. Certificates with constant phi whose bounds are known exactly:
     exchange only (Belavin-Polyakov), the critical family (exact infimum 4 pi (1 - 2 r^2) in degree -1),
     the simple bound for the Zeeman model, the wall of the infinite line.
B. The checker must REJECT every certificate that would claim more than a known infimum: the same fields with u
   lowered, with phi perturbed, and random ones.
C. Independent exact test of the pointwise inequality behind Theorem 4: at random rational points of the sphere and
   random rational tangent vectors, with the fields differentiated by hand in three-dimensional coordinates (no code
   shared with the checker),  f - N + Gamma >= 0  must hold for every certificate that the checker accepts.
"""
import random
import sys
from fractions import Fraction as Fr

import magnet as M

sys.stdout.reconfigure(encoding="utf-8")
fails = 0
rnd = random.Random(20261123)


def report(tag, ok, extra=""):
    global fails
    fails += not ok
    print(f"  {tag:66s} {'ok' if ok else '**** FAILED ****'} {extra}", flush=True)


# ---------------------------------------------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------------------------------------------
print("A. certificates with known values")
half = Fr(1, 2)
# exchange only
mod = M.Model(A=half)
c0 = M.constant_certificate(mod, 0, [0])
c4 = M.constant_certificate(mod, 0, [2])
report("exchange only, u = 0: accepted; bound 4 pi |Q| in degree -1", c0.check() and c0.bound(-1) == (4, 0), f"bound = {c0.bound(-1)[0]} pi")
report("exchange only, u = 4A: accepted; bound 4 pi Q in degree +2", c4.check() and c4.bound(2) == (8, 0), f"bound = {c4.bound(2)[0]} pi")
report("   pointwise inequality (exact, random points)", pointwise_ok(c0)[0] and pointwise_ok(c4)[0])
# the critical family
for r in (Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(1)):
    mod = M.Model(A=half, kappa=r, K=-half, h=1)
    c = M.constant_certificate(mod, r, [2 * r * r])
    ok = c.check() and c.bound(-1) == (4 * (1 - 2 * r * r), 0)
    report(f"critical family r = {r}: bound 4 pi (1 - 2 r^2) in degree -1", ok and pointwise_ok(c)[0], f"bound = {c.bound(-1)[0]} pi")
    # nothing better: lowering u by one part in 10^6 must be rejected (the bound is the exact infimum)
    c_bad = M.constant_certificate(mod, r, [2 * r * r * (1 - Fr(1, 10 ** 6))])
    report(f"   u lowered by 1e-6 is rejected", not c_bad.check())
# the simple bound for the Zeeman model
for h in (Fr(2), Fr(5), Fr(20)):
    mod = M.Model(A=half, kappa=1, h=h)
    c = M.constant_certificate(mod, 1, [1 / h, -1 / h])
    ok = c.check() and c.bound(-1) == (4 * (1 - 1 / h), 0)
    report(f"Zeeman model h = {h}: simple bound 4 pi (1 - 1/h)", ok and pointwise_ok(c)[0], f"bound = {c.bound(-1)[0]} pi")
# the wall
mod = M.Model(A=1, K=1)
wb = M.wall_bound(mod, [-1, 1], [[0, 2]], 3)
report("wall, A = K = 1, width 3: bound 4 sqrt(AK) w = 12", wb == 12, f"bound = {wb}")
report("   a steeper s is rejected", M.wall_bound(mod, [-1, 1], [[0, Fr(2000001, 1000000)]], 3) is None)

print("B. falsification: no accepted certificate may exceed a known infimum")
# critical family: the infimum on the plane, 4 pi (1 - 2 r^2), bounds from below the infimum on every domain that is
# not the plane only from below; a certificate with Gamma = 0 holds on every domain, hence on larger and larger
# squares, hence must not exceed the infimum on the plane.
accepted = 0
worst = None
for trial in range(400):
    r = Fr(rnd.randint(1, 8), 8)
    mod = M.Model(A=half, kappa=r, K=-half, h=1)
    n = rnd.choice([1, 2, 3])
    mesh = sorted({Fr(-1), Fr(1)} | {Fr(rnd.randint(-7, 7), 8) for _ in range(n - 1)})
    amp = rnd.choice([Fr(0), Fr(1, 1000), Fr(1, 50), Fr(1, 5)])
    vals = [r + amp * Fr(rnd.randint(-10, 10), 10) for _ in mesh]
    slopes = [amp * Fr(rnd.randint(-10, 10), 10) for _ in mesh]
    us = [[2 * r * r * (1 + rnd.choice([Fr(-1, 100), Fr(0), Fr(1, 100), Fr(1, 3)])) + amp * Fr(rnd.randint(-5, 5), 10),
           amp * Fr(rnd.randint(-5, 5), 10)] for _ in range(len(mesh) - 1)]
    c = M.Certificate.from_hermite(mod, mesh, vals, slopes, us)
    if c.check():
        accepted += 1
        margin = 4 * (1 - 2 * r * r) - c.bound(-1)[0]            # must be >= 0
        worst = margin if worst is None else min(worst, margin)
        ok_pt, _ = pointwise_ok(c, 60)
        if margin < 0 or not ok_pt:
            report(f"trial {trial}", False, f"margin {float(margin)}")
report(f"critical family, random certificates: {accepted} accepted of 400, none above the exact infimum", worst is not None and worst >= 0,
       f"smallest margin {float(worst) if worst is not None else None} pi")
# exchange only: nothing above 8 pi A |Q|
accepted = 0
worst = None
for trial in range(300):
    mod = M.Model(A=half)
    amp = rnd.choice([Fr(0), Fr(1, 1000), Fr(1, 20)])
    mesh = [Fr(-1), Fr(rnd.randint(-3, 3), 4), Fr(1)]
    vals = [amp * Fr(rnd.randint(-10, 10), 10) for _ in mesh]
    slopes = [amp * Fr(rnd.randint(-10, 10), 10) for _ in mesh]
    us = [[Fr(rnd.choice([0, 0, 1, 2, 2])) + amp * Fr(rnd.randint(-3, 3), 10)] for _ in range(2)]
    c = M.Certificate.from_hermite(mod, mesh, vals, slopes, us)
    if c.check():
        accepted += 1
        for Q in (-2, -1, 1, 3):
            margin = 4 * abs(Q) - c.bound(Q)[0]
            worst = margin if worst is None else min(worst, margin)
report(f"exchange only, random certificates: {accepted} accepted of 300, none above 8 pi A |Q|", worst is not None and worst >= 0)

print("C. the pointwise inequality for certificates with a positive Gamma (random fields)")
accepted = 0
allok = True
for trial in range(200):
    mod = M.Model(A=half, kappa=Fr(rnd.randint(0, 6), 4), K=Fr(rnd.randint(0, 4), 4), h=Fr(rnd.randint(1, 8), 2))
    mesh = [Fr(-1), Fr(rnd.randint(-3, 3), 4), Fr(1)]
    vals = [Fr(rnd.randint(-20, 20), 10) for _ in mesh]
    slopes = [Fr(rnd.randint(-20, 20), 10) for _ in mesh]
    us = [[Fr(rnd.randint(1, 19), 10), Fr(rnd.randint(-5, 5), 100)] for _ in range(2)]
    for Gamma in (Fr(1, 10), Fr(1), Fr(10), Fr(100), Fr(1000)):
        c = M.Certificate.from_hermite(mod, mesh, vals, slopes, us, Gamma)
        if c.check():
            accepted += 1
            ok_pt, w = pointwise_ok(c, 80)
            allok &= ok_pt
            break
report(f"random fields: {accepted} of 200 accepted with some Gamma; pointwise inequality exact at random points", allok and accepted > 20)

print("\nALL TESTS PASSED" if fails == 0 else f"\nFAILURES: {fails}")
sys.exit(1 if fails else 0)
