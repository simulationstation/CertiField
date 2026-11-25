"""The cases run by the harness: geometry, grip conditions and the two quantities, as in benchpack/problems/*.md.

Dirichlet data are lists of (plane, component, value): plane = (axis, coordinate); every degree of freedom of the
given component located on that plane is prescribed.  Components not listed are free (zero traction, or zero
shear for a frictionless grip, where only the normal component is listed).

Q1 = K is computed from the reaction on the moving grip:  K = scale * (sum of the reactions of the listed component
on the grip plane) / (prescribed grip displacement of the model), which for a conforming method equals
a(u_h, u_h) / delta^2 of the full body (both are recorded).
"""
import math

E = 1.0
NU = 0.3
TH1 = math.pi / 2 - 0.1

P1_PATCH_AREA = 0.5 * 0.2 * (0.55 ** 2 - 0.5 ** 2)       # |P| of the full patch = 0.00525


def in_p1_patch(c):
    """c: (2, n) element centroids.  The half patch of the quarter model (mesh conforms to its boundary)."""
    r = (c[0] ** 2 + c[1] ** 2) ** 0.5
    th = [math.atan2(y, x) for x, y in zip(c[0], c[1])]
    import numpy as np
    return (r > 0.5) & (r < 0.55) & (np.array(th) > TH1)


def in_p2_patch(c):
    return (c[0] > 0) & (c[0] < 0.1) & (c[1] > 0) & (c[1] < 0.1)


# The grip moves by d_model in the model; delta = 1 for the full body.
CASES = {
    # sanity case: P1 without the hole, frictionless grips; exact K = 1, Q2 = 0.5
    "P1solid": dict(dim=2, mesh="plate", hole=False,
                    bc=[((0, 0.0), 0, 0.0), ((1, 0.0), 1, 0.0), ((0, 1.0), 0, 0.5)],
                    grip=((0, 1.0), 0), d_model=0.5, kscale=2.0, energy_scale=4.0, q2="p1"),
    "P1": dict(dim=2, mesh="plate", hole=True,
               bc=[((0, 0.0), 0, 0.0), ((1, 0.0), 1, 0.0), ((0, 1.0), 0, 0.5)],
               grip=((0, 1.0), 0), d_model=0.5, kscale=2.0, energy_scale=4.0, q2="p1"),
    "P1c": dict(dim=2, mesh="plate", hole=True, grade=True,
                bc=[((0, 0.0), 0, 0.0), ((1, 0.0), 1, 0.0), ((0, 1.0), 0, 0.5), ((0, 1.0), 1, 0.0)],
                grip=((0, 1.0), 0), d_model=0.5, kscale=2.0, energy_scale=4.0, q2="p1"),
    # P2: full model, left edge clamped, right edge moved by (0, -1); K = -F_y
    "P2": dict(dim=2, mesh="lbracket",
               bc=[((0, -1.0), 0, 0.0), ((0, -1.0), 1, 0.0), ((0, 1.0), 0, 0.0), ((0, 1.0), 1, -1.0)],
               grip=((0, 1.0), 1), d_model=-1.0, kscale=1.0, energy_scale=1.0, q2="p2"),
    # P3: octant, mirror planes x, y, z = 0; grip x = 1 moved by 1/2
    "P3c": dict(dim=3, mesh="block", grade=True,
                bc=[((0, 0.0), 0, 0.0), ((1, 0.0), 1, 0.0), ((2, 0.0), 2, 0.0), ((0, 1.0), 0, 0.5),
                    ((0, 1.0), 1, 0.0), ((0, 1.0), 2, 0.0)],
                grip=((0, 1.0), 0), d_model=0.5, kscale=4.0, energy_scale=8.0, q2=None),
    "P3f": dict(dim=3, mesh="block",
                bc=[((0, 0.0), 0, 0.0), ((1, 0.0), 1, 0.0), ((2, 0.0), 2, 0.0), ((0, 1.0), 0, 0.5)],
                grip=((0, 1.0), 0), d_model=0.5, kscale=4.0, energy_scale=8.0, q2=None),
}


def make_mesh(case, h, order, curved):
    import meshes
    c = CASES[case]
    if c["mesh"] == "plate":
        return meshes.plate_quarter(h, order=order, curved=curved, hole=c["hole"], grade_grip_corner=c.get("grade", False))
    if c["mesh"] == "lbracket":
        return meshes.lbracket(h, order=order)
    if c["mesh"] == "block":
        return meshes.block_octant(h, order=order, curved=curved, grade_grip_edges=c.get("grade", False))
    raise KeyError(case)


# mesh levels (nominal size h) per problem
LEVELS = {
    "P1solid": [0.2, 0.1, 0.05, 0.025],
    "P1": [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625],
    "P1c": [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625],
    "P2": [0.4, 0.2, 0.1, 0.05, 0.025, 0.0125],
    "P3c": [0.4, 0.2, 0.1, 0.085, 0.07, 0.05],
    "P3f": [0.4, 0.2, 0.1, 0.085, 0.07, 0.05],
}
# quadratic elements: fewer levels where the finest linear level would exceed the memory of a direct solve
LEVELS_QUAD = {"P3c": [0.4, 0.2, 0.1, 0.085, 0.07]}
