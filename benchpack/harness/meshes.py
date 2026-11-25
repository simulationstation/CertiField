"""gmsh meshes for the benchmark problems P1, P2, P3 (definitions in benchpack/problems/).

Every function returns a dict
    points : (n, d) float array of node coordinates (d = 2 or 3)
    cells  : (m, k) int array, VTK / Abaqus node order
             tri3: 0 1 2 | tri6: 0 1 2, (0,1) (1,2) (2,0)
             tet4: 0 1 2 3 | tet10: 0 1 2 3, (0,1) (1,2) (2,0) (0,3) (1,3) (2,3)
    order  : 1 or 2
    curved : True if the mid-edge nodes of boundary edges lie on the curved boundary,
             False if they sit at the chord midpoints (polygonal / polyhedral boundary)
    h      : the nominal mesh size of the level

The quarter plate (P1) and the octant (P3) use the symmetry planes of the problem definitions.
The patch of Q2 is a separate surface in P1 and P2, so the mesh conforms to the patch boundary.
"""
import math
import numpy as np
import gmsh

TH1 = math.pi / 2 - 0.1          # the patch of P1 spans pi/2 - 0.1 <= theta <= pi/2 in the quarter model


def _init(name):
    if gmsh.isInitialized():
        gmsh.finalize()
    gmsh.initialize(["-noenv"])
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.model.add(name)


def _extract(dim, order, curved, h):
    gmsh.model.mesh.generate(dim)
    if order == 2:
        # gmsh places the new nodes on the CAD geometry (curved boundary); for the straight variant
        # they are moved back to the chord midpoints below.
        gmsh.model.mesh.setOrder(2)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    pts = coords.reshape(-1, 3)
    etype = {(2, 1): 2, (2, 2): 9, (3, 1): 4, (3, 2): 11}[(dim, order)]
    _, enodes = gmsh.model.mesh.getElementsByType(etype)
    nn = {2: 3, 9: 6, 4: 4, 11: 10}[etype]
    cells_tag = np.asarray(enodes, dtype=np.int64).reshape(-1, nn)
    # renumber to 0..n-1 over the nodes actually used
    used = np.unique(cells_tag)
    tag2idx = {int(t): i for i, t in enumerate(tags)}
    idx = np.array([tag2idx[int(t)] for t in used])
    pts = pts[idx, :dim].copy()
    remap = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    cells = remap[cells_tag]
    if etype == 11:
        cells = cells[:, [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]]      # gmsh tet10 -> VTK / Abaqus order
    gmsh.finalize()
    if order == 2:
        edges = [(0, 1), (1, 2), (2, 0)] if dim == 2 else [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
        nv = dim + 1
        # check the node order: every mid-edge node lies close to its edge midpoint
        for j, (a, b) in enumerate(edges):
            mid = 0.5 * (pts[cells[:, a]] + pts[cells[:, b]])
            L = np.linalg.norm(pts[cells[:, a]] - pts[cells[:, b]], axis=1)
            dev = np.linalg.norm(pts[cells[:, nv + j]] - mid, axis=1)
            if np.any(dev > 0.3 * L):
                raise RuntimeError("mid-edge node order check failed")
            if not curved:
                pts[cells[:, nv + j]] = mid
    return {"points": pts, "cells": cells, "order": order, "curved": bool(curved and order == 2), "h": h}


def plate_quarter(h, order=1, curved=True, hole=True, hole_ratio=0.25, grade_grip_corner=False):
    """P1 quarter model [0,1]^2 minus the disc r < 1/2 (or the solid square if hole=False).
    Size h on the outer boundary, h*hole_ratio on the hole and the patch.  With grade_grip_corner (clamped grips),
    the size is also limited by h * r^0.6, r = distance to the grip corner (1, 1), where the stress is singular."""
    _init("P1")
    g = gmsh.model.geo
    hh = h * hole_ratio
    c1, s1 = math.cos(TH1), math.sin(TH1)
    pO = g.addPoint(0, 0, 0, hh)
    pB = g.addPoint(1, 0, 0, h)
    pC = g.addPoint(1, 1, 0, h)
    pD = g.addPoint(0, 1, 0, h)
    pE = g.addPoint(0, 0.55, 0, hh)
    pF = g.addPoint(0, 0.5, 0, hh)
    pG = g.addPoint(0.5 * c1, 0.5 * s1, 0, hh)
    pH = g.addPoint(0.55 * c1, 0.55 * s1, 0, hh)
    lEF = g.addLine(pE, pF)
    aFG = g.addCircleArc(pF, pO, pG)
    lGH = g.addLine(pG, pH)
    aHE = g.addCircleArc(pH, pO, pE)
    lBC = g.addLine(pB, pC)
    lCD = g.addLine(pC, pD)
    lDE = g.addLine(pD, pE)
    if hole:
        pA = g.addPoint(0.5, 0, 0, hh)
        aGA = g.addCircleArc(pG, pO, pA)
        lAB = g.addLine(pA, pB)
        rest = g.addCurveLoop([lAB, lBC, lCD, lDE, -aHE, -lGH, aGA])
    else:
        pA = g.addPoint(0.5, 0, 0, hh)
        aGA = g.addCircleArc(pG, pO, pA)          # interior line, keeps the two meshes alike
        lOA = g.addLine(pO, pA)
        lAB = g.addLine(pA, pB)
        lFO = g.addLine(pF, pO)
        rest = g.addCurveLoop([lAB, lBC, lCD, lDE, -aHE, -lGH, aGA])
        inner = g.addCurveLoop([lOA, -aGA, -aFG, lFO])
    patch = g.addCurveLoop([lEF, aFG, lGH, aHE])
    g.addPlaneSurface([patch])
    g.addPlaneSurface([rest])
    if not hole:
        g.addPlaneSurface([inner])
    g.synchronize()
    if grade_grip_corner:
        gmsh.model.mesh.setSizeCallback(
            lambda dim, tag, x, y, z, lc: min(lc, h * max(math.hypot(1 - x, 1 - y), 1e-7) ** 0.6))
    return _extract(2, order, curved, h)


def lbracket(h, order=1, beta=0.75, beta_grip=0.6, rmin=1e-7):
    """P2 L-bracket (-1,1)^2 minus [0,1) x (-1,0], with the square [0,0.1]^2 as its own surface.
    Mesh size min(h r0^beta, h r_k^beta_grip, h): r0 = distance to the re-entrant corner, r_k = distance to the four
    corners where a clamped edge meets a free edge, (-1,-1), (-1,1), (1,0), (1,1) (graded refinement)."""
    _init("P2")
    g = gmsh.model.geo
    P = {k: g.addPoint(x, y, 0, h) for k, (x, y) in {
        "p0": (-1, -1), "p1": (0, -1), "c": (0, 0), "q1": (0.1, 0), "p3": (1, 0), "p4": (1, 1),
        "p5": (-1, 1), "q2": (0.1, 0.1), "q3": (0, 0.1)}.items()}
    L = lambda a, b: g.addLine(P[a], P[b])
    l01, l1c, lcq1, lq13, l34, l45, l50 = L("p0", "p1"), L("p1", "c"), L("c", "q1"), L("q1", "p3"), \
        L("p3", "p4"), L("p4", "p5"), L("p5", "p0")
    lq12, lq23, lq3c = L("q1", "q2"), L("q2", "q3"), L("q3", "c")
    sq = g.addCurveLoop([lcq1, lq12, lq23, lq3c])
    rest = g.addCurveLoop([l01, l1c, -lq3c, -lq23, -lq12, lq13, l34, l45, l50])
    g.addPlaneSurface([sq])
    g.addPlaneSurface([rest])
    g.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    corners = [(-1, -1), (-1, 1), (1, 0), (1, 1)]

    def size(dim, tag, x, y, z, lc):
        s0 = max(math.hypot(x, y), rmin) ** beta
        s1 = min(max(math.hypot(x - a, y - b), rmin) for a, b in corners) ** beta_grip
        return h * min(1.0, s0, s1)
    gmsh.model.mesh.setSizeCallback(size)
    return _extract(2, order, False, h)


def block_octant(h, order=1, curved=True, cavity_ratio=0.35, grade_grip_edges=False):
    """P3 octant [0,1]^3 minus the ball r < 1/2.  Size h*cavity_ratio on the sphere, growing linearly to h.
    With grade_grip_edges (clamped grips) the size is also limited by h * max(d, 1/64)^0.5, d = distance to the edges where the
    grip x = 1 meets the free faces y = 1 and z = 1 (singular stresses there)."""
    _init("P3")
    o = gmsh.model.occ
    box = o.addBox(0, 0, 0, 1, 1, 1)
    ball = o.addSphere(0, 0, 0, 0.5)
    o.cut([(3, box)], [(3, ball)])
    o.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    s = cavity_ratio

    def size(dim, tag, x, y, z, lc):
        r = math.sqrt(x * x + y * y + z * z)
        val = h * (s + (1 - s) * min(1.0, max(0.0, (r - 0.5) / 0.5)))
        if grade_grip_edges:
            d = min(math.hypot(1 - x, 1 - y), math.hypot(1 - x, 1 - z))
            val = min(val, h * max(d, 1.0 / 64) ** 0.5)      # at most 8 times finer than h
        return val
    gmsh.model.mesh.setSizeCallback(size)
    return _extract(3, order, curved, h)


def mesh_measure(m):
    """Area (2D) or volume (3D) of the meshed domain, exact for the element geometry (quadratic maps included)."""
    p, c = m["points"], m["cells"]
    d = p.shape[1]
    if d == 2:
        if m["order"] == 1:
            a, b, cc = p[c[:, 0]], p[c[:, 1]], p[c[:, 2]]
            return float(np.sum(0.5 * np.abs((b - a)[:, 0] * (cc - a)[:, 1] - (b - a)[:, 1] * (cc - a)[:, 0])))
        # quadratic triangle: area via the boundary integral of x dy over the 3 parabolic edges (exact)
        tot = 0.0
        for (i, j, k) in [(0, 1, 3), (1, 2, 4), (2, 0, 5)]:
            A, B, M = p[c[:, i]], p[c[:, j]], p[c[:, k]]
            # x(s) = A(1-s)(1-2s) + 4M s(1-s) + B s(2s-1); integrate x y' ds exactly with 3-point Gauss (degree 3)
            for s, w in zip([0.5 - math.sqrt(0.15), 0.5, 0.5 + math.sqrt(0.15)], [5 / 18, 8 / 18, 5 / 18]):
                N = np.array([(1 - s) * (1 - 2 * s), 4 * s * (1 - s), s * (2 * s - 1)])
                dN = np.array([4 * s - 3, 4 - 8 * s, 4 * s - 1])
                X = N[0] * A + N[1] * M + N[2] * B
                dX = dN[0] * A + dN[1] * M + dN[2] * B
                tot += w * np.sum(X[:, 0] * dX[:, 1])
        return float(abs(tot)) if tot > 0 else float(-tot)
    raise NotImplementedError


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    for curved in (False, True):
        m = plate_quarter(0.2, order=2, curved=curved)
        print("P1 quarter, quadratic, curved =", curved, "nodes", len(m["points"]), "cells", len(m["cells"]),
              "area", mesh_measure(m), "exact", 1 - math.pi / 16)
    m = plate_quarter(0.2, order=1)
    print("P1 quarter, linear: area", mesh_measure(m))
    m = lbracket(0.2, order=2)
    print("P2 nodes", len(m["points"]), "cells", len(m["cells"]))
    m = block_octant(0.2, order=2)
    print("P3 nodes", len(m["points"]), "cells", len(m["cells"]))
