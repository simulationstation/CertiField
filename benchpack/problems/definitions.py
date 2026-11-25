"""Machine-readable definitions of the benchmark problems P1 and P2 (exact rationals).

This is the only file that the certificate builder (benchpack/builder) and the independent checker
(benchpack/checker) may share.  It states the problems of P1.md and P2.md in exact form, plus the few fixed
conventions that the certificates rely on (the symmetry reduction of P1 and the weight polynomial rho of P1).
It contains no numerical method.

Conventions used everywhere:
- Plane stress, unit thickness, Voigt order (xx, yy, xy) with engineering shear strain gamma_xy = du_x/dy + du_y/dx.
  sigma = C eps,  eps = S sigma  (S = C^-1),  energy density (1/2) eps.C.eps.
- Airy stress function phi: sigma_xx = d2phi/dy2, sigma_yy = d2phi/dx2, sigma_xy = -d2phi/dxdy.
- delta (the grip displacement) is 1 in every problem.
"""
from fractions import Fraction as Fr

# ---------------------------------------------------------------------------------------------------------------------
# material (common to P1 and P2)
E = Fr(1)
NU = Fr(3, 10)
DELTA = Fr(1)


def plane_stress_C(E=E, nu=NU):
    """3x3 stiffness matrix (Voigt xx, yy, xy with engineering shear): sigma = C eps"""
    c = E / (1 - nu * nu)
    return [[c, nu * c, Fr(0)], [nu * c, c, Fr(0)], [Fr(0), Fr(0), E / (2 * (1 + nu))]]


def plane_stress_S(E=E, nu=NU):
    """3x3 compliance matrix S = C^-1: eps = S sigma"""
    return [[1 / E, -nu / E, Fr(0)], [-nu / E, 1 / E, Fr(0)], [Fr(0), Fr(0), 2 * (1 + nu) / E]]


C = plane_stress_C()          # [[100/91, 30/91, 0], [30/91, 100/91, 0], [0, 0, 5/13]]
S = plane_stress_S()          # [[1, -3/10, 0], [-3/10, 1, 0], [0, 0, 13/5]]
C1111 = C[0][0]               # 100/91: sigma_xx produced by a unit eps_xx at zero eps_yy

# ---------------------------------------------------------------------------------------------------------------------
# P1: gripped plate with a circular hole
P1 = dict(
    name='P1',
    domain='[-1,1]^2 minus the open disc of radius a centred at (0,0)',
    box=((Fr(-1), Fr(1)), (Fr(-1), Fr(1))),
    hole_radius=Fr(1, 2),
    grips=dict(faces=('x=-1', 'x=1'), u_x=(-DELTA / 2, DELTA / 2),
               frictionless='u_x prescribed, shear traction zero (the default variant)',
               clamped='u = (+-delta/2, 0) (variant P1c)'),
    free=('y=-1', 'y=1', 'the hole boundary x^2 + y^2 = 1/4'),
    # Q1 = K = (total x-force on the grip x = 1)/delta = a(u,u)/delta^2  (a(u,v) = int eps(u).C.eps(v))
    # Q2 = (1/|P|) int_P sigma_xx dA, P = {(r cos t, r sin t): r_min <= r <= r_max, |t - pi/2| <= t_half}
    patch=dict(r_min=Fr(1, 2), r_max=Fr(11, 20), theta_centre='pi/2', theta_half=Fr(1, 10)),
    patch_area=Fr(21, 4000),      # (2 t_half / 2)(r_max^2 - r_min^2) = (1/10)(21/400)
    known_answers_without_hole=dict(K=Fr(1), Q2=Fr(1, 2)),
)

# Symmetry reduction used by the P1 certificates (allowed by P1.md).  The solution is symmetric under x -> -x and
# y -> -y; the vertical rigid translation (free for frictionless grips) is fixed by u_y = 0 on y = 0.
P1_QUARTER = dict(
    domain='Omega_q = {(x,y) in [0,1]^2 : x^2 + y^2 >= 1/4}',
    box=((Fr(0), Fr(1)), (Fr(0), Fr(1))),
    displacement_conditions=dict(x0='u_x = 0', y0='u_y = 0', x1='u_x = delta/2 (frictionless)',
                                 x1_clamped='u = (delta/2, 0)'),
    stress_conditions=dict(x0='sigma_xy = 0', y0='sigma_xy = 0', x1='sigma_xy = 0 (frictionless only)',
                           y1='sigma.n = 0', hole='sigma.n = 0'),
    # K = 4 a_q(u,u) where a_q is the energy form on Omega_q; minimal energy on Omega_q is K delta^2 / 8.
    K_from_quarter_energy=Fr(4),
    # half patch P_h = P with x >= 0 = {r_min <= r <= r_max, pi/2 - t_half <= t <= pi/2};  Q2 = (1/|P_h|) int_{P_h} sigma_xx
    half_patch_area=Fr(21, 8000),
    # Force of a stress field through the quarter grip, F_q = int_0^1 sigma_xx(1,y) dy; full grip force = 2 F_q.
)

# Weight polynomial of the P1 stress fields: rho(x,y) = sum RHO[(i,j)] x^i y^j
#   rho = 2x^2 - x^4 + y^4 + (3/2) y^2 - 7/16 = (x^2 + y^2 - 1/4)(7/4 + y^2 - x^2).
# On [0,1]^2 the second factor is >= 3/4 > 0, so on the box rho = 0 exactly on the hole boundary.
# d rho/dx = 4x(1 - x^2) vanishes on x = 0 and x = 1; d rho/dy = 4y^3 + 3y vanishes on y = 0.
RHO = {(2, 0): Fr(2), (4, 0): Fr(-1), (0, 4): Fr(1), (0, 2): Fr(3, 2), (0, 0): Fr(-7, 16)}

# ---------------------------------------------------------------------------------------------------------------------
# P2: L-bracket
P2 = dict(
    name='P2',
    domain='(-1,1)^2 minus [0,1) x (-1,0]; re-entrant corner at (0,0)',
    box=((Fr(-1), Fr(1)), (Fr(-1), Fr(1))),
    removed=((Fr(0), Fr(1)), (Fr(-1), Fr(0))),
    clamped='x = -1, -1 <= y <= 1: u = (0, 0)',
    grip='x = 1, 0 <= y <= 1: u = (0, -delta) (rigid grip, clamped to the material)',
    free=('y = 1, -1 <= x <= 1', 'y = -1, -1 <= x <= 0', 'x = 0, -1 <= y <= 0', 'y = 0, 0 <= x <= 1'),
    # Q1 = K = (total vertical force on the right grip)/delta = a(u,u)/delta^2
    # Q2 = (1/|S|) int_S sigma_xx dA, S = [0, 1/10] x [0, 1/10]
    patch=((Fr(0), Fr(1, 10)), (Fr(0), Fr(1, 10))),
    patch_area=Fr(1, 100),
)


def p1_in_quarter(x, y):
    """exact test: (x, y) (Fractions) lies in the closed quarter domain Omega_q"""
    x, y = Fr(x), Fr(y)
    return 0 <= x <= 1 and 0 <= y <= 1 and x * x + y * y >= Fr(1, 4)


def p2_in_domain(x, y):
    """exact test: (x, y) lies in the closure of the L-bracket"""
    x, y = Fr(x), Fr(y)
    return -1 <= x <= 1 and -1 <= y <= 1 and not (x > 0 and y < 0)
