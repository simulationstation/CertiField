(* repro_ndsolve_gradient.wls -- minimal reproducer: NIntegrate over an ElementMesh region of a function that is
   piecewise constant on the mesh elements (the squared gradient of a linear finite element field).

   wolframscript -file benchpack/harness/repro_ndsolve_gradient.wls

   Used to explain why K computed as NIntegrate of the strain energy density of NDSolve's InterpolatingFunction
   differs from the finite element energy a(u_h, u_h) on the 2D harness meshes (benchpack/results/VERDICTS.md).
   FLOAT / DIAGNOSTIC; checked by the agent only.

   Part 1. Two linear triangles on the unit square, nodal values u(0,0) = 0, u(1,0) = 1, u(1,1) = 0, u(0,1) = 0.
   The interpolant is u = x - y on the lower triangle and u = 0 on the upper one, so du/dx = 1 at (0.7, 0.2),
   du/dx = 0 at (0.2, 0.7), and the integral of |grad u|^2 over the square is exactly 1.
   Part 2. The same on a mesh made by ToElementMesh[Rectangle[]] and on the same mesh rebuilt from its coordinates and
   elements, with pseudo-random nodal values; exact value = sum over the triangles of area * |grad u|^2. *)

Needs["NDSolve`FEM`"];
mesh = ToElementMesh["Coordinates" -> {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}},
   "MeshElements" -> {TriangleElement[{{1, 2, 3}, {1, 3, 4}}]}];
if = ElementMeshInterpolation[{mesh}, {0., 1., 0., 0.}];
Print["Part 1"];
Print["u(0.7, 0.2)            = ", if[0.7, 0.2], "   (piecewise linear interpolant: 0.5)"];
Print["du/dx at (0.7, 0.2)    = ", Derivative[1, 0][if][0.7, 0.2], "   (element gradient: 1)"];
Print["du/dx at (0.2, 0.7)    = ", Derivative[1, 0][if][0.2, 0.7], "   (element gradient: 0)"];
Do[Print["du/dx at (0.5 + d, 0.5 - d), inside the lower triangle, d = ", d, ":  ",
   Derivative[1, 0][if][0.5 + d, 0.5 - d], "   (element gradient: 1)"], {d, {1.*^-2, 1.*^-3, 1.*^-4, 1.*^-6}}];
g = Derivative[1, 0][if][x, y]^2 + Derivative[0, 1][if][x, y]^2;
Print["integral |grad u|^2, NIntegrate default               = ", NIntegrate[g, {x, y} \[Element] mesh],
  "   (exact: 1)"];
Print["integral |grad u|^2, NIntegrate PrecisionGoal -> 10    = ",
  NIntegrate[g, {x, y} \[Element] mesh, PrecisionGoal -> 10]];
Print["integral |grad u|^2, NIntegrate Method FiniteElement   = ",
  Quiet[Check[NIntegrate[g, {x, y} \[Element] mesh, Method -> "FiniteElement"], "method not accepted"]]];
Print["integral |grad u|^2, sum over the two triangles        = ",
  NIntegrate[g, {x, y} \[Element] Triangle[{{0, 0}, {1, 0}, {1, 1}}]] +
   NIntegrate[g, {x, y} \[Element] Triangle[{{0, 0}, {1, 1}, {0, 1}}]]];
Print["smooth integrand x^2 + y^2 over the mesh (exact 2/3)   = ", NIntegrate[x^2 + y^2, {x, y} \[Element] mesh]];

Print["Part 2"];
m1 = ToElementMesh[Rectangle[], MaxCellMeasure -> 0.01, "MeshOrder" -> 1, "MeshElementType" -> TriangleElement];
crd = m1["Coordinates"]; els = m1["MeshElements"][[1, 1]];
SeedRandom[1];
vals = RandomReal[{-1, 1}, Length[crd]];
f = ElementMeshInterpolation[{m1}, vals];
grads = Map[Module[{p = crd[[#]], v = vals[[#]], J}, J = {p[[2]] - p[[1]], p[[3]] - p[[1]]};
     LinearSolve[J, {v[[2]] - v[[1]], v[[3]] - v[[1]]}]] &, els];
Print["triangles: ", Length[els], ";  max |gradient of the interpolation at a centroid - element gradient| = ",
  Max[Table[Norm[{Derivative[1, 0][f] @@ Mean[crd[[els[[k]]]]], Derivative[0, 1][f] @@ Mean[crd[[els[[k]]]]]} -
      grads[[k]]], {k, Length[els]}]]];
ex = Total[Table[Module[{p = crd[[els[[k]]]], J}, J = {p[[2]] - p[[1]], p[[3]] - p[[1]]};
     Abs[Det[J]]/2 grads[[k]] . grads[[k]]], {k, Length[els]}]];
gg = Derivative[1, 0][f][x, y]^2 + Derivative[0, 1][f][x, y]^2;
ni = NIntegrate[gg, {x, y} \[Element] m1];
ni2 = Total[Table[NIntegrate[gg, {x, y} \[Element] Triangle[crd[[els[[k]]]]]], {k, Length[els]}]];
Print["sum of area * |element gradient|^2          = ", ex];
Print["NIntegrate over the ElementMesh             = ", ni, "   relative difference ", (ni - ex)/ex];
Print["sum of NIntegrate over the single triangles = ", ni2, "   relative difference ", (ni2 - ex)/ex];
Print[$Version];
