(* ndsolve_compare.wls -- stiffness of the same body from NDSolve's finite element solver, set against the certified range.

   Body: block (-1,1)^3, Young's modulus 1, Poisson ratio 3/10, spherical cavity of radius 1/2 at the centre.
   One octant is solved, with mirror conditions on x = 0, y = 0, z = 0 (normal displacement zero).
   Grip on x = 1, moved by 1/2:   clamped: u = (1/2, 0, 0);   frictionless: only u_x = 1/2.
   Faces y = 1, z = 1 and the sphere are free.
   K = 8 * (integral over the octant of strain : C : strain), which is the full-body stiffness for unit grip separation.

   NOT executed by the author of the certificate (no Wolfram kernel was available there). Fix syntax if needed,
   but do not change the body, the material or the boundary conditions. *)

Needs["NDSolve`FEM`"];
nu = 3/10; lam = nu/((1 + nu) (1 - 2 nu)); mu = 1/(2 (1 + nu));
region = RegionDifference[Cuboid[{0, 0, 0}, {1, 1, 1}], Ball[{0, 0, 0}, 1/2]];
certified = <|"clamped" -> {1.817322, 1.833487}, "frictionless" -> {1.745105, 1.758277}|>;

stiffness[grip_, maxCell_] := Module[{mesh, ops, bcs, uif, vif, wif, exx, eyy, ezz, gxy, gyz, gxz, dens, k, t},
  mesh = ToElementMesh[region, MaxCellMeasure -> maxCell, "MeshOrder" -> 2];
  ops = SolidMechanicsPDEComponent[{{u[x, y, z], v[x, y, z], w[x, y, z]}, {x, y, z}},
     <|"YoungModulus" -> 1, "PoissonRatio" -> nu|>];
  bcs = Join[
    {DirichletCondition[u[x, y, z] == 0, x == 0], DirichletCondition[v[x, y, z] == 0, y == 0],
     DirichletCondition[w[x, y, z] == 0, z == 0], DirichletCondition[u[x, y, z] == 1/2, x == 1]},
    If[grip === "clamped", {DirichletCondition[v[x, y, z] == 0, x == 1], DirichletCondition[w[x, y, z] == 0, x == 1]}, {}]];
  {t, {uif, vif, wif}} = AbsoluteTiming[NDSolveValue[{ops == {0, 0, 0}, bcs}, {u, v, w}, {x, y, z} \[Element] mesh]];
  exx = D[uif[x, y, z], x]; eyy = D[vif[x, y, z], y]; ezz = D[wif[x, y, z], z];
  gxy = D[uif[x, y, z], y] + D[vif[x, y, z], x];
  gyz = D[vif[x, y, z], z] + D[wif[x, y, z], y];
  gxz = D[uif[x, y, z], z] + D[wif[x, y, z], x];
  dens = lam (exx + eyy + ezz)^2 + 2 mu (exx^2 + eyy^2 + ezz^2) + mu (gxy^2 + gyz^2 + gxz^2);
  k = 8 NIntegrate[dens, {x, y, z} \[Element] mesh];
  {Length[mesh["Coordinates"]], t, k}];

verdict[k_, {lo_, hi_}] := Which[
  k > hi, "ABOVE the proven range by at least " <> ToString[NumberForm[100 (k - hi)/hi, {5, 3}]] <> " %",
  k < lo, "BELOW the proven range by at least " <> ToString[NumberForm[100 (lo - k)/lo, {5, 3}]] <> " %",
  True, "inside the proven range"];

Do[
  Print[grip, " grips; proven: ", certified[grip][[1]], " <= K <= ", certified[grip][[2]]];
  Do[
    Module[{r = stiffness[grip, mc]},
      Print["   max cell volume ", ScientificForm[mc, 2], "   nodes ", r[[1]], "   solve time ", NumberForm[r[[2]], {6, 2}],
        " s   K = ", NumberForm[r[[3]], {8, 6}], "   ", verdict[r[[3]], certified[grip]]]],
    {mc, {1.*^-2, 1.*^-3, 2.*^-4, 5.*^-5, 1.*^-5}}],
  {grip, {"clamped", "frictionless"}}];
