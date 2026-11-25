(* ndsolve_P3.wls -- copy of Reference/Curved3D/ndsolve_compare.wls for the benchmark harness.

   Body, material, boundary conditions and the stiffness formula are unchanged from the original.  Changes:
   - the certified ranges are those of benchpack/problems/P3.md (grid 24), to 8 decimals, rounded outward;
   - the mesh measure (volume of the meshed octant) and the run data are also written to
     benchpack/results/ndsolve-toelementmesh_P3.csv;
   - numbers are printed with CForm (NumberForm is not rendered by wolframscript -file);
   - added: K_system = 8 U.S.U from the assembled finite element system of the same problem, i.e. a(u_h, u_h)
     (the NIntegrate value of the original differs from it; see benchpack/results/VERDICTS.md);
   - loop order: levels outer, grips inner (the finest level stops the kernel on the harness machine).
   Run from the repository root:   wolframscript -file benchpack/harness/ndsolve_P3.wls

   Body: block (-1,1)^3, Young's modulus 1, Poisson ratio 3/10, spherical cavity of radius 1/2 at the centre.
   One octant is solved, with mirror conditions on x = 0, y = 0, z = 0 (normal displacement zero).
   Grip on x = 1, moved by 1/2:   clamped: u = (1/2, 0, 0);   frictionless: only u_x = 1/2.
   Faces y = 1, z = 1 and the sphere are free.
   K = 8 * (integral over the octant of strain : C : strain), which is the full-body stiffness for unit grip separation.
   Values are FLOAT (no error control); checked by the agent only. *)

Needs["NDSolve`FEM`"];
nu = 3/10; lam = nu/((1 + nu) (1 - 2 nu)); mu = 1/(2 (1 + nu));
region = RegionDifference[Cuboid[{0, 0, 0}, {1, 1, 1}], Ball[{0, 0, 0}, 1/2]];
certified = <|"clamped" -> {1.81732246, 1.83348613}, "frictionless" -> {1.74510549, 1.75827651}|>;
exactVolume = N[1 - Pi/48, 16];

stiffness[grip_, maxCell_] := Module[{mesh, ops, bcs, uif, vif, wif, exx, eyy, ezz, gxy, gyz, gxz, dens, k, t, tk, vol,
    ksys},
  mesh = ToElementMesh[region, MaxCellMeasure -> maxCell, "MeshOrder" -> 2];
  vol = Total[Flatten[mesh["MeshElementMeasure"]]];
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
  {tk, k} = AbsoluteTiming[8 NIntegrate[dens, {x, y, z} \[Element] mesh]];
  (* added for the harness: the finite element energy a(u_h, u_h) from the assembled system of the same problem *)
  Module[{state, femdata, sd, dpde, dbcs, load, stiff, damp, mass, S0, U},
    {state} = NDSolve`ProcessEquations[{ops == {0, 0, 0}, bcs}, {u, v, w}, {x, y, z} \[Element] mesh];
    femdata = state["FiniteElementData"]; sd = state["SolutionData"][[1]];
    dpde = DiscretizePDE[femdata["PDECoefficientData"], femdata["FEMMethodData"], sd];
    dbcs = DiscretizeBoundaryConditions[femdata["BoundaryConditionData"], femdata["FEMMethodData"], sd];
    {load, stiff, damp, mass} = dpde["SystemMatrices"]; S0 = stiff;
    DeployBoundaryConditions[{load, stiff}, dbcs];
    U = Flatten[LinearSolve[stiff, load]];
    ksys = 8 U . S0 . U];
  {Length[mesh["Coordinates"]], t, k, Length[mesh["MeshElements"][[1, 1]]], vol, tk, ksys}];

str[x_] := ToString[CForm[N[x]]];
verdict[k_, {lo_, hi_}] := Which[
  k > hi, "ABOVE the proven range by at least " <> str[Round[100 (k - hi)/hi, 0.001]] <> " %",
  k < lo, "BELOW the proven range by at least " <> str[Round[100 (lo - k)/lo, 0.001]] <> " %",
  True, "inside the proven range"];

rows = {{"code", "case", "element", "maxcell", "ndof", "nelem", "time_s", "integrate_time_s", "K_nintegrate",
   "K_system", "mesh_volume", "exact_volume"}};
csv = FileNameJoin[{Directory[], "benchpack", "results", "ndsolve-toelementmesh_P3.csv"}];
(* loop order changed (levels outer, grips inner): the finest level, 1e-5, stopped the kernel on the harness machine
   (about 750 000 unknowns), and it should not keep the frictionless runs from being done *)
Do[
  Print[grip, " grips; proven: ", certified[grip][[1]], " <= K <= ", certified[grip][[2]]];
    Module[{r = stiffness[grip, mc]},
      Print["   max cell volume ", str[mc], "   nodes ", r[[1]], "   solve time ", str[Round[r[[2]], 0.01]],
        " s   K = ", str[r[[3]]], "   ", verdict[r[[3]], certified[grip]], "   K_system = ", str[r[[7]]], "  ",
        verdict[r[[7]], certified[grip]],
        "   mesh volume - exact ", str[r[[5]] - exactVolume]];
      AppendTo[rows, {"NDSolve", If[grip === "clamped", "P3c", "P3f"], "quadratic tet (MeshOrder 2)", mc, 3 r[[1]],
        r[[4]], r[[2]], r[[6]], str[r[[3]]], str[r[[7]]], str[r[[5]]], str[exactVolume]}];
      Export[csv, Map[ToString, rows, {2}], "CSV"]],
  {mc, {1.*^-2, 1.*^-3, 2.*^-4, 5.*^-5, 1.*^-5}}, {grip, {"clamped", "frictionless"}}];
Print["Mathematica ", $Version];
