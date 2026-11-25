(* ndsolve_fem.wls -- NDSolve's finite element solver on the harness meshes (gmsh), driven by run_ndsolve.py.

   wolframscript -file benchpack/harness/ndsolve_fem.wls <jobs.wl>

   <jobs.wl> holds a list of associations, one per run, written by run_ndsolve.py:
     "name", "dim", "coords" (CSV), "cells" (CSV, 1-based, Mathematica = Abaqus node order), "patchcoords",
     "patchcells" (sub-mesh of the Q2 patch, or None), "bc" -> {{axis, value, component, prescribed}, ...},
     "energyScale", "out" (CSV with the results), "nodal" (CSV with the nodal displacements).
   Jobs whose "out" file exists already are skipped (so a crashed kernel can be restarted).
   Plane stress in 2D ("SolidMechanicsModelForm" -> "PlaneStress", thickness 1), 3D elasticity in 3D; E = 1, nu = 3/10.

   Two values of K are written:
     K_system    = energyScale * U.S.U with S the assembled stiffness matrix (before the Dirichlet rows are applied)
                   and U the solution vector of the same discretization (NDSolve`ProcessEquations, DiscretizePDE,
                   DeployBoundaryConditions, LinearSolve) -- this is a(u_h, u_h), the finite element value;
     K_nintegrate = energyScale * NIntegrate of the strain energy density of the InterpolatingFunction returned by
                   NDSolveValue (the usual post-processing; it differs from a(u_h, u_h), see results/VERDICTS.md).
   Q2_nintegrate = NIntegrate of sigma_xx of the InterpolatingFunction over the patch sub-mesh / its area.
   The nodal displacements are written so that run_ndsolve.py can evaluate Q2 of u_h exactly.
   Values are FLOAT (no error control); checked by the agent only. *)

Needs["NDSolve`FEM`"];
nu = 3/10; ee = 1;
str[x_] := ToString[CForm[N[x]]];

runJob[job_] := Module[{dim, crd, inc, elem, mesh, vars, X, ops, bcs, t, ts, sol, dens, kn, tk, q2 = "nan",
    area = "nan", pmesh, sxx, vol, lam, mu, state, femdata, sd, dpde, dbcs, load, stiff, damp, mass, S0, U, ks, n},
  dim = job["dim"];
  crd = N[Import[job["coords"], "CSV"]];
  inc = Round[Import[job["cells"], "CSV"]];
  n = Length[crd];
  elem = If[dim == 2, TriangleElement[inc], TetrahedronElement[inc]];
  mesh = ToElementMesh["Coordinates" -> crd, "MeshElements" -> {elem}];
  vol = Total[Flatten[mesh["MeshElementMeasure"]]];
  X = If[dim == 2, {x, y}, {x, y, z}];
  vars = If[dim == 2, {u[x, y], v[x, y]}, {u[x, y, z], v[x, y, z], w[x, y, z]}];
  ops = SolidMechanicsPDEComponent[{vars, X},
    If[dim == 2, <|"YoungModulus" -> ee, "PoissonRatio" -> nu, "SolidMechanicsModelForm" -> "PlaneStress",
      "Thickness" -> 1|>, <|"YoungModulus" -> ee, "PoissonRatio" -> nu|>]];
  bcs = Map[DirichletCondition[vars[[#[[3]] + 1]] == #[[4]], X[[#[[1]] + 1]] == #[[2]]] &, Rationalize[job["bc"], 0]];
  (* the finite element system; NDSolve`ProcessEquations reads the vector equation unevaluated, so the zero vector
     is written out literally (ConstantArray[0, dim] there, or a pre-evaluated equation, gives "underdetermined") *)
  {ts, U} = AbsoluteTiming[
    {state} = If[dim == 2,
      NDSolve`ProcessEquations[{ops == {0, 0}, bcs}, Head /@ vars, X \[Element] mesh],
      NDSolve`ProcessEquations[{ops == {0, 0, 0}, bcs}, Head /@ vars, X \[Element] mesh]];
    femdata = state["FiniteElementData"];
    sd = state["SolutionData"][[1]];
    dpde = DiscretizePDE[femdata["PDECoefficientData"], femdata["FEMMethodData"], sd];
    dbcs = DiscretizeBoundaryConditions[femdata["BoundaryConditionData"], femdata["FEMMethodData"], sd];
    {load, stiff, damp, mass} = dpde["SystemMatrices"];
    S0 = stiff;
    DeployBoundaryConditions[{load, stiff}, dbcs];
    Flatten[LinearSolve[stiff, load]]];
  ks = job["energyScale"] U . S0 . U;
  Export[job["nodal"], Transpose[Partition[U, n]], "CSV"];
  (* the usual route: NDSolveValue and NIntegrate of the InterpolatingFunction *)
  {t, sol} = AbsoluteTiming[NDSolveValue[{ops == ConstantArray[0, dim], bcs}, Head /@ vars, X \[Element] mesh]];
  If[dim == 2,
    Module[{exx, eyy, gxy},
      exx = D[sol[[1]] @@ X, x]; eyy = D[sol[[2]] @@ X, y];
      gxy = D[sol[[1]] @@ X, y] + D[sol[[2]] @@ X, x];
      dens = ee/(1 - nu^2) (exx^2 + eyy^2 + 2 nu exx eyy) + ee/(2 (1 + nu)) gxy^2;
      sxx = ee/(1 - nu^2) (exx + nu eyy)],
    Module[{exx, eyy, ezz, gxy, gyz, gxz},
      lam = nu/((1 + nu) (1 - 2 nu)); mu = 1/(2 (1 + nu));
      exx = D[sol[[1]] @@ X, x]; eyy = D[sol[[2]] @@ X, y]; ezz = D[sol[[3]] @@ X, z];
      gxy = D[sol[[1]] @@ X, y] + D[sol[[2]] @@ X, x];
      gyz = D[sol[[2]] @@ X, z] + D[sol[[3]] @@ X, y];
      gxz = D[sol[[1]] @@ X, z] + D[sol[[3]] @@ X, x];
      dens = lam (exx + eyy + ezz)^2 + 2 mu (exx^2 + eyy^2 + ezz^2) + mu (gxy^2 + gyz^2 + gxz^2)]];
  {tk, kn} = AbsoluteTiming[job["energyScale"] NIntegrate[dens, X \[Element] mesh]];
  If[StringQ[job["patchcells"]],
    pmesh = ToElementMesh["Coordinates" -> N[Import[job["patchcoords"], "CSV"]],
      "MeshElements" -> {TriangleElement[Round[Import[job["patchcells"], "CSV"]]]}];
    area = NIntegrate[1, X \[Element] pmesh];
    q2 = NIntegrate[sxx, X \[Element] pmesh]/area];
  Export[job["out"], {{"ndof", "nelem", "time_s", "ndsolvevalue_time_s", "integrate_time_s", "K_system",
      "K_nintegrate", "Q2_nintegrate", "patch_area_nintegrate", "mesh_measure"},
     {dim n, Length[inc], str[ts], str[t], str[tk], str[ks], str[kn], str[q2], str[area], str[vol]}}, "CSV"];
  Print[job["name"], "  ndof ", dim n, "  t ", str[ts], " s  K_system = ", str[ks], "  K_nintegrate = ", str[kn],
    "  Q2_nintegrate = ", str[q2], "  mesh measure ", str[vol]];
];

jobs = Get[$ScriptCommandLine[[2]]];
Print["Mathematica ", $Version, "; ", Length[jobs], " jobs"];
Do[If[! FileExistsQ[job["out"]], runJob[job]], {job, jobs}];
