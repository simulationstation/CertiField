#!/usr/bin/env wolframscript
(* Certified lower bounds for Skyrme fields of degree B from the grids in Data/SkyrmeDegree.

     wolframscript -file tools/skyrme_degree_claims.wls
         recomputes every stored claim of Data/SkyrmeDegree/claims.json from its grid and compares
     wolframscript -file tools/skyrme_degree_claims.wls --write b2_mu0_quick b2_mu1_quick
         runs the certified programme on the named grids and stores their claims

   A stored claim is the computed lower bound lowered by one part in 10^9 and rounded down to ten digits.
   With --write the run also performs the path test: at every level the certified function, at the volume of the
   floating-point extremal of the relaxed problem, must not exceed the cost of that extremal beyond the level. *)

root = DirectoryName[$InputFileName, 2];
PacletDirectoryLoad[root];
Needs["CertiField`"];
p = "CertiField`Private`";
sym[name_String] := Symbol[p <> name];
dataDir = FileNameJoin[{root, "Data", "SkyrmeDegree"}];
claimsFile = FileNameJoin[{dataDir, "claims.json"}];
claims = If[FileExistsQ[claimsFile], Import[claimsFile, "RawJSON"], <||>];
args = Rest[$ScriptCommandLine];
bad = 0;

If[MemberQ[args, "--write"],
  Do[
    Module[{meta, deg, muq, preset, grid, tab, res, lb, ctg, slack, claim},
      meta = Import[FileNameJoin[{dataDir, name, "meta.json"}], "RawJSON"];
      deg = meta["degree"]; muq = sym["skyrmeDegreeRational"][meta["mu_exact"]];
      preset = Capitalize[meta["preset"]];
      grid = sym["loadSkyrmeGrid"][FileNameJoin[{dataDir, name}]];
      tab = sym["minorantTable"][{2000, 900, 1/10^10}];
      res = sym["levelSetLowerBoundDegree"][deg, muq, grid, tab, "Trace" -> True];
      lb = res["LowerBound"];
      ctg = sym["readF64"][FileNameJoin[{dataDir, name, "ctg.f64"}]];
      slack = Flatten[Table[{ctg[[r[[1]] + 1]] - r[[4]], ctg[[r[[1]] + 1]] - r[[5]]}, {r, res["Trace"]}]]/(12 Pi^2);
      claim = sym["lowerString"][SetPrecision[lb, Infinity] (1 - 1/10^9), 10];
      Print[name, ":  degree ", deg, ", mu = ", muq, ", ", preset, ":  certified E/(12 pi^2) >= ", ToString[lb, InputForm],
        "   claim ", claim, "   float value along the extremal ", ToString[NumberForm[meta["extremal"]["value_float"], 8]],
        "   levels ", res["Levels"], ", transitions ", res["Pairs"], ", fallbacks ", res["Fallbacks"], ", far ", res["FarBinding"],
        "   path test: smallest slack ", ToString[ScientificForm[Min[slack], 3]], " over ", Length[slack]/2, " levels   ",
        ToString[Round[res["Seconds"]]], " s"];
      If[! (Min[slack] >= 0 && lb <= meta["extremal"]["value_float"]),
        Print["  **** PATH TEST FAILED: the certified function exceeds the cost of an admissible path ****"]; bad++];
      claims[name] = <|"degree" -> deg, "mu" -> meta["mu_exact"], "preset" -> preset, "lower" -> claim,
        "computed" -> ToString[lb, InputForm], "levels" -> res["Levels"], "transitions" -> res["Pairs"],
        "fallbacks" -> res["Fallbacks"], "far_binding" -> res["FarBinding"],
        "relaxed_value_float" -> meta["extremal"]["value_float"], "path_test_smallest_slack" -> Min[slack]|>],
    {name, DeleteCases[args, "--write"]}];
  If[bad == 0,
    Export[claimsFile, KeySort[claims], "RawJSON"];
    Print["wrote ", claimsFile]],
  (* verification of the stored claims *)
  Do[
    Module[{ok = sym["verifySkyrmeDegreeClaim"][name]},
      Print[name, ":  claim E/(12 pi^2) >= ", claims[name, "lower"], "   ", If[TrueQ[ok], "OK", "FAILED"]];
      If[! TrueQ[ok], bad++]],
    {name, Keys[claims]}];
  Print[If[bad == 0, "ALL CLAIMS VERIFIED", ToString[bad] <> " PROBLEM(S)"]]];
Exit[If[bad == 0, 0, 1]];
