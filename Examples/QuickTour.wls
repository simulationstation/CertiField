#!/usr/bin/env wolframscript
(* A two-minute tour of CertiField.  Run with:  wolframscript -file Examples/QuickTour.wls
   (from a checkout; with the paclet installed, the PacletDirectoryLoad line is unnecessary) *)

PacletDirectoryLoad[DirectoryName[$InputFileName, 2]];
Needs["CertiField`"];

Print["CertiField ", $CertiFieldVersion, "   models: ", $CertiFieldModels, "\n"];

(* 1. A textbook case with a known answer: the phi^4 kink has energy 2 Sqrt[2]/3 *)
kink = CertifiedVariationalBound["Phi4Kink"];
Print[kink["Summary"]];
Print["   exact value 2 Sqrt[2]/3 = ", N[2 Sqrt[2]/3, 16], "\n"];

(* 2. Your own candidate.  A rough guess gives an honest, wide gap; ask for a tolerance and it is Inconclusive *)
rough = CertifiedVariationalBound["Phi4Kink", Function[x, Tanh[x]], "GapTolerance" -> 10^-3];
Print[rough["Summary"], "\n"];

(* 3. A boundary-value problem solved by NDSolve, certified on a finite domain *)
sol = NDSolveValue[{u''[x] == u[x]^3 - u[x], u[0] == -1, u[3] == 1/2}, u, {x, 0, 3}];
bvp = CertifiedVariationalBound[CertiFieldModel["ScalarField1D",
   <|"Potential" -> Function[s, (s^2 - 1)^2/4], "Domain" -> {0, 3}, "BoundaryValues" -> {-1, 1/2}|>], sol];
Print[bvp["Summary"]];
Print["   theorem: ", TextString[bvp["Theorem"]], "\n"];

(* 4. The Skyrme model: a lower bound over every degree-one map (about 30 s) *)
sk = CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> 1];
Print[sk["Summary"]];
Print["   ", TextString[sk["Theorem"]]];
Print["   not claimed: ", sk["NotClaimed"], "\n"];

(* 5. Your own energy in two or three dimensions: E[u] = Integrate[|grad u|^2/2 + W(u)] on a box, u = 0 on the boundary.
      W is convex here, so the lower bound is the complementary-energy bound and takes seconds, also in 3D.
      (A W that is not convex is handled too; the field is then optimised, which takes from under a minute.) *)
box = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND",
   <|"Potential" -> Function[u, u^4/4 + u^2/2 - u], "Domain" -> {{0, 1}, {0, 1}, {0, 1}}, "Cells" -> 32, "BoundaryData" -> 0|>]];
Print[box["Summary"]];
Print["   ", TextString[box["Theorem"]], "\n"];

(* 6. Anyone can re-check a certificate *)
Print["VerifyCertificate[kink] = ", VerifyCertificate[kink]];
Print["a tampered certificate:   ", VerifyCertificate[CertifiedBound[Append[kink["Association"], "LowerBound" -> 1.]]]];
