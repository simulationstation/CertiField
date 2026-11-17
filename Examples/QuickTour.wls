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
Print["   theorem: ", ToString[bvp["Theorem"]], "\n"];

(* 4. The Skyrme model: a lower bound over every degree-one map (about 30 s) *)
sk = CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> 1];
Print[sk["Summary"]];
Print["   ", ToString[sk["Theorem"]]];
Print["   not claimed: ", sk["NotClaimed"], "\n"];

(* 5. Anyone can re-check a certificate *)
Print["VerifyCertificate[kink] = ", VerifyCertificate[kink]];
Print["a tampered certificate:   ", VerifyCertificate[CertifiedBound[Append[kink["Association"], "LowerBound" -> 1.]]]];
