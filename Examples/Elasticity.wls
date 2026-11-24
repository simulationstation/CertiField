#!/usr/bin/env wolframscript
(* Linear elasticity with CertiField (a few seconds).  Run with:  wolframscript -file Examples/Elasticity.wls
   (from a checkout; with the paclet installed, the PacletDirectoryLoad line is unnecessary) *)

PacletDirectoryLoad[DirectoryName[$InputFileName, 2]];
Needs["CertiField`"];

(* 1. A block in plane strain: clamped on the left face, pulled by 0.1 on the right face, no traction on the other
      two.  E[u] = Integrate[e(u).D.e(u)/2 - b.u], isotropic with Lame constants 3/2 and 1 (Poisson ratio 0.3).
      Without a candidate the package computes the finite-element solution on the grid. *)
block = CertiFieldModel["LinearElasticity", <|"Lame" -> {3/2, 1}, "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 32, "Degree" -> 2,
    "TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}}, "BoundaryData" -> Function[{x, y}, {x/10, 0}]|>];
cert = CertifiedVariationalBound[block];
Print[cert["Summary"]];
Print["   ", TextString[cert["Theorem"]]];
Print["   not claimed: ", cert["NotClaimed"], "\n"];

(* 2. Your own candidate: a finite-element solution from NDSolve.  The lower bound does not depend on it. *)
sol = With[{lam = 3/2, mu = 1},
   NDSolveValue[{
     Inactive[Div][{{-(lam + 2 mu), 0}, {0, -mu}} . Inactive[Grad][u1[x, y], {x, y}], {x, y}] +
       Inactive[Div][{{0, -lam}, {-mu, 0}} . Inactive[Grad][u2[x, y], {x, y}], {x, y}] == 0,
     Inactive[Div][{{0, -mu}, {-lam, 0}} . Inactive[Grad][u1[x, y], {x, y}], {x, y}] +
       Inactive[Div][{{-mu, 0}, {0, -(lam + 2 mu)}} . Inactive[Grad][u2[x, y], {x, y}], {x, y}] == 0,
     DirichletCondition[{u1[x, y] == 0, u2[x, y] == 0}, x == 0], DirichletCondition[{u1[x, y] == 1/10, u2[x, y] == 0}, x == 1]},
    {u1, u2}, {x, 0, 1}, {y, 0, 1}]];
fem = CertifiedVariationalBound[block, Function[{x, y}, {sol[[1]][x, y], sol[[2]][x, y]}]];
Print[fem["Summary"], "\n"];

(* 3. A rough candidate gets an honest, wide gap *)
rough = CertifiedVariationalBound[block, Function[{x, y}, {x/10, 0}]];
Print[rough["Summary"], "\n"];

(* 4. A problem with a known answer: a bar in pure bending (stress (y, 0, 0), free sides).  The minimum is 1/32. *)
bar = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", <|"Lame" -> {1, 1}, "Domain" -> {{0, 2}, {-1/2, 1/2}}, "Cells" -> {4, 2}, "Degree" -> 2,
     "TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}}, "BoundaryData" -> Function[{x, y}, {3 x y/8, -y^2/16 - 3 x^2/16}]|>]];
Print[bar["Summary"]];
Print["   1/32 inside: ", bar["LowerBound"] <= 1/32 <= bar["UpperBound"], "\n"];

(* 5. A column under its own weight: clamped at the bottom, three free faces, body force (0, -1) *)
column = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", <|"Lame" -> {3/2, 1}, "Domain" -> {{0, 1}, {0, 2}}, "Cells" -> {16, 32}, "Degree" -> 2,
     "BodyForce" -> {0, -1}, "TractionFreeFaces" -> {{1, "Lower"}, {1, "Upper"}, {2, "Upper"}}, "BoundaryData" -> {0, 0}|>]];
Print[column["Summary"], "\n"];

(* 6. Anyone can re-check a certificate *)
Print["VerifyCertificate[cert] = ", VerifyCertificate[cert]];
Print["a tampered certificate:   ", VerifyCertificate[CertifiedBound[Append[cert["Association"], "LowerBound" -> 0.01503]]]];
