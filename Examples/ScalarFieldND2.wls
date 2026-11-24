#!/usr/bin/env wolframscript
(* Three certified examples of ScalarFieldND with a source term, a prescribed mean, periodic directions and
   potentials on a value range (docs/SCALAR_ND_2.md).  Run from the repository root:
     wolframscript -file Examples/ScalarFieldND2.wls                                  (a few minutes) *)

PacletDirectoryLoad[DirectoryName[$InputFileName, 2]];
Needs["CertiField`"];

show[name_, c_] := Print[name, "\n    ", c["Summary"], "\n    ", Round[c["Seconds"], 0.1], " s, certificate re-verified: ", VerifyCertificate[c], "\n"];

(* 1. Phase separation of a binary mixture: regular-solution free energy with chi = 3 on [0, 1], mean composition 1/2,
      composition 1/10 held on the face x1 = 0, the face x1 = 4 free, periodic in x2. *)
show["Phase separation (regular-solution free energy, prescribed mean)",
  CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|
    "Potential" -> Function[u, u Log[u] + (1 - u) Log[1 - u] + 3 u (1 - u)], "ValueRange" -> {0, 1}, "Mean" -> 1/2,
    "Domain" -> {{0, 4}, {0, 1/2}}, "Cells" -> {32, 4}, "BoundaryData" -> 1/10, "NaturalFaces" -> {{1, "Upper"}}, "PeriodicDirections" -> {2}|>]]];

(* 2. Image restoration with quadratic smoothing and a robust data term tabulated per pixel: a 16 x 16 image of a
      bright square on a dark background with noise and a few outliers; labels 0, 1/8, ..., 1; cost
      4 min((t - f)^2, 1/16).  Total-variation smoothing is not covered, and uncertified global methods (functional
      lifting) already solve problems of this kind; what is added here is a certified interval for the minimum. *)
img = BlockRandom[SeedRandom[3];
  Table[If[5 <= i <= 12 && 5 <= j <= 12, 0.8, 0.2] + RandomReal[{-0.1, 0.1}] + If[RandomReal[] < 0.05, 0.6 RandomChoice[{-1, 1}], 0.], {i, 16}, {j, 16}]];
show["Image restoration (tabulated data term, quadratic smoothing)",
  CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|
    "Stiffness" -> 1/100, "ValueRange" -> {0, 1}, "DataTerm" -> Table[4 Min[(t - img[[i, j]])^2, 1/16], {i, 16}, {j, 16}, {t, 0, 1, 1/8}],
    "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> {16, 16}, "NaturalFaces" -> {{1, "Lower"}, {1, "Upper"}, {2, "Lower"}, {2, "Upper"}}|>],
    Method -> {"Calibration", "KnotIntervals" -> 8}]];

(* 3. A uniaxial ferroelectric film in an applied field: Landau energy -P^2/2 + P^4/4, gradient energy |grad P|^2/2,
      the field E(x) = (x1 - 1)/5 pointing down on the left half and up on the right half (energy density -E P), all
      faces free: the minimiser has a domain wall in the middle. *)
show["Ferroelectric domains in an applied field (source term)",
  CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|
    "Potential" -> Function[p, -p^2/2 + p^4/4], "Source" -> Function[{x, y}, (x - 1)/5],
    "Domain" -> {{0, 2}, {0, 1}}, "Cells" -> {32, 16}, "NaturalFaces" -> {{1, "Lower"}, {1, "Upper"}, {2, "Lower"}, {2, "Upper"}}|>]]];
