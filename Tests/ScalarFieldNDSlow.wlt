(* ::Package:: *)

(* ScalarFieldND, long runs (about ten minutes).  Run with:  wolframscript -file tools/run_tests.wls --nd-slow
   1. Every certificate stored by the Python reference in Reference/ScalarND/results is recomputed by this
      implementation (the large ones only if their files are present).
   2. The acceptance runs of calibration mode (Reference/ScalarND/PORTING.md, Section 6): the wall (96 x 4 cells,
      32 knot intervals) and the two-dimensional ramp (64 x 64, 32), with the fields found by the Wolfram Language
      optimiser, against the reference's stored bounds. *)

Needs["CertiField`"];

ndRef = FileNameJoin[{CertiField`Private`$CertiFieldDirectory, "Reference", "ScalarND", "results"}];
ndClaims[name_String] := With[{m = Import[FileNameJoin[{ndRef, name <> ".json"}], "RawJSON"]},
  CertiField`Private`ndRational /@ {m["lower"], m["upper"]}];

(* ---------- the reference's certificates ---------- *)

ndStored = Select[FileBaseName /@ FileNames["*.json", ndRef], FileExistsQ[FileNameJoin[{ndRef, # <> ".npz"}]] &];

VerificationTest[Length[ndStored] >= 7, True, TestID -> "nd-slow-stored-certificates-present"]

VerificationTest[
  Select[ndStored, ! TrueQ[Quiet[Catch[CertiField`Private`ndVerifyReferenceCertificate[FileNameJoin[{ndRef, #}]], "CertiField"]["Verified"]]] &],
  {}, TestID -> "nd-slow-all-stored-reference-certificates-verified"]

(* ---------- acceptance: the wall ---------- *)
(* exact minimum 0.9449674266 (cross-section times the one-dimensional value); reference: 0.9442106791 <= inf E <= 0.9451112063 *)

ndWall2d = CertifiedVariationalBound[
  CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, (u^2 - 1)^2/4], "Domain" -> {{-3, 3}, {0, 1}}, "Cells" -> {96, 4},
    "BoundaryData" -> Function[{x, y}, If[x < 0, -1, 1]], "NaturalFaces" -> {{2, "Lower"}, {2, "Upper"}}|>],
  Function[{x, y}, Tanh[x/Sqrt[2]]], "RefineCandidate" -> True,
  Method -> {"Calibration", "TRange" -> {-1, 1}, "Levels" -> {{{24, 1}, 8}, {{48, 2}, 8}, {{48, 2}, 16}, {{48, 2}, 32}, {{96, 4}, 32}}}];

VerificationTest[ndWall2d["LowerBound"] <= 0.9449674266 <= ndWall2d["UpperBound"], True, TestID -> "nd-slow-wall-encloses-exact-value"]
VerificationTest[
  With[{ref = ndClaims["wall2d"]},
    {ndWall2d["LowerBound"] >= (1 - 1/1000) ref[[1]], ndWall2d["LowerBound"] <= ref[[2]], ndWall2d["UpperBound"] <= (1 + 1/1000) ref[[2]], ndWall2d["UpperBound"] >= ref[[1]]}],
  {True, True, True, True}, TestID -> "nd-slow-wall-within-0.1-percent-of-reference"]
VerificationTest[ndWall2d["RelativeGap"] < 0.0011, True, TestID -> "nd-slow-wall-gap"]
VerificationTest[VerifyCertificate[ndWall2d], True, TestID -> "nd-slow-wall-certificate-verifies"]

(* ---------- acceptance: the ramp ---------- *)
(* no closed form; reference: 8.275401140 <= inf E <= 8.291701687.  The two certificates are for the same problem,
   so each lower bound must lie below the other's upper bound. *)

ndRamp2d = CertifiedVariationalBound[
  CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, 4 (u^2 - 1)^2], "Domain" -> {{-1, 1}, {-1, 1}}, "Cells" -> {64, 64},
    "BoundaryData" -> Function[{x, y}, x]|>],
  Method -> {"Calibration", "KnotIntervals" -> 32}];

VerificationTest[
  With[{ref = ndClaims["ramp2d"]},
    {ndRamp2d["LowerBound"] >= (1 - 1/1000) ref[[1]], ndRamp2d["LowerBound"] <= ref[[2]], ndRamp2d["UpperBound"] <= (1 + 1/1000) ref[[2]], ndRamp2d["UpperBound"] >= ref[[1]]}],
  {True, True, True, True}, TestID -> "nd-slow-ramp-within-0.1-percent-of-reference"]
VerificationTest[ndRamp2d["Certificate"]["Field"]["TRange"], {-1, 1}, TestID -> "nd-slow-ramp-automatic-range"]
VerificationTest[ndRamp2d["RelativeGap"] < 0.0025, True, TestID -> "nd-slow-ramp-gap"]
VerificationTest[VerifyCertificate[ndRamp2d], True, TestID -> "nd-slow-ramp-certificate-verifies"]
