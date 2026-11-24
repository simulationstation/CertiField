(* ::Package:: *)

(* Skyrme fields of degree two: the level-set lower bound with the flux of degree B (docs/SKYRME_DEGREE_TWO.md).
   Five runs of the dynamic programme on Quick grids, about five minutes.
   Run with:  wolframscript -file tools/run_tests.wls --skyrme-degree *)

Needs["CertiField`"];

sdTab = CertiField`Private`minorantTable[{2000, 900, 1/10^10}];
sdClaims = CertiField`Private`skyrmeDegreeClaims[];
sdRat = CertiField`Private`skyrmeDegreeRational;

(* ---------- degree one: the new function performs the operations of the existing one ---------- *)

VerificationTest[
  Module[{g = CertiField`Private`loadSkyrmeGrid[CertiField`Private`skyrmeGridDir[1, "Quick"]], a, b},
    a = CertiField`Private`levelSetLowerBound[1, g, sdTab];
    b = CertiField`Private`levelSetLowerBoundDegree[1, 1, g, sdTab];
    {a["LowerBound"] === b["LowerBound"], a["Pairs"] === b["Pairs"], a["Fallbacks"] === b["Fallbacks"],
     a["FarBinding"] === b["FarBinding"], b["Pairs"]}],
  {True, True, True, True, 5025904}, TestID -> "degree-one-same-numbers-as-the-existing-function"]

(* ---------- degree two, massless ---------- *)

sd2 = CertifiedVariationalBound["SkyrmeDegreeTwo", "PionMass" -> 0];

VerificationTest[sd2["Status"], "Certified", TestID -> "degree-two-status"]

(* the floating-point value of the relaxed problem (the extremal of Reference/SkyrmeDegree/relaxed.py) is 2.2084398;
   a certified lower bound of the same problem must lie below it, and on the Quick grid it lies less than 0.001 below *)
VerificationTest[2.2075 < sd2["LowerBound"] < 2.2084398, True, TestID -> "degree-two-below-and-near-the-relaxed-value"]

(* above Faddeev's bound 2 by more than one per cent: the bound is not a negative result *)
VerificationTest[sd2["LowerBound"] > 2.02, True, TestID -> "degree-two-exceeds-faddeev-bound"]

(* constraints that must hold or there is an error in the proof or in the code:
   certified upper bound for I_2 of the degree-two note (2.358835); twice the value of the degree-one relaxed problem
   (1.19216, Proposition E of the note); twice the certified hedgehog upper bound *)
VerificationTest[
  {sd2["LowerBound"] <= 2.358835, sd2["LowerBound"] <= 2 1.19216, sd2["LowerBound"] <= sd2["UpperBound"]},
  {True, True, True}, TestID -> "degree-two-does-not-exceed-known-upper-bounds"]

(* upper bound: two disjoint copies of the hedgehog trial map (1.231447 each) *)
VerificationTest[2.46289 < sd2["UpperBound"] < 2.4629, True, TestID -> "degree-two-upper-bound-is-twice-the-hedgehog"]

VerificationTest[StringContainsQ[ToString[sd2["Theorem"]], "degree 2"], True, TestID -> "degree-two-theorem-names-the-degree"]

VerificationTest[VerifyCertificate[sd2], True, TestID -> "degree-two-verifies"]

VerificationTest[VerifyCertificate[CertifiedBound[Append[sd2["Association"], "LowerBound" -> 2.2085]]], False,
  TestID -> "degree-two-rejects-a-tampered-lower-bound"]

VerificationTest[VerifyCertificate[CertifiedBound[Append[sd2["Association"], "UpperBound" -> 2.46]]], False,
  TestID -> "degree-two-rejects-a-tampered-upper-bound"]

(* the stored claim is the computed bound lowered by one part in 10^9 *)
VerificationTest[
  CertiField`Private`exactOf[sd2["Certificate"]["LevelSet"]["LowerBound"]] >= sdRat[sdClaims["b2_mu0_quick", "lower"]] > 2,
  True, TestID -> "degree-two-stored-claim-mu0"]

(* ---------- degree two, pion mass one: path test and stored claim ---------- *)
(* At every level the certified function, at the volume of the floating-point extremal, must not exceed the cost of
   that extremal beyond the level: the extremal is an admissible path, so its cost bounds the value function from
   above.  A violation would mean an error in the programme (or in the floating-point path). *)

sdRun = Module[{dir = CertiField`Private`skyrmeDegreeGridDir[2, 1, "Quick"], grid, res, ctg},
  grid = CertiField`Private`loadSkyrmeDegreeGrid[2, 1, "Quick"];
  res = CertiField`Private`levelSetLowerBoundDegree[2, 1, grid, sdTab, "Trace" -> True];
  ctg = CertiField`Private`readF64[FileNameJoin[{dir, "ctg.f64"}]];
  <|"Result" -> res, "Float" -> grid["Meta"]["extremal"]["value_float"],
    "Slack" -> Flatten[Table[{ctg[[r[[1]] + 1]] - r[[4]], ctg[[r[[1]] + 1]] - r[[5]]}, {r, res["Trace"]}]]/(12 Pi^2)|>];

VerificationTest[Length[sdRun["Slack"]] > 1000 && Min[sdRun["Slack"]] >= 0, True, TestID -> "degree-two-path-test-mu1"]

VerificationTest[sdRun["Float"] - 0.002 < sdRun["Result"]["LowerBound"] < sdRun["Float"], True,
  TestID -> "degree-two-mu1-below-and-near-the-relaxed-value"]

VerificationTest[
  CertiField`Private`exactOf[sdRun["Result"]["LowerBound"]] >= sdRat[sdClaims["b2_mu1_quick", "lower"]] > 2,
  True, TestID -> "degree-two-stored-claim-mu1"]

(* certified upper bound for I_2 at mu = 1 of the degree-two note (2.721257); twice the degree-one relaxed value 1.41444 *)
VerificationTest[{sdRun["Result"]["LowerBound"] <= 2.721257, sdRun["Result"]["LowerBound"] <= 2 1.41444}, {True, True},
  TestID -> "degree-two-mu1-does-not-exceed-known-upper-bounds"]

(* ---------- interface ---------- *)

VerificationTest[MemberQ[$CertiFieldModels, "SkyrmeDegreeTwo"] && CertiFieldModelInfo["SkyrmeDegreeTwo"]["Status"] === "Supported",
  True, TestID -> "degree-two-registered"]

VerificationTest[FailureQ[CertifiedVariationalBound["SkyrmeDegreeTwo", "PionMass" -> 1/2]], True,
  TestID -> "degree-two-unsupported-mass-refused"]

VerificationTest[FailureQ[CertifiedVariationalBound["SkyrmeDegreeTwo", "Resolution" -> "Publication"]], True,
  TestID -> "degree-two-unsupported-resolution-refused"]

VerificationTest[
  {sdRat["2.2079"], sdRat["-0.5"], sdRat["3"], FailureQ[Catch[sdRat["2e3"], "CertiField"]]},
  {22079/10000, -1/2, 3, True}, TestID -> "degree-two-decimal-strings"]
