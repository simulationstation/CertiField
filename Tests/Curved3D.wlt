(* ::Package:: *)

(* ElasticCavity3D suite (about three minutes).  Run with:  wolframscript -file tools/run_tests.wls --curved3d
   The stored grid-8 certificates of the Python reference (Reference/Curved3D/certificates) are re-evaluated from their
   fields; the two implementations round differently, by about 10^-12, so their bounds must agree to 10^-6. *)

Needs["CertiField`"];

c3Dir = FileNameJoin[{CertiField`Private`$CertiFieldDirectory, "Reference", "Curved3D", "certificates"}];
(* a random coefficient array that satisfies the ties of the model (grid n, frictionless grips) *)
c3RandomField[n_, kind_] := Module[{t = CertiField`Private`c3Ties[n, kind, False]},
  SeedRandom[11];
  Developer`ToPackedArray[ArrayReshape[Normal[t[[1]] . RandomReal[{-1, 1}, Dimensions[t[[1]]][[2]]]], {If[kind === "Psi", 6, 3], n + 2, n + 2, n + 2}]]];
c3SolidModel = Sequence["Grips" -> "Frictionless", "Cells" -> 4, "Cavity" -> False, "Depths" -> {0, 0, 0, 0}];

(* ---------- known answer: the solid block between frictionless grips has K = 2 (uniform tension) ---------- *)

c3Solid = CertifiedVariationalBound["ElasticCavity3D", c3SolidModel];

VerificationTest[c3Solid["Status"], "Certified", TestID -> "c3-solid-block-status"]
VerificationTest[c3Solid["LowerBound"] <= 2 <= c3Solid["UpperBound"] && c3Solid["UpperBound"] - c3Solid["LowerBound"] < 10^-8, True,
  TestID -> "c3-solid-block-encloses-2"]
VerificationTest[VerifyCertificate[c3Solid], True, TestID -> "c3-verify-genuine-certificate"]

(* ---------- tampering: a raised bound, a changed field that keeps the ties, a field that breaks them ---------- *)

VerificationTest[
  With[{a = c3Solid["Association"], c = c3Solid["Certificate"]},
    {VerifyCertificate[CertifiedBound[Append[a, "LowerBound" -> 2.000001]]],
     VerifyCertificate[CertifiedBound[Append[a, "Certificate" -> Append[c, "Psi" -> c["Psi"] + 0.01 c3RandomField[4, "Psi"]]]]],
     VerifyCertificate[CertifiedBound[Append[a, "Certificate" -> Append[c, "Psi" -> c["Psi"] + 0.01]]]],
     VerifyCertificate[CertifiedBound[Append[a, "Certificate" -> Append[c, "V" -> ReplacePart[c["V"], {1, 6, 2, 2} -> 1.]]]]]}],
  {False, False, False, False}, TestID -> "c3-verify-rejects-tampered-certificates"]

(* ---------- machine overflow ends in "Inconclusive", never in a number ---------- *)

VerificationTest[
  {CertifiedVariationalBound["ElasticCavity3D", c3SolidModel, "Field" -> 1.*^200 c3RandomField[4, "Psi"]],
   CertifiedVariationalBound["ElasticCavity3D", 1.*^200 c3RandomField[4, "V"], c3SolidModel]},
  {Failure["Inconclusive", _], Failure["Inconclusive", _]}, SameTest -> MatchQ, TestID -> "c3-overflow-is-inconclusive"]

(* ---------- refused inputs ---------- *)

VerificationTest[
  {FailureQ[CertifiedVariationalBound["ElasticCavity3D", "Cells" -> 6]],
   FailureQ[CertifiedVariationalBound["ElasticCavity3D", "Grips" -> "Glued"]],
   CertifiedVariationalBound["ElasticCavity3D", c3SolidModel, "Field" -> c3RandomField[4, "Psi"] + 0.5],
   CertifiedVariationalBound["ElasticCavity3D", RandomReal[1, {3, 5, 5, 5}], c3SolidModel]},
  {True, True, Failure["BadField", _], Failure["BadCandidate", _]}, SameTest -> MatchQ, TestID -> "c3-bad-inputs-refused"]

(* ---------- the cavity on a coarse grid, solver included: consistent with the certified grid-24 range ---------- *)

VerificationTest[
  With[{c = CertifiedVariationalBound["ElasticCavity3D", "Grips" -> "Clamped", "Cells" -> 4, "Depths" -> {1, 2, 2, 3}]},
    {c["Status"], c["LowerBound"] <= 1.833487, c["UpperBound"] >= 1.817322, c["LowerBound"] > 1.5,
     StringContainsQ[ToString[c["Theorem"]], "closed ball"] && StringContainsQ[c["Summary"], "<= K <="]}],
  {"Certified", True, True, True, True}, TestID -> "c3-coarse-cavity-consistent"]

(* ---------- grid 8: the stored fields of the reference give the stored bounds ---------- *)

VerificationTest[
  Table[
    With[{f = FileNameJoin[{c3Dir, "cert_" <> ToLowerCase[g] <> "_N8_L2335.npz"}]},
      With[{z = CertiField`Private`c3ReadNPZ[f],
            c = CertifiedVariationalBound["ElasticCavity3D", f, "Grips" -> g, "Cells" -> 8, "Depths" -> {2, 3, 3, 5}, "Field" -> f]},
        {c["Status"], Abs[c["LowerBound"] - z["K_low"]] < 10^-6, Abs[c["UpperBound"] - z["K_up"]] < 10^-6}]],
    {g, {"Clamped", "Frictionless"}}],
  {{"Certified", True, True}, {"Certified", True, True}}, TestID -> "c3-grid-8-reproduces-stored-bounds"]
