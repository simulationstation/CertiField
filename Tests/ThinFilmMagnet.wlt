(* ::Package:: *)

(* ThinFilmMagnet suite.  Run with:  wolframscript -file tools/run_tests.wls --magnet *)

Needs["CertiField`"];

tfmRoot = CertiField`Private`$CertiFieldDirectory;
tfmTV = FileNameJoin[{tfmRoot, "Reference", "Magnetisation", "testvectors", "small"}];
tfmExp = Import[FileNameJoin[{tfmTV, "expected.json"}], "RawJSON"];
tfmRead[f_] := Developer`ToPackedArray[BinaryReadList[FileNameJoin[{tfmTV, f}], "Real64", ByteOrdering -> -1]];

(* the certified energy and degree of a perturbed field agree with the Python reference on identical inputs *)
VerificationTest[
  Module[{x = tfmRead["x1.f64"], y = tfmRead["x2.f64"], Mn = ArrayReshape[tfmRead["Mn.f64"], tfmExp["dims"]], ok = True, m, mod, ub, deg},
    Do[m = tfmExp[nm, "model"];
      mod = <|"A" -> ToExpression[m["A"]], "DMI" -> ToExpression[m["kappa"]], "K" -> ToExpression[m["K"]], "Field" -> ToExpression[m["h"]]|>;
      ub = CertiField`Private`tfmUpperBound[mod, x, y, Mn];
      ok = ok && Abs[ub - tfmExp[nm, "upper"]] <= 10^-9 Abs[tfmExp[nm, "upper"]], {nm, {"zeeman", "all", "positiveK"}}];
    deg = CertiField`Private`tfmDegree[x, y, Mn];
    ok && deg[[1]] == tfmExp["degree"] && Abs[N[deg[[2, 1]]] - tfmExp["degree_enclosure"][[1]]] < 10^-9 && Abs[N[deg[[2, 2]]] - tfmExp["degree_enclosure"][[2]]] < 10^-9],
  True, TestID -> "tfm-competitor-matches-reference"]

(* the shipped dual certificate passes the exact check and gives the reference's exact bound; a lowered u is rejected *)
VerificationTest[
  Module[{file = FileNameJoin[{tfmRoot, "Data", "Magnetisation", "zeeman_h2.json"}], dual, bad},
    dual = CertiField`Private`tfmImportDual[file];
    bad = dual; bad["U"] = bad["U"] (1 - 1/10^6);
    {CertiField`Private`tfmCheckDual[dual], CertiField`Private`tfmBoundCoefficients[dual, -1, 1][[1]] == ToExpression[Import[file, "RawJSON"]["bound_coefficient_of_pi"]],
     CertiField`Private`tfmCheckDual[bad]}],
  {True, True, False}, TestID -> "tfm-dual-certificate-exact-check"]

(* known answer: exchange only, degree -1: inf E = 8 pi A, not attained; the lower bound is exactly 4 pi *)
tfmEx = CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet", <|"A" -> 1/2, "Degree" -> -1, "Domain" -> {{-1, 1}, {-1, 1}}, "Cells" -> 120|>],
  Function[r, 2 ArcTan[r, 1/10 (1 - r^2)^2]]];
VerificationTest[
  {tfmEx["Status"], tfmEx["LowerBound"] <= 4 Pi <= tfmEx["LowerBound"] + 10^-14, 4 Pi < tfmEx["UpperBound"] < 1.05 4 Pi},
  {"Certified", True, True}, TestID -> "tfm-exchange-only-belavin-polyakov"]

(* known answer: critical coupling (A = 1/2, kappa = 1/2, K = -1/2, h = 1): inf E = 2 pi on the plane *)
tfmCrit = CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet",
    <|"A" -> 1/2, "DMI" -> 1/2, "K" -> -1/2, "Field" -> 1, "Domain" -> {{-32, 32}, {-32, 32}}, "Cells" -> 200|>],
  Function[r, 2 ArcTan[r, (1 - (r/32)^2)^2]]];
VerificationTest[
  {tfmCrit["Status"], tfmCrit["LowerBound"] <= 2 Pi <= tfmCrit["LowerBound"] + 10^-14, 2 Pi < tfmCrit["UpperBound"] < 1.03 2 Pi},
  {"Certified", True, True}, TestID -> "tfm-critical-coupling"]

(* M1: the chiral skyrmion class of the Li-Melcher model at h = 2 on [-4, 4]^2 *)
tfmM1 = CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet",
  <|"A" -> 1/2, "DMI" -> 1, "Field" -> 2, "Degree" -> -1, "Domain" -> {{-4, 4}, {-4, 4}}, "Cells" -> 240|>]];
VerificationTest[
  {tfmM1["Status"], Abs[tfmM1["LowerBound"] - 8.44861324354419] < 10^-10, tfmM1["LowerBound"] < tfmM1["UpperBound"], tfmM1["RelativeGap"] < 0.005},
  {"Certified", True, True, True}, TestID -> "tfm-skyrmion-interval"]

VerificationTest[VerifyCertificate[tfmM1], True, TestID -> "tfm-verify-genuine-certificate"]

VerificationTest[
  {VerifyCertificate[CertifiedBound[Append[tfmM1["Association"], "LowerBound" -> tfmM1["LowerBound"] + 0.01]]],
   VerifyCertificate[CertifiedBound[Append[tfmM1["Association"], "UpperBound" -> tfmM1["UpperBound"] - 0.01]]],
   VerifyCertificate[CertifiedBound[ReplacePart[tfmM1["Association"], {"Certificate", "Dual", "U"} -> tfmM1["Association"]["Certificate"]["Dual"]["U"] (1 - 1/10^6)]]]},
  {False, False, False}, TestID -> "tfm-verify-rejects-tampering"]

(* a field value that is not in the cache: the package finds its own dual field (Python reference: 10.830954365) *)
tfmH4 = CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet", <|"DMI" -> 1, "Field" -> 4, "Domain" -> {{-3, 3}, {-3, 3}}, "Cells" -> 240|>]];
VerificationTest[
  {tfmH4["Status"], tfmH4["Certificate"]["Dual"]["Source"], tfmH4["LowerBound"] >= 10.830954365 (1 - 10^-3), tfmH4["RelativeGap"] < 0.005},
  {"Certified", "search", True, True}, TestID -> "tfm-dual-search-field-not-in-cache"]

(* kappa < 0: the cached field for |kappa| with phi -> -phi, checked exactly for the model itself; degree +1 *)
VerificationTest[
  Module[{neg = CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet", <|"DMI" -> -1, "Field" -> 2, "Domain" -> {{-4, 4}, {-4, 4}}, "Cells" -> 240|>]],
      plus = CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet", <|"DMI" -> 1, "Field" -> 2, "Degree" -> 1, "Domain" -> {{-4, 4}, {-4, 4}}, "Cells" -> 240|>]]},
    {neg["LowerBound"] == tfmM1["LowerBound"], neg["RelativeGap"] < 0.005, plus["LowerBound"] <= 4 Pi <= plus["LowerBound"] + 10^-14, plus["RelativeGap"] < 0.02}],
  {True, True, True, True}, TestID -> "tfm-negative-kappa-and-degree-plus-one"]

(* the stray field of a film of thickness t: known answer first (degree 0: the uniform state is the minimiser) *)
VerificationTest[
  With[{c = CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet", <|"K" -> 1, "Field" -> 1, "Degree" -> 0, "Thickness" -> 1/10,
      "Domain" -> {{-2, 2}, {-2, 2}}, "Cells" -> 20|>], Function[{x, y}, {0, 0, 1}]]},
    {c["Status"], c["LowerBound"] == 0, 0 <= c["UpperBound"] < 10^-12}],
  {"Certified", True, True}, TestID -> "tfm-stray-field-uniform-state"]

tfmStray = CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet",
  <|"A" -> 1/2, "DMI" -> 1, "K" -> 1/2, "Field" -> 2, "Thickness" -> 1/20, "Domain" -> {{-4, 4}, {-4, 4}}, "Cells" -> 240|>]];
VerificationTest[
  {tfmStray["Status"], 0 < tfmStray["LowerBound"] < tfmStray["UpperBound"], tfmStray["RelativeGap"] < 0.03, VerifyCertificate[tfmStray],
   VerifyCertificate[CertifiedBound[Append[tfmStray["Association"], "LowerBound" -> tfmStray["LowerBound"] + 0.001]]]},
  {"Certified", True, True, True, False}, TestID -> "tfm-stray-field-skyrmion"]

(* failure modes *)
VerificationTest[
  {FailureQ[CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet", <|"DMI" -> 1, "Field" -> 2, "Domain" -> {{-1/3, 1/3}, {-1, 1}}|>]]],
   FailureQ[CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet", <|"DMI" -> 1, "Field" -> 2, "Degree" -> -2, "Domain" -> {{-4, 4}, {-4, 4}}, "Cells" -> 40|>]]],
   FailureQ[CertifiedVariationalBound[CertiFieldModel["ThinFilmMagnet", <|"DMI" -> 1, "Field" -> 2, "Domain" -> {{-4, 4}, {-4, 4}}, "Cells" -> 40|>],
     Function[r, 0.1 r]]]},
  {True, True, True}, TestID -> "tfm-failure-modes"]
