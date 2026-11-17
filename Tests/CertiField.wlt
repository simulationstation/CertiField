(* ::Package:: *)

(* Fast test suite (about half a minute).  Run with:  wolframscript -file tools/run_tests.wls *)

Needs["CertiField`"];

VerificationTest[StringQ[$CertiFieldVersion], True, TestID -> "version"]

(* ---------- rigorous primitives ---------- *)

VerificationTest[
  Module[{xs, ys, iv, ex},
    SeedRandom[7]; xs = RandomReal[{-10, 10}, 5000]; ys = RandomReal[{-10, 10}, 5000];
    ex = SetPrecision[#, Infinity] &;
    iv = CertiField`Private`ivMul[CertiField`Private`ivPt[xs], CertiField`Private`ivPt[ys]];
    And @@ MapThread[ex[#1] <= ex[#3] ex[#4] <= ex[#2] &, {iv[[1]], iv[[2]], xs, ys}]],
  True, TestID -> "outward-rounded-multiplication-encloses-exact-product"]

VerificationTest[
  Module[{v = Join[{0., 1.*^-300, 1.*^-201, 2.}, RandomReal[{0, 1000}, 2000]], ex = SetPrecision[#, Infinity] &},
    And @@ MapThread[ex[#1]^2 <= ex[#3] <= ex[#2]^2 &, {CertiField`Private`sqrtLo[v], CertiField`Private`sqrtHi[v], v}]],
  True, TestID -> "verified-square-roots"]

VerificationTest[
  With[{b = CertiField`Private`sinBr[7/5]}, b[[1]] <= Sin[7/5] <= b[[2]] && b[[2]] - b[[1]] < 10^-39],
  True, TestID -> "sine-bracket"]

VerificationTest[
  With[{q = 1/3}, SetPrecision[CertiField`Private`fLo[q], Infinity] <= q <= SetPrecision[CertiField`Private`fHi[q], Infinity]],
  True, TestID -> "directed-rounding"]

(* a minorant entry above h_kappa must be rejected; one below must be accepted *)
VerificationTest[
  {CertiField`Private`certifyMinorantMask[{2.1}, {0.}, 1.], CertiField`Private`certifyMinorantMask[{1.9}, {0.}, 1.]},
  {{0}, {1}}, TestID -> "minorant-certification-rejects-invalid-entries"]

VerificationTest[
  AssociationQ[Catch[CertiField`Private`buildMinorantTable[60, 50, 1.*^-6], "CertiFieldTable"]],
  True, TestID -> "small-minorant-table-certifies"]

VerificationTest[
  Catch[CertiField`Private`rigorousIntegralFn[Function[x, Sin[x]], Function[x, Sin[x]], {0, 1}, "Tolerance" -> 10^-20],
    "CertiField"] /. Interval[{a_, b_}] :> (a <= 1 - Cos[1] <= b && b - a < 10^-19),
  True, TestID -> "rigorous-quadrature-encloses-integral"]

(* ---------- one-dimensional models: exact minima are enclosed ---------- *)

phi4 = CertifiedVariationalBound["Phi4Kink"];

VerificationTest[phi4["Status"], "Certified", TestID -> "phi4-status"]
VerificationTest[phi4["LowerBound"] <= 2 Sqrt[2]/3 <= phi4["UpperBound"], True, TestID -> "phi4-encloses-exact-minimum"]
VerificationTest[phi4["RelativeGap"] < 10^-12, True, TestID -> "phi4-gap-tight"]

VerificationTest[
  With[{c = CertifiedVariationalBound["SineGordonKink"]}, c["LowerBound"] <= 8 <= c["UpperBound"] && c["RelativeGap"] < 10^-11],
  True, TestID -> "sine-gordon-encloses-8"]

VerificationTest[
  With[{c = CertifiedVariationalBound["MagneticDomainWall", "A" -> 2, "K" -> 1/2]}, c["LowerBound"] <= 4 <= c["UpperBound"]],
  True, TestID -> "domain-wall-encloses-4-sqrt-AK"]

(* a poor candidate gives an honest, wide gap; a gap tolerance turns it into Inconclusive *)
VerificationTest[
  With[{c = CertifiedVariationalBound["Phi4Kink", Function[x, Tanh[x]], "GapTolerance" -> 10^-3]},
    {c["Status"], c["LowerBound"] <= 2 Sqrt[2]/3 <= c["UpperBound"], c["RelativeGap"] > 0.05}],
  {"Inconclusive", True, True}, TestID -> "rough-candidate-inconclusive"]

(* finite domain, candidate from NDSolve *)
VerificationTest[
  Module[{u, x, sol, c},
    sol = NDSolveValue[{u''[x] == u[x]^3 - u[x], u[0] == -1, u[3] == 1/2}, u, {x, 0, 3}];
    c = CertifiedVariationalBound[CertiFieldModel["ScalarField1D",
      <|"Potential" -> Function[s, (s^2 - 1)^2/4], "Domain" -> {0, 3}, "BoundaryValues" -> {-1, 1/2}|>], sol];
    c["Status"] === "Certified" && c["RelativeGap"] < 10^-9],
  True, TestID -> "finite-domain-ndsolve-candidate"]

(* ---------- verification and tamper detection ---------- *)

VerificationTest[VerifyCertificate[phi4], True, TestID -> "verify-genuine-certificate"]
VerificationTest[VerifyCertificate[CertifiedBound[Append[phi4["Association"], "LowerBound" -> 1.]]], False,
  TestID -> "verify-rejects-tampered-lower-bound"]

(* ---------- errors are reported as Failure objects ---------- *)

VerificationTest[FailureQ[CertifiedVariationalBound["NoSuchModel"]], True,
  {CertifiedVariationalBound::nomodel}, TestID -> "unknown-model"]
VerificationTest[
  CertifiedVariationalBound[CertiFieldModel["ScalarField1D", <|"Potential" -> Function[s, Abs[s^2 - 1]], "BoundaryValues" -> {-1, 1}|>],
    Function[x, Tanh[x]]] // FailureQ,
  True, TestID -> "non-smooth-potential-refused"]
VerificationTest[
  CertifiedVariationalBound[CertiFieldModel["ScalarField1D", <|"Potential" -> Function[s, (s^2 - 1)^2], "BoundaryValues" -> {0, 1}|>],
    Function[x, Tanh[x]]] // FailureQ,
  True, TestID -> "non-vacuum-boundary-value-refused"]
VerificationTest[FailureQ[CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> 2]], True, TestID -> "unsupported-pion-mass"]

(* ---------- presentation ---------- *)

VerificationTest[StringQ[phi4["Summary"]] && StringContainsQ[phi4["Summary"], "Certified"], True, TestID -> "summary-string"]
VerificationTest[Head[ToBoxes[phi4]], InterpretationBox | TagBox | TemplateBox, SameTest -> MatchQ, TestID -> "summary-box-renders"]
