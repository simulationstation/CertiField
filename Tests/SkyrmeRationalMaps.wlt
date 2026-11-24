(* ::Package:: *)

(* SkyrmeRationalMap: certified upper bounds from the rational-map ansatz (docs/SKYRME_RATIONAL_MAPS.md).
   A few minutes.  Run with:  wolframscript -file tools/run_tests.wls --skyrme-rational-maps *)

Needs["CertiField`"];

(* the interval arithmetic on second-order Taylor data encloses the exact derivatives of the angular integrand *)
VerificationTest[
  Module[{z = \[FormalZ], p, q, rho, t, zz, P, Q, W, WW, g, ex, jet},
    p = z^2 - I z/3; q = 1 + z/2;
    zz = rho (1 + I t)/(1 - I t);
    W = Expand[D[p, z] q - p D[q, z]];
    P = Together[(1 - I t)^2 (p /. z -> zz)]; Q = Together[(1 - I t)^2 (q /. z -> zz)]; WW = Together[(1 - I t)^2 (W /. z -> zz)];
    g = 2 rho (1 + rho^2)^2 (1 + t^2)^3 ComplexExpand[WW Conjugate[WW]]^2/ComplexExpand[P Conjugate[P] + Q Conjugate[Q]]^4;
    ex = Flatten[Table[{g, D[g, {rho, 2}], D[g, {t, 2}]} /. {rho -> 3/8 + i/256, t -> 1/4 + j/256}, {i, 0, 4}, {j, 0, 4}], 1];
    jet = First[CertiField`Private`rmEvalG[CertiField`Private`rmPieces[p, q, z, 2][[1]],
      {{{0.375}, {0.390625}}, {{1.}, {1.}}, {{0.}, {0.}}, {{0.}, {0.}}, {{0.}, {0.}}, {{0.}, {0.}}},
      {{{0.25}, {0.265625}}, {{0.}, {0.}}, {{1.}, {1.}}, {{0.}, {0.}}, {{0.}, {0.}}, {{0.}, {0.}}}]];
    And @@ Flatten[Table[SetPrecision[jet[[k, 1, 1]], Infinity] <= e[[i]] <= SetPrecision[jet[[k, 2, 1]], Infinity], {e, ex}, {i, 3}, {k, {{1, 4, 6}[[i]]}}]]],
  True, TestID -> "rational-map-jets-enclose-exact-derivatives"]

(* known angular integrals: 1 for R = z, Pi + 8/3 for R = z^2 (Houghton, Manton, Sutcliffe) *)
VerificationTest[
  Module[{a = CertiField`Private`rationalMapIntegral[\[FormalZ], 1, \[FormalZ], 1], b = CertiField`Private`rationalMapIntegral[\[FormalZ]^2, 1, \[FormalZ], 2]},
    {a["Lower"] <= 1 <= a["Upper"], a["Upper"] - a["Lower"] < 10^-4,
     b["Lower"] <= CertiField`Private`$piLo + 8/3 && CertiField`Private`$piHi + 8/3 <= b["Upper"], b["Upper"] - b["Lower"] < 10^-3}],
  {True, True, True, True}, TestID -> "rational-map-known-angular-integrals"]

(* degree one with II = 1 is exactly the hedgehog computation of SkyrmeDegreeOne *)
VerificationTest[
  Module[{r = {0, 1/4, 1/2, 3/4, 1, 3/2, 2, 3, 5, 8}, rat, tr},
    rat[x_] := Rationalize[SetPrecision[x, 22], 10^-18];
    tr = <|"r" -> r, "p" -> Join[{0}, rat /@ N[Pi - 4 ArcTan[Exp[-3/2 #]] & /@ r[[2 ;; -2]], 22]],
      "dp" -> (rat /@ N[6 Exp[-3/2 #]/(1 + Exp[-3 #]) & /@ Most[r], 22])|>;
    CertiField`Private`rationalMapUpperBound[tr, 1, 1, 1/2]["UpperBound"] == CertiField`Private`hedgehogUpperBound[tr, 1/2]["UpperBound"]],
  True, TestID -> "rational-map-degree-one-is-the-hedgehog"]

(* B = 1: the hedgehog energy 1.2314457 (massless) from above *)
rm1 = CertifiedVariationalBound["SkyrmeRationalMap", "Degree" -> 1];
VerificationTest[{rm1["Status"], rm1["LowerBound"], 1.2314452 < rm1["UpperBound"] < 1.2316}, {"Certified", 1., True},
  TestID -> "rational-map-degree-one"]

(* B = 2, R = z^2: 1.208 per unit of degree in Houghton, Manton, Sutcliffe, Table 1 *)
rm2 = CertifiedVariationalBound["SkyrmeRationalMap", "Degree" -> 2];
VerificationTest[{rm2["Status"], rm2["LowerBound"], 1.2075 <= rm2["UpperBound"]/2 <= 1.209}, {"Certified", 2., True},
  TestID -> "rational-map-degree-two-published-value"]
VerificationTest[VerifyCertificate[rm2], True, TestID -> "rational-map-verifies"]
VerificationTest[VerifyCertificate[CertifiedBound[Append[rm2["Association"], "UpperBound" -> rm2["UpperBound"] - 10.^-6]]], False,
  TestID -> "rational-map-verify-rejects-tampered-upper-bound"]

(* failure modes *)
VerificationTest[FailureQ[CertifiedVariationalBound["SkyrmeRationalMap", "Degree" -> 0]], True, TestID -> "rational-map-degree-zero-refused"]
VerificationTest[FailureQ[CertifiedVariationalBound["SkyrmeRationalMap", "Degree" -> 2, "RationalMap" -> Function[z, (z^2 - 1)/(z - 1)]]], True,
  TestID -> "rational-map-wrong-degree-refused"]
VerificationTest[FailureQ[CertifiedVariationalBound["SkyrmeRationalMap", "Degree" -> 2, "RationalMap" -> Function[z, Sqrt[2] z^2]]], True,
  TestID -> "rational-map-inexact-coefficients-refused"]
