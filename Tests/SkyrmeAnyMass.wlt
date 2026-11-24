(* ::Package:: *)

(* SkyrmeDegreeOne for a pion mass that is not one of the four shipped ones (docs/SKYRME_ANY_MASS.md).
   About eight minutes: two certificates and three verifications at Quick resolution.
   Run with:  wolframscript -file tools/run_tests.wls --skyrme-any-mass *)

Needs["CertiField`"];

(* ---------- the heuristics against the original construction (no certified number depends on these) ---------- *)

(* the tangent point of the envelope: closed-form equation against the bisection of Minorants.wl *)
VerificationTest[
  Module[{k1 = (3 - Sqrt[5.])/2, k2 = (3 + Sqrt[5.])/2, closed},
    closed[k_] := If[k1 <= k <= k2, 0., Module[{lo = 0., hi = Sqrt[1. + k], mid},
        Do[mid = 0.5 (lo + hi); If[mid mid - 1. - k + Sqrt[k (4. mid + 5.)] < 0., lo = mid, hi = mid], {80}]; 0.5 (lo + hi)]];
    Max[Abs[closed[#] - CertiField`Private`dTangent[#]] & /@ {1.*^-8, 0.001, 0.1, 0.3, 0.38, 1., 2.7, 5., 40., 300.}] < 10^-9],
  True, TestID -> "any-mass-tangent-point"]

(* the right-hand side on the three branches, against direct formulas *)
VerificationTest[
  Module[{h, hD, hK, tau = 1.2, r = 0.5, kap = 30., mu2 = 0.49, d, v, out},
    h[k_, x_] := 2 Sqrt[(1 + 2 x) (x^2 + 2 x + k)]; hK[k_, x_] := Sqrt[(1 + 2 x)/(x^2 + 2 x + k)];
    hD[k_, x_] := (2 (x^2 + 2 x + k) + (1 + 2 x) (2 x + 2))/Sqrt[(1 + 2 x) (x^2 + 2 x + k)];
    d = Sin[tau]^2/r^2; v = 2 mu2 (1 + Cos[tau]);
    out = CertiField`Private`skRhs[tau, {r, kap, 0.}, mu2, True];
    Max[Abs[out - {hK[kap, d], -2 mu2 Sin[tau] - 2/r (h[kap, d] - d hD[kap, d]), 4 Pi r^2 (h[kap, d] - (kap - v) hK[kap, d])}]] < 10^-10],
  True, TestID -> "any-mass-extremal-system-hedgehog-branch"]

VerificationTest[
  Module[{h, hK, tau = 2.9, r = 3., kap = 0.05, mu2 = 0., d, dt, lam, out},
    h[k_, x_] := 2 Sqrt[(1 + 2 x) (x^2 + 2 x + k)]; hK[k_, x_] := Sqrt[(1 + 2 x)/(x^2 + 2 x + k)];
    d = Sin[tau]^2/r^2; dt = CertiField`Private`dTangent[kap]; lam = d/dt;
    out = CertiField`Private`skRhs[tau, {r, kap, 0.}, mu2, True];
    d < dt && Max[Abs[out[[;; 2]] - {(1 - lam)/Sqrt[kap] + lam hK[kap, dt], -2/r 2 Sqrt[kap]}]] < 10^-9],
  True, TestID -> "any-mass-extremal-system-two-phase-branch"]

(* critical slopes of the relaxed problem and of the hedgehog at mu = 1: 2.5729512863 and 2.5741735436 in the original *)
VerificationTest[
  {Abs[CertiField`Private`skCriticalSlope[1, True] - 2.5729512863262576] < 10^-7,
   Abs[CertiField`Private`skCriticalSlope[1, False] - 2.5741735436336066] < 10^-7},
  {True, True}, TestID -> "any-mass-critical-slopes-match-original"]

(* the grid generated for a shipped mass is the shipped grid, up to the accuracy of the differential-equation solver *)
VerificationTest[
  Module[{g = CertiField`Private`skyrmeGeneratedGrid[1, "Quick"], s = CertiField`Private`loadSkyrmeGrid[CertiField`Private`skyrmeGridDir[1, "Quick"]]},
    {Length[g["taus"]] == Length[s["taus"]], Normal[g["offs"]] === Normal[s["offs"]], Max[Abs[g["taus"] - s["taus"]]] < 10^-7,
     Max[Abs[Log[g["vref"]/s["vref"]]]] < 10^-6, Max[Abs[g["eug"]/s["eug"] - 1]] < 10^-14, CertiField`Private`skyrmeGridQ[g]}],
  {True, True, True, True, True, True}, TestID -> "any-mass-generated-grid-matches-shipped-grid"]

(* monotone interpolation: reproduces a cubic's monotone data at the nodes and stays monotone in between *)
VerificationTest[
  Module[{x = {0., 1., 2., 4., 7.}, y = {0., 1., 1.5, 1.6, 5.}, d, v},
    d = CertiField`Private`skPchipSlopes[x, y];
    v = First[CertiField`Private`skPchipEval[x, y, d, Subdivide[0., 7., 200]]];
    {Max[Abs[First[CertiField`Private`skPchipEval[x, y, d, x]] - y]] < 10^-14, Min[Differences[v]] >= 0.}],
  {True, True}, TestID -> "any-mass-monotone-interpolation"]

(* ---------- a mass without shipped data ---------- *)

am = CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> 3/4];

VerificationTest[am["Status"], "Certified", TestID -> "any-mass-status"]
VerificationTest[am["Certificate"]["GridSource"], "Generated", TestID -> "any-mass-grid-is-generated"]
(* between the certified lower bound at 13/20 (1.33062, Quick) and the hedgehog bound at 1 (1.41623): Lemma B *)
VerificationTest[1.33062 < am["LowerBound"] <= am["UpperBound"] < 1.41623, True, TestID -> "any-mass-between-neighbouring-masses"]
VerificationTest[am["RelativeGap"] < 0.004, True, TestID -> "any-mass-gap"]
VerificationTest[VerifyCertificate[am], True, TestID -> "any-mass-verifies"]

(* the verification is a recomputation: a raised lower bound, a lower bound taken from a larger mass, and a grid that is
   not increasing are all rejected *)
VerificationTest[VerifyCertificate[CertifiedBound[Append[am["Association"], "LowerBound" -> am["LowerBound"] + 10.^-6]]], False,
  TestID -> "any-mass-verify-rejects-tampered-lower-bound"]
VerificationTest[
  VerifyCertificate[CertifiedBound[Append[am["Association"], "Certificate" -> Append[am["Certificate"], "LowerBoundMass" -> 1]]]], False,
  TestID -> "any-mass-verify-rejects-lower-bound-from-a-larger-mass"]
VerificationTest[
  Module[{c = am["Certificate"], g},
    g = c["Grid"]; g["eug"] = Reverse[g["eug"]];
    VerifyCertificate[CertifiedBound[Append[am["Association"], "Certificate" -> Append[c, "Grid" -> g]]]]], False,
  TestID -> "any-mass-verify-rejects-invalid-grid"]

(* ---------- the monotonicity fallback ---------- *)

mono = CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> 3/4, Method -> "Monotonicity"];

VerificationTest[{mono["Status"], mono["Certificate"]["LowerBoundMass"], mono["Certificate"]["GridSource"]}, {"Certified", 13/20, "Monotonicity"},
  TestID -> "any-mass-monotone-uses-largest-shipped-mass-below"]
(* it is the Quick lower bound of the mass 13/20 on its shipped grid *)
VerificationTest[Abs[mono["LowerBound"] - 1.3306232116384036] < 10^-12, True, TestID -> "any-mass-monotone-bound-is-the-shipped-one"]
VerificationTest[mono["LowerBound"] <= am["LowerBound"] && mono["UpperBound"] == am["UpperBound"], True, TestID -> "any-mass-monotone-bound-is-weaker"]
VerificationTest[StringContainsQ[ToString[mono["Theorem"]], "13/20"] && ! StringContainsQ[ToString[am["Theorem"]], "13/20"], True,
  TestID -> "any-mass-theorem-names-the-mass-of-the-lower-bound"]
VerificationTest[VerifyCertificate[mono], True, TestID -> "any-mass-monotone-verifies"]

(* ---------- failure modes ---------- *)

VerificationTest[FailureQ[CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> -1/2]], True, TestID -> "any-mass-negative-mass-refused"]
VerificationTest[FailureQ[CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> 5/2]], True, TestID -> "any-mass-large-mass-refused"]
VerificationTest[FailureQ[CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> Sqrt[2]/2]], True, TestID -> "any-mass-irrational-mass-refused"]
VerificationTest[FailureQ[CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> 3/4, Method -> "Nonsense"]], True, TestID -> "any-mass-unknown-method-refused"]
