(* ::Package:: *)

(* Extensions of ScalarFieldND (docs/SCALAR_ND_2.md): source term, prescribed mean, periodic directions, potentials on a
   value range.  Run with:  wolframscript -file tools/run_tests.wls --nd2
   Known minima must lie inside the certified intervals; the inequalities are never relaxed. *)

Needs["CertiField`"];

nd2Model[p_] := CertiFieldModel["ScalarFieldND", p];
nd2Raw[c_CertifiedBound] := Module[{a = c[[1]]},
  CertiField`Private`nd2LowerBound[CertiField`Private`nd2Spec[a["Parameters"]], a["Certificate", "Field"], "Subintervals" -> a["Certificate", "Subintervals"]]];

(* ---------- potentials without growth: E = Integrate[k |grad u|^2 - f u], k = 1/2 ---------- *)

(* f = 4k = 2, u = 0 on x1 = 0 and x1 = 1: u* = x1 (1 - x1), min E = -k/3 per unit cross-section; the flux of the
   minimiser is a Raviart-Thomas field, so the lower bound is exact up to rounding *)
nd2Src = CertifiedVariationalBound[nd2Model[<|"Potential" -> Function[u, 0], "Domain" -> {{0, 1}, {0, 1/4}}, "Cells" -> {8, 2},
  "Source" -> 2, "NaturalFaces" -> {{2, "Lower"}, {2, "Upper"}}|>]];
VerificationTest[nd2Src["LowerBound"] <= -1/24 <= nd2Src["UpperBound"], True, TestID -> "nd2-constant-source-encloses-minimum"]
VerificationTest[Abs[nd2Raw[nd2Src] + 1/24] < 10^-14, True, TestID -> "nd2-constant-source-lower-bound-exact-to-rounding"]
VerificationTest[VerifyCertificate[nd2Src], True, TestID -> "nd2-constant-source-certificate-verifies"]

VerificationTest[
  With[{c = CertifiedVariationalBound[nd2Model[<|"Potential" -> Function[u, 0], "Domain" -> {{0, 1}, {0, 1/4}}, "Cells" -> {8, 2},
      "Source" -> 2, "PeriodicDirections" -> {2}|>]]},
    {c["LowerBound"] <= -1/24 <= c["UpperBound"], Abs[nd2Raw[c] + 1/24] < 10^-14, VerifyCertificate[c]}],
  {True, True, True}, TestID -> "nd2-constant-source-periodic-sides"]

VerificationTest[
  With[{c = CertifiedVariationalBound[nd2Model[<|"Potential" -> Function[u, 0], "Domain" -> {{0, 1}, {0, 1/4}, {0, 1/4}}, "Cells" -> {8, 2, 2},
      "Source" -> 2, "NaturalFaces" -> {{2, "Lower"}, {2, "Upper"}}, "PeriodicDirections" -> {3}|>]]},
    {c["LowerBound"] <= -1/96 <= c["UpperBound"], Abs[nd2Raw[c] + 1/96] < 10^-14}],
  {True, True}, TestID -> "nd2-constant-source-3d"]

(* f = c x1 with c = 3/2: u* = c x1 (1 - x1^2)/(12 k), min E = -c^2/(180 k) = -1/40 *)
nd2Lin = CertifiedVariationalBound[nd2Model[<|"Potential" -> Function[u, 0], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> {16, 16},
  "Source" -> Function[{x, y}, 3/2 x], "NaturalFaces" -> {{2, "Lower"}, {2, "Upper"}}|>]];
VerificationTest[{nd2Lin["LowerBound"] <= -1/40 <= nd2Lin["UpperBound"], nd2Lin["RelativeGap"] < 0.02}, {True, True},
  TestID -> "nd2-linear-source-encloses-minimum"]

(* the mass-constrained strip: E = Integrate[k |grad u|^2], u = 0 on x1 = 0 and x1 = 1, mean m: u* = 6 m x1 (1 - x1),
   min E = 12 k m^2 per unit cross-section; m = 1/3, cross-section 1/4: 1/6 *)
nd2Strip = CertifiedVariationalBound[nd2Model[<|"Potential" -> Function[u, 0], "Domain" -> {{0, 1}, {0, 1/4}}, "Cells" -> {16, 2},
  "Mean" -> 1/3, "NaturalFaces" -> {{2, "Lower"}, {2, "Upper"}}|>]];
VerificationTest[{nd2Strip["LowerBound"] <= 1/6 <= nd2Strip["UpperBound"], Abs[nd2Raw[nd2Strip] - 1/6] < 10^-13, nd2Strip["RelativeGap"] < 0.01,
    VerifyCertificate[nd2Strip]}, {True, True, True, True}, TestID -> "nd2-mass-constrained-strip"]

(* ---------- the lower bound against an exact evaluation of its formula on a tiny grid ---------- *)

(* the right-hand side of Theorem 8 with Propositions 10, Lemmas 3 and 11, for a two-dimensional grid, a spline field
   with K knot intervals, and data that are constant on each Dirichlet face; evaluated in exact (or 40-digit)
   arithmetic, with no code shared with nd2LowerBound except exactRange (the bound for W'' on a subinterval) *)
nd2Exact[spec_, f_, qs_, prec_] := Module[
  {kk = f["K"], t0, t1, dl, w, cf, lam = SetPrecision[f["Lambda"], Infinity], k = spec["k"], h = spec["h"], P = spec["P"], cells = spec["Cells"], vr = spec["VR"],
   c1 = spec["C1"], c2 = spec["C2"], lin = SetPrecision[spec["Lin"], Infinity], g = SetPrecision[spec["G"], Infinity], wq, pAt, sval, sd, s2, flux = 0, intM = 0, m},
  {t0, t1} = f["TRange"]; dl = (t1 - t0)/kk; w = dl/qs;
  cf = SetPrecision[f["C"], Infinity];
  pAt[p_, t_] := Fold[#1 t + #2 &, 0, Reverse[p]];
  wq[t_] := N[pAt[P, t] + If[c1 == 0, 0, c1 (t - vr[[1]]) If[t == vr[[1]], 0, Log[t - vr[[1]]]]] +
    If[c2 == 0, 0, c2 (vr[[2]] - t) If[t == vr[[2]], 0, Log[vr[[2]] - t]]], prec];
  sval[c_, t_] := Module[{z = Clip[(t - t0)/dl, {0, kk}], kap, th}, kap = Min[Floor[z], kk - 1]; th = z - kap;
    c[[kap + 1]] (1 - th)^2/2 + c[[kap + 2]] (-2 th^2 + 2 th + 1)/2 + c[[kap + 3]] th^2/2];
  sd[c_, kap_, th_] := ((c[[kap + 2]] - c[[kap + 1]]) (1 - th) + (c[[kap + 3]] - c[[kap + 2]]) th)/dl;
  s2[c_, kap_] := (c[[kap + 1]] - 2 c[[kap + 2]] + c[[kap + 3]])/dl^2;
  Do[
    Module[{faces = {{cf[[1, a, b]], cf[[1, a + 1, b]]}, {cf[[2, a, b]], cf[[2, a, b + 1]]}}},
      Do[
        m = -Infinity;
        Do[
          Module[{ph, d2, h2, ca = 0, cb = 0, ta, tb},
            ph[t_, th_] := Sum[(sval[faces[[i, 2]], t] - sval[faces[[i, 1]], t])/h[[i]], {i, 2}] +
              Sum[sd[faces[[i, v[[i]] + 1]], kap, th]^2, {i, 2}]/(4 k) - wq[t] + (lam + lin[[a + v[[1]], b + v[[2]]]]) t;
            d2 = Sum[(s2[faces[[i, 2]], kap] - s2[faces[[i, 1]], kap])/h[[i]], {i, 2}];
            ta = t0 + (kap + (j - 1)/qs) dl; tb = t0 + (kap + j/qs) dl;
            h2 = Last[CertiField`Private`exactRange[CertiField`Private`polyD[CertiField`Private`polyD[P]], ta, tb]];
            If[c1 != 0, If[ta > vr[[1]], h2 += c1/(ta - vr[[1]]), ca = c1 w]];
            If[c2 != 0, If[tb < vr[[2]], h2 += c2/(vr[[2]] - tb), cb = c2 w]];
            m = Max[m, Max[ph[ta, (j - 1)/qs] + ca, ph[tb, j/qs] + cb] + w^2/8 Max[0, h2 - d2]]],
          {kap, 0, kk - 1}, {j, 1, qs}];
        intM += Times @@ h/4 m,
        {v, Tuples[{0, 1}, 2]}]],
    {a, cells[[1]]}, {b, cells[[2]]}];
  (* the Dirichlet faces x1 = const, where the data are constant *)
  Do[flux += -h[[2]] sval[cf[[1, 1, b]], g[[1, b]]] + h[[2]] sval[cf[[1, -1, b]], g[[-1, b]]], {b, cells[[2]]}];
  flux - intM + lam spec["Mean"] spec["Volume"]];

nd2TinyParams = <|"Potential" -> Function[u, u^4/4 - u^2/2], "ValueRange" -> {-3/2, 3/2}, "Source" -> Function[{x, y}, y - x/2], "Mean" -> 1/5,
  "PeriodicDirections" -> {2}, "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> {2, 2}, "BoundaryData" -> Function[{x, y}, 1/2 - x]|>;
nd2TinyField[spec_, kk_, tr_, seed_] := BlockRandom[SeedRandom[seed];
  Module[{c1 = RandomReal[{-1, 1}, {3, 2, kk + 2}], c2 = RandomReal[{-1, 1}, {2, 3, kk + 2}]},
    c2[[All, 3]] = c2[[All, 1]];
    <|"TRange" -> tr, "K" -> kk, "Affine" -> False, "S0" -> {ConstantArray[0., {3, 2}], ConstantArray[0., {2, 3}]}, "C" -> {c1, c2}, "Lambda" -> 0.3|>]];

VerificationTest[
  Module[{spec = CertiField`Private`nd2Spec[nd2TinyParams], f, lb, ex},
    f = nd2TinyField[spec, 2, {-3/2, 3/2}, 7];
    lb = CertiField`Private`nd2LowerBound[spec, f, "Subintervals" -> 4];
    ex = nd2Exact[spec, f, 4, Infinity];
    {lb <= ex, ex - lb < 10^-10}],
  {True, True}, TestID -> "nd2-lower-bound-equals-exact-formula-polynomial"]

VerificationTest[
  Module[{spec = CertiField`Private`nd2Spec[Join[nd2TinyParams, <|"Potential" -> Function[u, u Log[u] + (1 - u) Log[1 - u] + 3 u (1 - u)],
      "ValueRange" -> {0, 1}, "BoundaryData" -> Function[{x, y}, 1/2 - x/4]|>]], f, lb, ex},
    f = nd2TinyField[spec, 2, {0, 1}, 11];
    lb = CertiField`Private`nd2LowerBound[spec, f, "Subintervals" -> 4];
    ex = nd2Exact[spec, f, 4, 40];
    {lb <= ex, ex - lb < 10^-10}],
  {True, True}, TestID -> "nd2-lower-bound-equals-exact-formula-logarithmic"]

VerificationTest[
  Module[{x = Developer`ToPackedArray[{1., 0.5, 2., 3., 0.1, 1.*^-5, 7.3*^10, 1.*^-310, 0.7071067811865476, 1.4142135623730951}], e},
    e = CertiField`Private`nd2LnPt[x];
    And @@ MapThread[SetPrecision[#1, Infinity] <= Log[SetPrecision[#3, Infinity]] <= SetPrecision[#2, Infinity] &, {e[[1]], e[[2]], x}]],
  True, TestID -> "nd2-logarithm-encloses"]

(* ---------- the logarithmic free energy with a prescribed mean: certified, falsification, tampering ---------- *)
(* the two phases held at the ends, mean 1/2: one interface; minimum about 0.0277 (one-dimensional floating-point
   computation, for orientation only) *)

nd2PhaseParams = <|"Potential" -> Function[u, u Log[u] + (1 - u) Log[1 - u] + 3 u (1 - u)], "ValueRange" -> {0, 1}, "Mean" -> 1/2,
  "Domain" -> {{0, 4}, {0, 1/2}}, "Cells" -> {8, 2}, "BoundaryData" -> Function[{x, y}, If[x < 2, 1/10, 9/10]], "PeriodicDirections" -> {2}|>;
nd2Phase = CertifiedVariationalBound[nd2Model[nd2PhaseParams], Method -> {"Calibration", "KnotIntervals" -> 8}];

VerificationTest[{nd2Phase["Status"], nd2Phase["LowerBound"] <= nd2Phase["UpperBound"], VerifyCertificate[nd2Phase]}, {"Certified", True, True},
  TestID -> "nd2-logarithmic-potential-certified"]

(* no admissible field may give a lower bound above the certified upper bound, and no admissible candidate an energy
   below the certified lower bound *)
VerificationTest[
  Module[{spec = CertiField`Private`nd2Spec[nd2PhaseParams], f0 = nd2Phase[[1]]["Certificate", "Field"], lbs, ubs},
    lbs = Table[
      BlockRandom[SeedRandom[s];
        Module[{c = (# + RandomReal[{-0.3, 0.3}, Dimensions[#]]) & /@ f0["C"]},
          c[[2, All, -1]] = c[[2, All, 1]];                (* direction 2 is periodic: equal values on the identified faces *)
          CertiField`Private`nd2LowerBound[spec, Join[f0, <|"C" -> c, "Lambda" -> f0["Lambda"] + RandomReal[{-0.5, 0.5}]|>], "Subintervals" -> 16]]],
      {s, 12}];
    ubs = Table[
      BlockRandom[SeedRandom[100 + s];
        Module[{u = RandomReal[{0.2, 0.8}, {9, 3}]}, u[[All, -1]] = u[[All, 1]]; u[[1]] = ConstantArray[0.1, 3]; u[[-1]] = ConstantArray[0.9, 3];
          CertiField`Private`nd2UpperBound[spec, u][[2]]]],
      {s, 12}];
    {Max[lbs] <= nd2Phase["UpperBound"], Min[ubs] >= nd2Phase["LowerBound"]}],
  {True, True}, TestID -> "nd2-falsification"]

VerificationTest[
  {VerifyCertificate[CertifiedBound[Append[nd2Phase[[1]], "LowerBound" -> nd2Phase["LowerBound"] + 10^-3]]],
   VerifyCertificate[CertifiedBound[Append[nd2Phase[[1]], "UpperBound" -> nd2Phase["UpperBound"] - 10^-3]]]},
  {False, False}, TestID -> "nd2-verify-rejects-tampered-certificates"]

(* ---------- refusals ---------- *)

VerificationTest[FailureQ[CertifiedVariationalBound[nd2Model[<|"Potential" -> Function[u, 0], "Domain" -> {{0, 1}, {0, 1}}, "Source" -> 1,
  "NaturalFaces" -> {{1, "Lower"}, {1, "Upper"}, {2, "Lower"}, {2, "Upper"}}|>]]], True, TestID -> "nd2-unbounded-linear-energy-refused"]
VerificationTest[FailureQ[CertifiedVariationalBound[nd2Model[<|"Potential" -> Function[u, u Log[u]], "Domain" -> {{0, 1}, {0, 1}}, "Mean" -> 1/2|>]]],
  True, TestID -> "nd2-logarithm-without-value-range-refused"]
VerificationTest[FailureQ[CertifiedVariationalBound[nd2Model[<|"Potential" -> Function[u, u^2], "Domain" -> {{0, 1}, {0, 1}}, "PeriodicDirections" -> {1},
  "BoundaryData" -> Function[{x, y}, x]|>]]], True, TestID -> "nd2-non-periodic-data-refused"]
