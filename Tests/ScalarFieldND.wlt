(* ::Package:: *)

(* ScalarFieldND suite (about two minutes).  Run with:  wolframscript -file tools/run_tests.wls --nd
   The expected numbers come from the Python reference (Reference/ScalarND).  The inequalities
   "lower <= reference value <= upper" must never be relaxed; tolerances on agreement with the reference are stated
   with their reason. *)

Needs["CertiField`"];

ndRef = FileNameJoin[{CertiField`Private`$CertiFieldDirectory, "Reference", "ScalarND"}];
ndTV = FileNameJoin[{ndRef, "testvectors"}];
ndIndex = Import[FileNameJoin[{ndTV, "index.json"}], "RawJSON"];
ndClose[a_, b_, tol_] := Abs[a - b] <= tol (1 + Abs[b]);
ndExact[x_] := SetPrecision[x, Infinity];

(* a test vector: the model built in the Wolfram Language, with the reference's boundary data, candidate and field *)
ndCase[name_String] := ndCase[name] = Module[{dir = FileNameJoin[{ndTV, name}], spec, g},
  spec = CertiField`Private`scalarFieldNDSpec[ToExpression[ndIndex[name, "wolfram_model"]]];
  g = CertiField`Private`ndImportArray[spec, dir, "gval.f64"];
  <|"SameData" -> Max[Abs[g - spec["G"]]] <= 10.^-14,
    "Spec" -> Append[spec, "G" -> g],
    "Field" -> CertiField`Private`ndImportField[spec, dir],
    "U" -> CertiField`Private`ndImportArray[spec, dir, "U.f64"],
    "Expected" -> ndIndex[name, "expected"]|>];

(* The exact flux through the Dirichlet faces, in rational arithmetic, for data that vary along at most one direction
   of each face (every 2D problem; in 3D, data such as g = x).  Independent of ndFluxLower: the clamped spline is
   integrated in closed form.  p^(t) = spline(clamp(t)); its antiderivative from tlo is piecewise cubic. *)
ndSplineValue[c_, kk_, {t0_, t1_}, t_] := Module[{z = Clip[(t - t0) kk/(t1 - t0), {0, kk}], kap, th},
  kap = Min[Floor[z], kk - 1]; th = z - kap;
  c[[kap + 1]] (1 - th)^2/2 + c[[kap + 2]] (-2 th^2 + 2 th + 1)/2 + c[[kap + 3]] th^2/2];
ndSplinePrimitive[c_, kk_, {t0_, t1_}, t_] := Module[{dl = (t1 - t0)/kk, z, kap, th, piece},
  piece[k_, x_] := dl (c[[k + 1]] (1 - (1 - x)^3)/6 + c[[k + 2]] (-2 x^3/3 + x^2 + x)/2 + c[[k + 3]] x^3/6);
  Which[
    t <= t0, (t - t0) ndSplineValue[c, kk, {t0, t1}, t0],
    t >= t1, Sum[piece[k, 1], {k, 0, kk - 1}] + (t - t1) ndSplineValue[c, kk, {t0, t1}, t1],
    True, z = (t - t0)/dl; kap = Min[Floor[z], kk - 1]; th = z - kap; Sum[piece[k, 1], {k, 0, kap - 1}] + piece[kap, th]]];
ndExactFlux[spec_, f_] := Module[{n = spec["n"], kk = f["K"], total = 0, vv, sign, area, s0, cb, nf, ga, gb, avg, cf},
  Do[
    If[spec["Dirichlet"][[i, side]],
      vv = Flatten /@ CertiField`Private`faceVertexValues[spec, i, side];
      sign = If[side == 1, -1, 1]; area = spec["CellVolume"]/spec["h"][[i]];
      s0 = Flatten[CertiField`Private`sliceAt[f["S0"][[i]], n, i, If[side == 1, 1, -1]]];
      nf = Length[s0];
      If[kk > 0, cb = ArrayReshape[CertiField`Private`sliceAt[f["C"][[i]], n, i, If[side == 1, 1, -1]], {nf, kk + 2}]];
      Do[
        (* the two values between which the datum varies linearly on this face *)
        {ga, gb} = Which[
          Length[vv] == 2, {vv[[1, j]], vv[[2, j]]},
          vv[[1, j]] == vv[[2, j]] && vv[[3, j]] == vv[[4, j]], {vv[[1, j]], vv[[3, j]]},
          vv[[1, j]] == vv[[3, j]] && vv[[2, j]] == vv[[4, j]], {vv[[1, j]], vv[[2, j]]},
          True, Throw[$Failed, "ndExactFlux"]];
        {ga, gb} = ndExact /@ {ga, gb};
        avg = ndExact[s0[[j]]] (ga + gb)/2;
        If[kk > 0,
          cf = ndExact[cb[[j]]];
          avg += If[ga == gb, ndSplineValue[cf, kk, f["TRange"], ga],
            (ndSplinePrimitive[cf, kk, f["TRange"], gb] - ndSplinePrimitive[cf, kk, f["TRange"], ga])/(gb - ga)]];
        total += sign area avg,
        {j, nf}]],
    {i, n}, {side, 2}];
  total];

(* ---------- primitives ---------- *)

VerificationTest[
  Module[{a = RandomReal[{-3, 3}, 2000], b, iv, ex = SetPrecision[#, Infinity] &},
    b = a + RandomReal[{0, 1}, 2000];
    iv = CertiField`Private`ivSqr[{a, b}];
    And @@ MapThread[ex[#3] <= ex[#1]^2 <= ex[#4] && ex[#3] <= ex[#2]^2 <= ex[#4] && (! (#1 <= 0 <= #2) || #3 == 0) &, {a, b, iv[[1]], iv[[2]]}]],
  True, TestID -> "nd-interval-square"]

VerificationTest[
  CertiField`Private`exactRange[{1, 0, -2, 0, 1}, -1/2, 3/4][[1]] <= 49/256 && CertiField`Private`exactRange[{1, 0, -2, 0, 1}, -1/2, 3/4][[2]] >= 1,
  True, TestID -> "nd-exact-range-encloses"]

(* restricted conjugates of W = (t^2-1)^2/4: rigorous upper bounds, equal to the reference to rounding *)
VerificationTest[
  Module[{p = ndIndex["primitives", "conjUpper"], W = {1/4, 0, -1/2, 0, 1/4}},
    {And @@ MapThread[ndClose[#1, #2, 10^-9] &, {CertiField`Private`conjUpper[W, N[p["d"]], 1, +1], p["plus_at_1"]}],
     And @@ MapThread[ndClose[#1, #2, 10^-9] &, {CertiField`Private`conjUpper[W, N[p["d"]], -1, -1], p["minus_at_m1"]}],
     And @@ MapThread[ndClose[#1, #2, 10^-9] &, {CertiField`Private`conjUpper[W, N[p["d"]], 0, +1], p["plus_at_0"]}]}],
  {True, True, True}, TestID -> "nd-conjugate-tables-match-reference"]

(* the conjugate of a quadratic is known exactly: W = t^2 - t, W*(d) = (d + 1)^2/4 = Max of the two half-line suprema *)
VerificationTest[
  Module[{d = {-2., -1., -0.3, 0., 0.5, 2.}, up},
    up = MapThread[Max, {CertiField`Private`conjUpper[{0, -1, 1}, d, 0, +1], CertiField`Private`conjUpper[{0, -1, 1}, d, 0, -1]}];
    And @@ MapThread[(SetPrecision[#1, Infinity] >= (SetPrecision[#2, Infinity] + 1)^2/4 && #1 - (#2 + 1)^2/4 < 10^-6) &, {up, d}]],
  True, TestID -> "nd-conjugate-of-quadratic"]

VerificationTest[
  Module[{p = ndIndex["primitives", "splineEval"], r},
    r = CertiField`Private`splineEval[N[p["C"]], p["K"], {-1, 1}, {N[p["t"]], N[p["t"]]}];
    And @@ MapThread[(#1 <= #4 + 10^-12 && #2 >= #3 - 10^-12 && #2 - #1 < 10^-10) &, {r[[1]], r[[2]], p["lo"], p["hi"]}]],
  True, TestID -> "nd-clamped-spline-matches-reference"]

(* the enclosure of the clamped spline contains its exact value (rational arithmetic), also at knots and clamps *)
VerificationTest[
  Module[{kk = 5, c, t, r},
    SeedRandom[3]; c = RandomReal[{-2, 2}, {40, kk + 2}];
    t = Join[RandomReal[{-1.4, 1.4}, 31], {-1., -0.6, 0.2, 1., 1.0000000000000002, -0.9999999999999999, 0., 3., -3.}];
    r = CertiField`Private`splineEval[c, kk, {-1, 1}, {t, t}];
    And @@ Table[ndExact[r[[1, j]]] <= ndSplineValue[ndExact[c[[j]]], kk, {-1, 1}, ndExact[t[[j]]]] <= ndExact[r[[2, j]]], {j, 40}]],
  True, TestID -> "nd-clamped-spline-encloses-exact-value"]

(* a machine overflow must stop the evaluation: the kernel would otherwise continue in software arithmetic, which the
   outward rounding does not cover *)
VerificationTest[
  With[{big = Developer`ToPackedArray[{1.*^200, 2.}], ok = Developer`ToPackedArray[{3., 2.}]},
    {FailureQ[Catch[CertiField`Private`ndMul[{big, big}, {big, big}], "CertiField"]],
     FailureQ[Catch[CertiField`Private`ndMul[{ok, ok}, {ok, ok}], "CertiField"]],
     FailureQ[CertifiedVariationalBound[CertiFieldModel["ScalarFieldND",
       <|"Potential" -> Function[u, u^4/4 + u^2/2 - u], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 8, "BoundaryData" -> 0|>],
       Function[{x, y}, 1.*^200 Sin[Pi x] Sin[Pi y]]]]}],
  {True, False, True}, TestID -> "nd-machine-overflow-is-refused"]

(* exact zeros (zero data, zero flux) lead to subnormal intermediate results; the bound must come out without messages.
   The minimum is 0; the bound is about -2 10^-7, the slack of the conjugate table (see PORTING.md, findings). *)
VerificationTest[
  Module[{spec = CertiField`Private`scalarFieldNDSpec[<|"Potential" -> Function[u, u^2], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 3, "BoundaryData" -> 0|>], f, lb},
    f = CertiField`Private`ndConvexField[spec, ConstantArray[0., {3, 3}]];
    lb = CertiField`Private`ndLowerBound[spec, f];
    -10^-6 < lb <= 0],
  True, TestID -> "nd-zero-data-no-underflow-messages"]

(* ---------- the rigorous evaluators against the reference, on identical inputs ---------- *)

Do[
  With[{name = nm},
    VerificationTest[ndCase[name]["SameData"], True, TestID -> "nd-" <> name <> "-boundary-data-agree"];
    VerificationTest[
      Module[{c = ndCase[name], ub},
        ub = CertiField`Private`ndUpperBound[c["Spec"], c["U"]];
        {ub[[1]] <= ub[[2]], ndClose[N[ub[[1]]], c["Expected", "upper_enclosure"][[1]], 10^-9], ndClose[N[ub[[2]]], c["Expected", "upper_enclosure"][[2]], 10^-9]}],
      {True, True, True}, TestID -> "nd-" <> name <> "-competitor-energy"]],
  {nm, {"convex2d", "convex3d", "wall2d", "ramp3d", "mixed2d"}}]

(* Boundary flux.  Where the datum is constant on each face, or the field has no spline part, there is no penalty
   and the port equals the reference to rounding. *)
Do[
  With[{name = nm},
    VerificationTest[
      Module[{c = ndCase[name]},
        ndClose[N[CertiField`Private`ndFluxLower[c["Spec"], c["Field"]]], c["Expected", "flux_lower"], 10^-9]],
      True, TestID -> "nd-" <> name <> "-boundary-flux"];
    VerificationTest[
      Module[{c = ndCase[name], lb},
        lb = CertiField`Private`ndLowerBound[c["Spec"], c["Field"], "Subintervals" -> c["Expected", "Subintervals"]];
        ndClose[N[lb], c["Expected", "lower"], If[c["Field"]["Affine"], 10^-6, 10^-9]]],
      True, TestID -> "nd-" <> name <> "-lower-bound"]],
  {nm, {"convex2d", "convex3d", "wall2d"}}]

(* Where the datum varies on the faces, the port takes the penalty constants of Lemma 5 over the knot intervals the
   datum meets, found by exact comparisons; the reference takes floor(...) -1/+1 in floating point, a larger set.
   The port is therefore at least as tight as the reference, and it must still lie below the exact flux. *)
Do[
  With[{name = nm},
    VerificationTest[
      Module[{c = ndCase[name], fl, ex},
        fl = CertiField`Private`ndFluxLower[c["Spec"], c["Field"]];
        ex = ndExactFlux[c["Spec"], c["Field"]];
        {fl <= ex, N[fl] >= c["Expected", "flux_lower"] - 10^-9 (1 + Abs[c["Expected", "flux_lower"]])}],
      {True, True}, TestID -> "nd-" <> name <> "-boundary-flux-below-exact-flux-and-above-reference"];
    (* the cell integral is unchanged: lower bound minus flux equals the reference's to rounding (10^-6 with an
       affine part, because the two implementations place the conjugate table slightly differently) *)
    VerificationTest[
      Module[{c = ndCase[name], lb, fl},
        fl = CertiField`Private`ndFluxLower[c["Spec"], c["Field"]];
        lb = CertiField`Private`ndLowerBound[c["Spec"], c["Field"], "Subintervals" -> c["Expected", "Subintervals"]];
        {ndClose[N[lb - fl], c["Expected", "lower"] - c["Expected", "flux_lower"], If[c["Field"]["Affine"], 10^-6, 10^-9]],
         N[lb] >= c["Expected", "lower"] - 10^-6 (1 + Abs[c["Expected", "lower"]])}],
      {True, True}, TestID -> "nd-" <> name <> "-lower-bound"]],
  {nm, {"ramp3d", "mixed2d"}}]

(* The same from random inputs: random grids, t-ranges, numbers of knot intervals, coefficients (affine and spline
   part) and nodal data, the data taken from the knots themselves (rounded to machine numbers), from the midpoints
   between them, from beyond the range and from random reals.  In 2D the data are arbitrary; in 3D they depend on
   the first coordinate only, so that the exact flux can be computed.  The certified value must never exceed it. *)
VerificationTest[
  Module[{ok = True, tight = 0, n, cells, kk, tr, pool, vals, g, spec, f, fl, ex},
    SeedRandom[20261123];
    Do[
      n = If[trial <= 30, 2, 3];
      cells = RandomInteger[{1, 4}, n];
      kk = RandomInteger[{1, 6}];
      tr = RandomChoice[{1, 1/2, 3/4}] {-1, 1} + RandomChoice[{0, 1/4, -1/3}];
      pool = Join[N[tr[[1]] + (tr[[2]] - tr[[1]]) Range[-2, 2 kk + 2]/(2 kk)], RandomReal[{-1.5, 1.5}, 8]];
      g = If[n == 2, RandomChoice[pool, cells + 1],
        With[{v = RandomChoice[pool, cells[[1]] + 1]}, Array[v[[#1]] &, cells + 1]]];
      spec = CertiField`Private`scalarFieldNDSpec[<|"Potential" -> Function[u, u^4 - u^2], "Domain" -> Table[{0, RandomChoice[{1, 2, 3/2}]}, n],
        "Cells" -> cells, "BoundaryData" -> g, "NaturalFaces" -> RandomChoice[{{}, {{n, "Upper"}}}]|>];
      f = <|"TRange" -> tr, "K" -> kk, "Affine" -> True,
        "S0" -> Table[RandomReal[{-1, 1}, CertiField`Private`faceShape[spec, i]], {i, n}],
        "C" -> Table[RandomReal[{-2, 2}, Append[CertiField`Private`faceShape[spec, i], kk + 2]], {i, n}]|>;
      fl = CertiField`Private`ndFluxLower[spec, f]; ex = ndExactFlux[spec, f];
      ok = ok && fl <= ex;
      If[ex - fl < 10^-12, tight++],
      {trial, 40}];
    {ok, tight < 40}],
  {True, True}, TestID -> "nd-boundary-flux-below-exact-flux-random-inputs"]

(* the flux bound is exact (to rounding) when the datum is constant on each face or the field is affine *)
VerificationTest[
  Table[
    With[{c = ndCase[name]},
      With[{fl = CertiField`Private`ndFluxLower[c["Spec"], c["Field"]], ex = ndExactFlux[c["Spec"], c["Field"]]},
        fl <= ex && ex - fl < 10^-12]],
    {name, {"convex2d", "wall2d"}}],
  {True, True}, TestID -> "nd-boundary-flux-exact-without-penalty"]

(* the finite-volume field built in the Wolfram Language is at least nearly as good as the reference's *)
VerificationTest[
  Module[{c = ndCase["convex2d"], p, f},
    p = CertiField`Private`ndFVSolve[c["Spec"], Mean[CertiField`Private`vertexPart[c["U"], 2, #] & /@ Tuples[{0, 1}, 2]]];
    f = CertiField`Private`ndConvexField[c["Spec"], p];
    N[CertiField`Private`ndLowerBound[c["Spec"], f]] >= c["Expected", "lower"] - 10^-6],
  True, TestID -> "nd-finite-volume-field"]

(* the reference's stored certificates, re-verified by this implementation: both bounds recomputed from the stored
   problem, candidate and field are at least as strong as the stored claims *)
VerificationTest[
  Table[
    With[{r = Catch[CertiField`Private`ndVerifyReferenceCertificate[FileNameJoin[{ndRef, "results", name}]], "CertiField"]},
      AssociationQ[r] && r["Verified"]],
    {name, {"quartic2d", "screened2d", "wall2d", "wall3d"}}],
  {True, True, True, True}, TestID -> "nd-reference-certificates-verified"]

(* ---------- the optimiser of calibration mode (heuristic part) ---------- *)

(* value and gradient of the smoothed functional at the reference's point *)
VerificationTest[
  Module[{o = ndIndex["optimiser"], spec, st, x, nc},
    spec = CertiField`Private`scalarFieldNDSpec[ToExpression[o["wolfram_model"]]];
    st = CertiField`Private`ndCalSetup[spec, {-1, 1}, o["K"], o["Qo"]];
    nc = Last[st["Off"]]; x = N[o["x"]];
    Table[
      With[{r = CertiField`Private`ndCalValueGradient[st, x[[;; nc]], cs["tau"]]},
        ndClose[r[[1]], cs["minus_J"], 10^-9] && Max[Abs[r[[2]] - cs["grad_minus_J"][[;; nc]]]] <= 10^-9 && Max[Abs[cs["grad_minus_J"][[nc + 1 ;;]]]] == 0],
      {cs, o["cases"]}]],
  {True, True, True}, TestID -> "nd-optimiser-objective-and-gradient-match-reference"]

(* the embeddings are exact: the same spline after doubling the knot intervals, the same divergence and the same
   normal components on the old faces after halving the cells *)
VerificationTest[
  Module[{kk = 4, c, f, g, t},
    SeedRandom[5]; c = RandomReal[{-1, 1}, {6, 3, kk + 2}];
    f = <|"TRange" -> {-1, 1}, "K" -> kk, "Affine" -> False, "S0" -> {ConstantArray[0., {6, 3}]}, "C" -> {c}|>;
    g = CertiField`Private`ndRefineT[f];
    t = RandomReal[{-1.2, 1.2}, 18];
    {g["K"], Dimensions[g["C"][[1]]],
     Max[Abs[Table[ndSplineValue[Flatten[c, 1][[j]], kk, {-1, 1}, t[[j]]] - ndSplineValue[Flatten[g["C"][[1]], 1][[j]], 2 kk, {-1, 1}, t[[j]]], {j, 18}]]] < 10^-13}],
  {8, {6, 3, 10}, True}, TestID -> "nd-refinement-in-t-is-exact"]

VerificationTest[
  Module[{sc, sf, f, g, div},
    sc = CertiField`Private`scalarFieldNDSpec[<|"Potential" -> Function[u, u^2], "Domain" -> {{0, 1}, {0, 2}}, "Cells" -> {3, 2}, "NaturalFaces" -> {{2, "Upper"}}|>];
    sf = CertiField`Private`scalarFieldNDSpec[<|"Potential" -> Function[u, u^2], "Domain" -> {{0, 1}, {0, 2}}, "Cells" -> {6, 4}, "NaturalFaces" -> {{2, "Upper"}}|>];
    SeedRandom[6];
    f = <|"TRange" -> {0, 0}, "K" -> 0, "Affine" -> True, "S0" -> {RandomReal[{-1, 1}, {4, 2}], RandomReal[{-1, 1}, {3, 3}]}, "C" -> {{}, {}}|>;
    f["S0"] = {f["S0"][[1]], f["S0"][[2]] ConstantArray[{1., 1., 0.}, 3]};
    g = CertiField`Private`ndProlongX[f, sf];
    div[s_, h_] := (s[[1]][[2 ;;, All]] - s[[1]][[;; -2, All]])/h[[1]] + (s[[2]][[All, 2 ;;]] - s[[2]][[All, ;; -2]])/h[[2]];
    {Dimensions /@ g["S0"],
     Max[Abs[div[g["S0"], N[sf["h"]]] - Flatten[Map[ConstantArray[#, {2, 2}] &, div[f["S0"], N[sc["h"]]], {2}], {{1, 3}, {2, 4}}]]] < 10^-13,
     Max[Abs[g["S0"][[1]][[1 ;; ;; 2, 1 ;; ;; 2]] - f["S0"][[1]]]] == 0.,
     CertiField`Private`ndRespectsNaturalFaces[sf, g]}],
  {{{7, 4}, {6, 5}}, True, True, True}, TestID -> "nd-prolongation-in-x-is-exact"]

(* ---------- public interface: known minima are enclosed ---------- *)
(* W = u^2 - u, u = 0 on the boundary: the minimum is a sine series; reference values -0.016025058075473 (unit square)
   and -0.009498616716544 (unit cube).  The Python reference gets gaps of 0.24 % (32^2) and 1.7 % (16^3). *)

ndSq = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND",
  <|"Potential" -> Function[u, u^2 - u], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 32, "BoundaryData" -> 0|>]];

VerificationTest[ndSq["Status"], "Certified", TestID -> "nd-square-status"]
VerificationTest[ndSq["LowerBound"] <= -0.016025058075473 <= ndSq["UpperBound"], True, TestID -> "nd-square-encloses-series-value"]
VerificationTest[ndSq["RelativeGap"] < 0.01, True, TestID -> "nd-square-gap"]
VerificationTest[ndSq["Certificate"]["Mode"], "Convex", TestID -> "nd-convex-potential-uses-convex-mode"]
VerificationTest[VerifyCertificate[ndSq], True, TestID -> "nd-verify-genuine-certificate"]
VerificationTest[VerifyCertificate[CertifiedBound[Append[ndSq["Association"], "LowerBound" -> -0.016]]], False,
  TestID -> "nd-verify-rejects-tampered-lower-bound"]
(* the datum is part of the certified problem: a certificate whose stored boundary values were changed is rejected *)
VerificationTest[
  With[{a = ndSq["Association"]},
    {VerifyCertificate[CertifiedBound[Append[a, "Certificate" -> Append[a["Certificate"], "BoundaryNodalValues" -> a["Certificate"]["BoundaryNodalValues"] + 0.001]]]],
     VerifyCertificate[CertifiedBound[Append[a, "Certificate" -> KeyDrop[a["Certificate"], "BoundaryNodalValues"]]]]}],
  {False, False}, TestID -> "nd-verify-rejects-tampered-boundary-values"]

VerificationTest[
  With[{c = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND",
      <|"Potential" -> Function[u, u^2 - u], "Domain" -> {{0, 1}, {0, 1}, {0, 1}}, "Cells" -> 16, "BoundaryData" -> 0|>]]},
    c["LowerBound"] <= -0.009498616716544 <= c["UpperBound"] && c["RelativeGap"] < 0.05],
  True, TestID -> "nd-cube-encloses-series-value"]

(* a user's own candidate: any function of the coordinates.  0.08 Sin[Pi x] Sin[Pi y] is the best multiple of the
   first sine mode to two digits; its energy is about 6 % above the minimum. *)
VerificationTest[
  With[{c = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND",
      <|"Potential" -> Function[u, u^4/4 + u^2/2 - u], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 32, "BoundaryData" -> 0|>],
      Function[{x, y}, 0.08 Sin[Pi x] Sin[Pi y]]]},
    c["Status"] === "Certified" && c["LowerBound"] <= c["UpperBound"] && 0.01 < c["RelativeGap"] < 0.5],
  True, TestID -> "nd-rough-candidate-gives-honest-wide-gap"]

(* a finite-element solution from NDSolve as the candidate: as good as the package's own *)
VerificationTest[
  Module[{u, x, y, sol, m, c},
    sol = NDSolveValue[{-Laplacian[u[x, y], {x, y}] + u[x, y]^3 + u[x, y] == 1, DirichletCondition[u[x, y] == 0, True]}, u, {x, 0, 1}, {y, 0, 1}];
    m = CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[s, s^4/4 + s^2/2 - s], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 32, "BoundaryData" -> 0|>];
    c = CertifiedVariationalBound[m, sol];
    c["Status"] === "Certified" && c["RelativeGap"] < 0.003 && c["UpperBound"] <= CertifiedVariationalBound[m]["UpperBound"] + 10^-6],
  True, TestID -> "nd-ndsolve-candidate"]

(* a very poor candidate still gives a valid certificate, and "RefineCandidate" turns it into a tight one *)
VerificationTest[
  Module[{m = CertiFieldModel["ScalarFieldND",
      <|"Potential" -> Function[u, u^4/4 + u^2/2 - u], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 32, "BoundaryData" -> 0|>], f = Function[{x, y}, 0.26 Sin[Pi x] Sin[Pi y]], a, b},
    a = CertifiedVariationalBound[m, f]; b = CertifiedVariationalBound[m, f, "RefineCandidate" -> True];
    {a["RelativeGap"] > 1, StringContainsQ[a["Summary"], "relative gap <= 125%"], b["RelativeGap"] < 0.01, b["UpperBound"] < a["UpperBound"], a["LowerBound"] <= b["UpperBound"]}],
  {True, True, True, True, True}, TestID -> "nd-refine-candidate"]

(* When W is not convex the energy has saddle points.  u = xy on the boundary of (-1,1)^2 with a double well: the
   symmetric critical point is a saddle, and plain Newton steps converge to it.  The package's own candidate and
   the candidate refined from an asymmetric start must both reach a minimiser (the two have the same energy;
   the Python reference gets 11.057755 on this grid), far below the function xy itself. *)
VerificationTest[
  Module[{g = Function[{x, y}, x y], spec, e},
    spec = CertiField`Private`scalarFieldNDSpec[<|"Potential" -> Function[u, 4 (u^2 - 1)^2], "Domain" -> {{-1, 1}, {-1, 1}}, "Cells" -> 16, "BoundaryData" -> g|>];
    e = N[CertiField`Private`ndUpperBound[spec, #][[2]]] & /@ {
      CertiField`Private`ndCompetitor[spec, Automatic, False, g],
      CertiField`Private`ndCompetitor[spec, Function[{x, y}, x y + 0.3 (1 - x^2) (1 - y^2)], True, None],
      CertiField`Private`ndCompetitor[spec, g, False, None]};
    {Abs[e[[1]] - e[[2]]] < 10^-6, Abs[e[[1]] - 11.057755] < 10^-5, e[[3]] > e[[1]] + 1}],
  {True, True, True}, TestID -> "nd-candidate-is-a-minimiser-not-a-saddle-point"]

(* only the boundary values of the boundary function are used: it may be singular inside the box *)
VerificationTest[
  With[{c = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND",
      <|"Potential" -> Function[u, u^2], "Domain" -> {{-1, 1}, {-1, 1}}, "Cells" -> 8, "BoundaryData" -> Function[{x, y}, Log[x^2 + y^2]]|>]]},
    ! FailureQ[c] && c["Status"] === "Certified" && 0 < c["LowerBound"] <= c["UpperBound"] && VerifyCertificate[c]],
  True, TestID -> "nd-boundary-function-singular-inside-the-box"]

(* grids with a single cell in some or all directions, in both modes; a constant potential, whose minimum is
   3 |box| = 6 exactly *)
VerificationTest[
  Module[{a, b, c},
    a = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND",
      <|"Potential" -> Function[u, u^2 - u], "Domain" -> {{0, 1}, {0, 1}, {0, 2}}, "Cells" -> {2, 1, 1}, "BoundaryData" -> Function[{x, y, z}, x y], "NaturalFaces" -> {{3, "Upper"}, {2, "Lower"}}|>]];
    b = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND",
      <|"Potential" -> Function[u, (u^2 - 1)^2], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 1, "BoundaryData" -> Function[{x, y}, x]|>], Method -> {"Calibration", "KnotIntervals" -> 1}];
    c = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, 3], "Domain" -> {{0, 1}, {0, 2}}, "Cells" -> 2, "BoundaryData" -> 1|>]];
    {a["LowerBound"] <= a["UpperBound"] && VerifyCertificate[a], b["LowerBound"] <= b["UpperBound"] && VerifyCertificate[b],
     c["LowerBound"] <= 6 <= c["UpperBound"] && c["RelativeGap"] < 10^-8 && VerifyCertificate[c]}],
  {True, True, True}, TestID -> "nd-single-cell-grids-and-constant-potential"]

(* a non-convex energy with a known minimum: the wall, with the field computed by the reference *)
VerificationTest[
  With[{c = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND",
      Append[ToExpression[ndIndex["wall2d", "wolfram_model"]], "Field" -> FileNameJoin[{ndTV, "wall2d"}]]],
      Function[{x, y}, Tanh[x/Sqrt[2]]]]},
    c["LowerBound"] <= 0.9449674266 <= c["UpperBound"] && c["LowerBound"] > 0.89],
  True, TestID -> "nd-wall-encloses-exact-value"]

(* the same wall with the field found by the Wolfram Language optimiser (calibration mode is chosen automatically
   because the potential is not convex); one level, as in the reference's test vector, whose bound is 0.8963 *)
ndWall = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", ToExpression[ndIndex["wall2d", "wolfram_model"]]],
  Method -> {"Calibration", "TRange" -> {-1, 1}, "Levels" -> {{{12, 2}, 4}}}];

VerificationTest[ndWall["LowerBound"] <= 0.9449674266 <= ndWall["UpperBound"], True, TestID -> "nd-calibration-wall-encloses-exact-value"]
VerificationTest[ndWall["LowerBound"] > 0.89 && ndWall["Certificate"]["Mode"] === "Calibration", True, TestID -> "nd-calibration-wall-as-tight-as-reference"]
VerificationTest[VerifyCertificate[ndWall], True, TestID -> "nd-calibration-certificate-verifies"]
VerificationTest[
  With[{a = ndWall["Association"]},
    VerifyCertificate[CertifiedBound[Append[a, "Certificate" -> Append[a["Certificate"],
      "Field" -> Append[a["Certificate"]["Field"], "C" -> 1.5 a["Certificate"]["Field"]["C"]]]]]]],
  False, TestID -> "nd-verify-rejects-tampered-field"]

(* automatic choices: mode, t-range and coarse-to-fine levels; three dimensions, data varying on the faces *)
VerificationTest[
  Module[{c = CertifiedVariationalBound[CertiFieldModel["ScalarFieldND",
      <|"Potential" -> Function[u, (u^2 - 1)^2], "Domain" -> {{-1, 1}, {-1, 1}, {-1, 1}}, "Cells" -> {4, 4, 2}, "BoundaryData" -> Function[{x, y, z}, x]|>],
      Method -> {"Calibration", "KnotIntervals" -> 4}]},
    {c["Status"], c["Certificate"]["Field"]["TRange"], c["LowerBound"] <= c["UpperBound"], c["LowerBound"] > 7.4, VerifyCertificate[c]}],
  {"Certified", {-1, 1}, True, True, True}, TestID -> "nd-calibration-3d-automatic-range"]

VerificationTest[
  {CertiField`Private`ndAutoLevels[{64, 64}, 32], CertiField`Private`ndAutoLevels[{96, 4}, 32], CertiField`Private`ndAutoLevels[{16, 16, 16}, 16]},
  {{{{8, 8}, 8}, {{16, 16}, 8}, {{16, 16}, 16}, {{32, 32}, 16}, {{32, 32}, 32}, {{64, 64}, 32}},
   {{{24, 1}, 8}, {{48, 2}, 8}, {{48, 2}, 16}, {{96, 4}, 16}, {{96, 4}, 32}},
   {{{8, 8, 8}, 8}, {{16, 16, 16}, 8}, {{16, 16, 16}, 16}}},
  TestID -> "nd-automatic-levels"]

(* ---------- presentation ---------- *)

(* numbers of modulus 10 or more and gaps of 10 % or more keep their leading digits *)
VerificationTest[
  {CertiField`Private`lowerString[10.679217, 6], CertiField`Private`upperString[10.952917, 6], CertiField`Private`lowerString[-12.5, 6],
   CertiField`Private`upperString[1234567.25, 6], CertiField`Private`lowerString[0.9375, 6], CertiField`Private`upperString[-0.015625, 6],
   CertiField`Private`percentString[124043/100000], CertiField`Private`percentString[124/1000], CertiField`Private`percentString[57/1000]},
  {"10.6792", "10.9530", "-12.5000", "1234570", "0.937500", "-0.0156250", "125%", "12.4%", "5.70%"},
  TestID -> "decimal-strings-keep-leading-digits"]

VerificationTest[
  StringQ[ndSq["Summary"]] && StringContainsQ[ndSq["Summary"], "Certified"] && StringContainsQ[ToString[ndSq["Theorem"]], "multilinear interpolant"] &&
    StringContainsQ[ToString[ndWall["Theorem"]], "multilinear interpolant"] && StringContainsQ[ToString[ndSq["NotClaimed"]], "minimiser"],
  True, TestID -> "nd-theorem-states-the-certified-datum"]
VerificationTest[Head[ToBoxes[ndWall]], InterpretationBox | TagBox | TemplateBox, SameTest -> MatchQ, TestID -> "nd-summary-box-renders"]

(* ---------- failure modes ---------- *)

VerificationTest[
  FailureQ[CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, 1 - Cos[u]], "Domain" -> {{0, 1}, {0, 1}}|>]]],
  True, TestID -> "nd-non-polynomial-potential-refused"]
VerificationTest[
  FailureQ[CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, u^2], "Domain" -> {{0, 1}}|>]]],
  True, TestID -> "nd-one-dimensional-domain-refused"]
VerificationTest[
  FailureQ[CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, -u], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 4|>]]],
  True, TestID -> "nd-potential-unbounded-below-is-inconclusive"]
VerificationTest[
  {FailureQ[CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, u^2], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 4|>], Method -> "NoSuchMethod"]],
   FailureQ[CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, (u^2 - 1)^2], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 4|>],
     Method -> {"Calibration", "Levels" -> {{{3, 3}, 4}, {{4, 4}, 4}}}]],
   FailureQ[CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, (u^2 - 1)^2], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 4,
     "Field" -> <|"TRange" -> {-1, 1}, "K" -> 2, "Affine" -> False, "S0" -> {}, "C" -> {}|>|>]]]},
  {True, True, True}, TestID -> "nd-bad-method-levels-and-field-refused"]
