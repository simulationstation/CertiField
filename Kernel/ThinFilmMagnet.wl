(* ::Package:: *)

(* ::Title:: *)
(*CertiField: thin-film magnetisation (unit-vector fields of a prescribed degree)*)

(* ::Text:: *)
(*E[m] = Integral over a rectangle of  A |grad m|^2 + kappa m.(curl m) + K (1 - m3^2) + h (1 - m3),  m : Omega -> S^2,*)
(*m = e3 on the boundary, degree Q.  docs/MAGNETISATION.md has the statements and proofs.*)
(**)
(*Lower bound (Proposition 5): a dual certificate (phi, u, Gamma) of piecewise polynomials in t = m3 with rational*)
(*coefficients, checked in exact rational arithmetic (Bernstein coefficients with bisection);*)
(*    E[m] >= 2 pi Q (Integral_{-1}^{1} u dt - 4A) - Gamma |Omega|.*)
(*Upper bound (Section 5): the normalised bilinear interpolant m_h = M_h/|M_h| of nodal vectors on a graded grid; its*)
(*energy and its degree are enclosed with the outward-rounded interval helpers of ScalarFieldND.wl.*)
(*Heuristic, never trusted: the radial profile of the competitor (Newton) and the choice of the dual certificate.*)

Begin["CertiField`Private`"];

(* ================================================================================================ *)
(* the model                                                                                        *)
(* ================================================================================================ *)

tfmExactQ[x_] := MatchQ[x, _Integer | _Rational];
tfmDyadicQ[x_] := tfmExactQ[x] && exactNumber[N[x]] == x;

tfmSpec[params_Association] := Module[{p = exactify[params], A, k, K, h, Q, dom, cells, beta},
  A = Lookup[p, "A", 1/2]; k = Lookup[p, "DMI", 0]; K = Lookup[p, "K", 0]; h = Lookup[p, "Field", 0]; Q = Lookup[p, "Degree", -1];
  dom = Lookup[p, "Domain", Missing[]]; cells = Lookup[p, "Cells", 200]; beta = Lookup[p, "Grading", Automatic];
  If[! (tfmExactQ[A] && A > 0 && tfmExactQ[k] && k >= 0 && tfmExactQ[K] && tfmExactQ[h] && h >= 0),
    Throw[failure["BadModel", "\"A\" > 0, \"DMI\" >= 0, \"K\" and \"Field\" >= 0 must be rational numbers."], "CertiField"]];
  If[! IntegerQ[Q], Throw[failure["BadModel", "\"Degree\" must be an integer."], "CertiField"]];
  If[! (MatchQ[dom, {{_, _}, {_, _}}] && And @@ (tfmDyadicQ /@ Flatten[dom]) && dom[[1, 1]] < dom[[1, 2]] && dom[[2, 1]] < dom[[2, 2]]),
    Throw[failure["BadModel", "\"Domain\" must be a rectangle {{a1, b1}, {a2, b2}} whose corners are dyadic rationals (exact machine numbers), e.g. {{-4, 4}, {-4, 4}}."], "CertiField"]];
  If[! (IntegerQ[cells] && cells >= 4 && EvenQ[cells]), Throw[failure["BadModel", "\"Cells\" must be an even integer >= 4 (cells per direction)."], "CertiField"]];
  If[! (beta === Automatic || (NumericQ[beta] && beta >= 0)), Throw[failure["BadModel", "\"Grading\" must be Automatic or a number >= 0."], "CertiField"]];
  <|"A" -> A, "DMI" -> k, "K" -> K, "Field" -> h, "Degree" -> Q, "Domain" -> dom, "Cells" -> cells, "Grading" -> beta,
    "Area" -> (dom[[1, 2]] - dom[[1, 1]]) (dom[[2, 2]] - dom[[2, 1]])|>];

tfmModelOf[spec_] := <|"A" -> spec["A"], "DMI" -> spec["DMI"], "K" -> spec["K"], "Field" -> spec["Field"]|>;

(* ================================================================================================ *)
(* the lower bound: exact check of a dual certificate                                               *)
(* ================================================================================================ *)

tfmPoly[c_List] := c . tfmT^Range[0, Length[c] - 1];
tfmV[mod_] := mod["K"] (1 - tfmT^2) + mod["Field"] (1 - tfmT);

(* Bernstein coefficients of the polynomial p (in tfmT) on [a, b], exact *)
tfmBernstein[p_, a_, b_] := Module[{q = Expand[p /. tfmT -> a + (b - a) tfmS], n, c},
  n = Max[Exponent[q, tfmS], 0];
  c = CoefficientList[q, tfmS, n + 1];
  Table[Sum[Binomial[k, i]/Binomial[n, i] c[[i + 1]], {i, 0, k}], {k, 0, n}]];

tfmSplit[bz_List] := Module[{cur = bz, left = {First[bz]}, right = {Last[bz]}},
  While[Length[cur] > 1, cur = (Most[cur] + Rest[cur])/2; AppendTo[left, First[cur]]; AppendTo[right, Last[cur]]];
  {left, Reverse[right]}];

(* True only if p >= 0 on [a, b] (exact; False means not proved) *)
tfmNonneg[p_, a_, b_] := Module[{stack = {{tfmBernstein[p, a, b], 0}}, bz, d, lr, ok = True},
  While[ok && stack =!= {},
    {bz, d} = Last[stack]; stack = Most[stack];
    If[Min[bz] < 0,
      If[First[bz] < 0 || Last[bz] < 0 || d >= 60, ok = False,
        lr = tfmSplit[bz]; stack = Join[stack, {{lr[[1]], d + 1}, {lr[[2]], d + 1}}]]]];
  ok];

tfmCellPolys[dual_, i_] := Module[{mod = dual["Model"], k, phi, u, a, b, P, R},
  k = mod["DMI"]; phi = tfmPoly[dual["Phi"][[i]]]; u = Expand[tfmPoly[dual["U"][[i]]]];
  a = (1 + tfmT) (phi + k); b = (1 - tfmT) (phi - k);
  P = Expand[((1 - tfmT) D[a, tfmT])^2/4]; R = Expand[((1 + tfmT) D[b, tfmT])^2/4];
  <|"P" -> P, "R" -> R, "VG" -> Expand[tfmV[mod] + dual["Gamma"]], "u" -> u, "w" -> Expand[4 mod["A"] - u]|>];

(* the conditions of Proposition 5, in exact arithmetic *)
tfmCheckDual[dual_Association] := Catch[Module[{mesh = dual["Mesh"], A = dual["Model", "A"], n, c, t0, t1, p1, p2, t},
  n = Length[mesh] - 1;
  If[! (VectorQ[mesh, tfmExactQ] && n >= 1 && mesh[[1]] == -1 && mesh[[-1]] == 1 && And @@ Thread[Most[mesh] < Rest[mesh]] &&
      Length[dual["Phi"]] == n && Length[dual["U"]] == n && And @@ (VectorQ[#, tfmExactQ] & /@ Join[dual["Phi"], dual["U"]]) && tfmExactQ[dual["Gamma"]]),
    Throw[False, "tfm"]];
  Do[p1 = tfmPoly[dual["Phi"][[i]]]; p2 = tfmPoly[dual["Phi"][[i + 1]]]; t = mesh[[i + 1]];
    If[(p1 /. tfmT -> t) != (p2 /. tfmT -> t) || (D[p1, tfmT] /. tfmT -> t) != (D[p2, tfmT] /. tfmT -> t), Throw[False, "tfm"]], {i, n - 1}];
  Do[c = tfmCellPolys[dual, i]; t0 = mesh[[i]]; t1 = mesh[[i + 1]];
    If[! Which[
        c["u"] === 0, c["P"] === 0 && tfmNonneg[4 A c["VG"] - c["R"], t0, t1],
        c["w"] === 0, c["R"] === 0 && tfmNonneg[4 A c["VG"] - c["P"], t0, t1],
        True, tfmNonneg[c["u"], t0, t1] && tfmNonneg[c["w"], t0, t1] &&
          tfmNonneg[-(c["P"] c["w"] + c["R"] c["u"] - c["VG"] c["u"] c["w"]), t0, t1]],
      Throw[False, "tfm"]], {i, n}];
  True], "tfm"];

tfmIntegralU[dual_] := Module[{m = dual["Mesh"]},
  Sum[With[{c = dual["U"][[i]]}, Sum[c[[j]] (m[[i + 1]]^j - m[[i]]^j)/j, {j, Length[c]}]], {i, Length[m] - 1}]];

(* the bound c pi + d, and an exact rational below it *)
tfmBoundCoefficients[dual_, Q_, area_] := {2 Q (tfmIntegralU[dual] - 4 dual["Model", "A"]), -dual["Gamma"] area};
tfmLowerExact[dual_, Q_, area_] := With[{cd = tfmBoundCoefficients[dual, Q, area]}, cd[[1]] If[cd[[1]] >= 0, $piLo, $piHi] + cd[[2]]];

tfmMinV[mod_] := Module[{K = mod["K"], h = mod["Field"], ts = {-1, 1}},
  If[K != 0 && -1 < -h/(2 K) < 1, AppendTo[ts, -h/(2 K)]];
  Min[(K (1 - #^2) + h (1 - #)) & /@ ts]];

(* closed-form certificates (Section 3 of the note): constant phi *)
tfmSimpleDuals[mod_, Q_] := Module[{A = mod["A"], k = mod["DMI"], K = mod["K"], h = mod["Field"], g0, v0, out = {}, mk},
  g0 = Max[0, -tfmMinV[mod]];
  mk = Function[{phi, u, gam}, <|"Model" -> mod, "Mesh" -> {-1, 1}, "Phi" -> {phi}, "U" -> {u}, "Gamma" -> gam, "Source" -> "closed form"|>];
  If[k == 0,
    AppendTo[out, mk[{0}, {If[Q > 0, 4 A, 0]}, g0]],
    (* phi = kappa: P = kappa^2 (1 - t)^2, R = 0, V = (1 - t)(h + K (1 + t)) *)
    If[h == -2 K && h > 0 && 2 k^2/h <= 4 A, AppendTo[out, mk[{k}, {2 k^2/h}, 0]]];
    v0 = Min[h, h + 2 K];
    If[v0 > 0 && 2 k^2/v0 <= 4 A, AppendTo[out, mk[{k}, {k^2/v0, -k^2/v0}, 0]]];
    AppendTo[out, mk[{k}, {4 A}, k^2/A + g0]]];
  out];

(* certificates computed by the Python reference and shipped with the paclet (Data/Magnetisation) *)
tfmImportDual[file_String] := Module[{d = Quiet@Import[file, "RawJSON"], mod},
  If[! AssociationQ[d] || d["format"] =!= "CertiField-magnet-dual-1", Throw[failure["MissingData", "Not a dual certificate: " <> file], "CertiField"]];
  mod = <|"A" -> ndRational[d["model", "A"]], "DMI" -> ndRational[d["model", "kappa"]], "K" -> ndRational[d["model", "K"]], "Field" -> ndRational[d["model", "h"]]|>;
  <|"Model" -> mod, "Mesh" -> ndRational /@ d["mesh"], "Phi" -> Map[ndRational, d["phi"], {2}], "U" -> Map[ndRational, d["u"], {2}],
    "Gamma" -> ndRational[d["Gamma"]], "Source" -> FileNameTake[file]|>];

tfmShippedDuals[mod_] := Module[{files = FileNames["*.json", FileNameJoin[{$CertiFieldDirectory, "Data", "Magnetisation"}]]},
  Select[tfmImportDual /@ files, #["Model"] === mod &]];

(* the best certificate that passes the exact check *)
tfmChooseDual[spec_, given_] := Module[{mod = tfmModelOf[spec], cands, valid, vals},
  cands = Which[
    given === Automatic, Join[tfmShippedDuals[mod], tfmSimpleDuals[mod, spec["Degree"]]],
    StringQ[given], {tfmImportDual[given]},
    AssociationQ[given], {Join[<|"Gamma" -> 0, "Source" -> "given"|>, given, <|"Model" -> mod|>]},
    True, Throw[failure["BadModel", "\"DualField\" must be Automatic, a file written by the Python reference, or an association with \"Mesh\", \"Phi\", \"U\" (and \"Gamma\")."], "CertiField"]];
  If[given =!= Automatic && cands[[1]]["Model"] =!= mod, Throw[failure["BadModel", "The dual certificate belongs to a different model."], "CertiField"]];
  valid = Select[cands, tfmCheckDual];
  If[valid === {}, Throw[failure["Inconclusive", "No dual certificate passed the exact check."], "CertiField"]];
  vals = N[tfmLowerExact[#, spec["Degree"], spec["Area"]]] & /@ valid;
  valid[[First[Ordering[vals, -1]]]]];

(* ================================================================================================ *)
(* the upper bound: energy and degree of the normalised bilinear interpolant (rigorous)             *)
(* ================================================================================================ *)

tfmDot[u_List, v_List] := Fold[ndAdd, MapThread[ndMul, {u, v}]];
tfmCross[{a1_, a2_, a3_}, {b1_, b2_, b3_}] := {ndSub[ndMul[a2, b3], ndMul[a3, b2]], ndSub[ndMul[a3, b1], ndMul[a1, b3]], ndSub[ndMul[a1, b2], ndMul[a2, b1]]};
tfmHalfSum[u_List, v_List] := MapThread[ivK[ndAdd[#1, #2], 1/2] &, {u, v}];
tfmInv[x_] := ndDiv[{1., 1.}, x];

(* corner vectors of every cell as point intervals; the boundary nodes must be exactly (0, 0, 1) *)
tfmCells[x1_, x2_, Mn_] := Module[{n1 = Length[x1] - 1, n2 = Length[x2] - 1, edge, h1, h2, cor},
  If[! (Developer`PackedArrayQ[Mn, Real] && Dimensions[Mn] == {n1 + 1, n2 + 1, 3} && Developer`PackedArrayQ[x1, Real] && Developer`PackedArrayQ[x2, Real] &&
      Min[Differences[x1]] > 0 && Min[Differences[x2]] > 0),
    Throw[failure["BadCandidate", "The grid and the nodal vectors do not fit."], "CertiField"]];
  edge = Join[Mn[[1]], Mn[[-1]], Mn[[All, 1]], Mn[[All, -1]]];
  If[! And @@ (# === {0., 0., 1.} & /@ edge), Throw[failure["BadCandidate", "The boundary nodes must be exactly (0, 0, 1)."], "CertiField"]];
  cor[i_, j_] := Table[ivPt[Mn[[i ;; i + n1 - 1, j ;; j + n2 - 1, c]]], {c, 3}];
  h1 = ndChk[ivSub[ivPt[Rest[x1]], ivPt[Most[x1]]]]; h2 = ndChk[ivSub[ivPt[Rest[x2]], ivPt[Most[x2]]]];
  <|"00" -> cor[1, 1], "10" -> cor[2, 1], "01" -> cor[1, 2], "11" -> cor[2, 2],
    "h1" -> {Transpose[ConstantArray[h1[[1]], n2]], Transpose[ConstantArray[h1[[2]], n2]]},
    "h2" -> {ConstantArray[h2[[1]], n1], ConstantArray[h2[[2]], n1]}|>];

(* Bernstein coefficients c[[k + 1, l + 1]] of sum over comps of M_c^2 on every cell *)
tfmDCoefficients[B_, comps_List] := Module[{d},
  d[p_, q_] := tfmDot[B[p][[comps]], B[q][[comps]]];
  {{d["00", "00"], d["00", "01"], d["01", "01"]},
   {d["00", "10"], ivK[ndAdd[d["00", "11"], d["10", "01"]], 1/2], d["01", "11"]},
   {d["10", "10"], d["10", "11"], d["11", "11"]}}];

tfmRange[ivs_List] := {Fold[vMin, #[[1]] & /@ ivs], Fold[vMax, #[[2]] & /@ ivs]};

$tfmG2 = Table[Binomial[2, i] Binomial[2, j]/(5 Binomial[4, i + j]), {i, 0, 2}, {j, 0, 2}];
$tfmT2 = Table[Binomial[2, i] Binomial[2, j] Binomial[2, l]/(7 Binomial[6, i + j + l]), {i, 0, 2}, {j, 0, 2}, {l, 0, 2}];
$tfmG12 = Table[Binomial[1, i] Binomial[2, k]/(4 Binomial[3, i + k]), {i, 0, 1}, {k, 0, 2}];

(* RIGOROUS: an exact rational >= E[m_h] *)
tfmUpperBound[mod_, x1_, x2_, Mn_] := ndQuiet@Module[
  {B = tfmCells[x1, x2, Mn], c, alo, bhi, a, b, Cl, Ck, v, w, a2b2, al2, be2, quartic, X, Y, ex, total, d0, gap, linear, dm, area, q, q0, q1, an,
   m3, m0, m1, neg, gz, ze, h1, h2},
  h1 = B["h1"]; h2 = B["h2"];
  c = tfmDCoefficients[B, {1, 2, 3}];
  {alo, bhi} = tfmRange[Flatten[c, 1]];
  If[Min[alo] <= 0, Throw[failure["Inconclusive", "The interpolant of the candidate may vanish: |M_h|^2 > 0 could not be shown on every cell."], "CertiField"]];
  a = ivPt[alo]; b = ivPt[bhi];
  Cl = Table[ivK[ivSum[c[[All, l]]], 1/3], {l, 3}];
  Ck = Table[ivK[ivSum[c[[k]]], 1/3], {k, 3}];
  (* exchange *)
  v = {tfmCross[B["00"], B["10"]], tfmHalfSum[tfmCross[B["00"], B["11"]], tfmCross[B["01"], B["10"]]], tfmCross[B["01"], B["11"]]};
  w = {tfmCross[B["00"], B["01"]], tfmHalfSum[tfmCross[B["00"], B["11"]], tfmCross[B["10"], B["01"]]], tfmCross[B["10"], B["11"]]};
  a2b2 = ndMul[ndMul[a, a], ndMul[b, b]];
  al2 = ndDiv[ivSum[{ndMul[a, a], ndMul[a, b], ndMul[b, b]}], a2b2];
  be2 = ndDiv[ndAdd[a, b], a2b2];
  quartic[vv_, Cd_] := Module[{s0 = {}, s1 = {}, dij},
    Do[dij = ivK[tfmDot[vv[[i]], vv[[j]]], If[i == j, 1, 2]];
      AppendTo[s0, ivK[dij, $tfmG2[[i, j]]]];
      Do[AppendTo[s1, ivK[ndMul[dij, Cd[[l]]], $tfmT2[[i, j, l]]]], {l, 3}], {i, 3}, {j, i, 3}];
    {ivSum[s0], ivSum[s1]}];
  X = quartic[v, Cl]; Y = quartic[w, Ck];
  ex = ivK[ndAdd[ndMul[ndSub[ndMul[al2, X[[1]]], ndMul[be2, X[[2]]]], ndDiv[h2, h1]], ndMul[ndSub[ndMul[al2, Y[[1]]], ndMul[be2, Y[[2]]]], ndDiv[h1, h2]]], mod["A"]];
  total = ex;
  area = ndMul[h1, h2];
  (* Dzyaloshinskii-Moriya: -kappa [h2 Int e1.v/D + h1 Int e2.w/D], 1/D = 2/d0 - D/d0^2 + [0, gap] *)
  If[mod["DMI"] != 0,
    d0 = ivPt[(alo + bhi)/2];
    gap = {0. alo, vMax[ndDiv[ivSqr[ndSub[a, d0]], ndMul[a, ndMul[d0, d0]]][[2]], ndDiv[ivSqr[ndSub[b, d0]], ndMul[b, ndMul[d0, d0]]][[2]]]};
    linear[vv_, comp_, Cd_] := Module[{n = vv[[All, comp]], n0, n1, nmax},
      n0 = ivK[ivSum[n], 1/3];
      n1 = ivSum[Flatten[Table[ivK[ndMul[n[[i]], Cd[[l]]], $tfmG2[[i, l]]], {i, 3}, {l, 3}], 1]];
      nmax = Fold[vMax, (vMax[Abs[#[[1]]], Abs[#[[2]]]] &) /@ n];
      ndAdd[ndSub[ndDiv[ivK[n0, 2], d0], ndDiv[n1, ndMul[d0, d0]]], ndMul[{-nmax, nmax}, gap]]];
    dm = ivK[ndAdd[ndMul[linear[v, 1, Cl], h2], ndMul[linear[w, 2, Ck], h1]], -mod["DMI"]];
    total = ndAdd[total, dm]];
  (* anisotropy: K (M1^2 + M2^2)/D *)
  If[mod["K"] != 0,
    q = tfmDCoefficients[B, {1, 2}];
    q0 = ivK[ivSum[Flatten[q, 1]], 1/9];
    q1 = ivSum[Flatten[Table[ivK[ndMul[q[[k, l]], c[[k2, l2]]], $tfmG2[[k, k2]] $tfmG2[[l, l2]]], {k, 3}, {l, 3}, {k2, 3}, {l2, 3}], 3]];
    an = If[mod["K"] > 0,
      ndDiv[ndSub[ndMul[q0, ndAdd[a, b]], q1], ndMul[a, b]],
      With[{dd = ivPt[(alo + bhi)/2]}, ndSub[ndDiv[ivK[q0, 2], dd], ndDiv[q1, ndMul[dd, dd]]]]];
    total = ndAdd[total, ivK[ndMul[an, area], mod["K"]]]];
  (* Zeeman: h (1 - M3/|M|),  1/|M| = (3 - D)/2 + r,  0 <= r <= (1 - D)^2/(2D) *)
  If[mod["Field"] != 0,
    m3 = #[[3]] & /@ {B["00"], B["10"], B["01"], B["11"]};
    m0 = ivK[ivSum[m3], 1/4];
    m1 = ivSum[Flatten[Table[ivK[ndMul[{B["00"], B["01"], B["10"], B["11"]}[[1 + 2 i + j, 3]], c[[k, l]]], $tfmG12[[i + 1, k]] $tfmG12[[j + 1, l]]], {i, 0, 1}, {j, 0, 1}, {k, 3}, {l, 3}], 3]];
    neg = Ramp[-Fold[vMin, #[[1]] & /@ m3]];
    gz = vMax[ndDiv[ivSqr[ndSub[{1., 1.}, a]], ivK[a, 2]][[2]], ndDiv[ivSqr[ndSub[{1., 1.}, b]], ivK[b, 2]][[2]]];
    ze = ivSum[{{1., 1.}, ivK[m0, -3/2], ivK[m1, 1/2], ndMul[ivPt[neg], ivPt[gz]]}];
    total = ndAdd[total, ivK[ndMul[ze, area], mod["Field"]]]];
  exactTotal[total[[2]]]];

(* RIGOROUS: {Q or None, {lo, hi}} with lo <= Integral m.(d1 m x d2 m) <= hi (exact rationals) *)
tfmDegree[x1_, x2_, Mn_] := ndQuiet@Module[{B = tfmCells[x1, x2, Mn], c, alo, bhi, a, b, det, corners, d, wlo, whi, cell, lo, hi, cands},
  c = tfmDCoefficients[B, {1, 2, 3}];
  {alo, bhi} = tfmRange[Flatten[c, 1]];
  If[Min[alo] <= 0, Return[{None, {-Infinity, Infinity}}, Module]];
  a = ivPt[alo]; b = ivPt[bhi];
  det[p_, q_, r_] := tfmDot[B[p], tfmCross[B[q], B[r]]];
  corners = {det["00", "10", "01"], det["00", "10", "11"], det["00", "11", "01"], det["10", "11", "01"]};
  d = tfmRange[corners];
  wlo = ndDiv[{2., 2.}, ndMul[b, ndAdd[b, {1., 1.}]]][[1]];
  whi = ndDiv[ndAdd[a, {1., 1.}], ivK[ndMul[a, a], 2]][[2]];
  cell = ndMul[d, {wlo, whi}];
  lo = exactTotal[cell[[1]]]; hi = exactTotal[cell[[2]]];
  cands = Select[Range[Floor[Min[lo/12, lo/13]] - 2, Ceiling[Max[hi/12, hi/13]] + 2],
    ! (Max[4 $piLo #, 4 $piHi #] < lo || Min[4 $piLo #, 4 $piHi #] > hi) &];
  {If[Length[cands] == 1, First[cands], None], {lo, hi}}];

(* ================================================================================================ *)
(* the competitor (heuristic: any nodal values give a valid bound)                                  *)
(* ================================================================================================ *)

(* 2n cells on [a, b], symmetric about the midpoint, finer near it *)
tfmGradedAxis[a_, b_, n_, beta_] := Module[{tau = N[Range[0, n]/n], half, mid = N[(a + b)/2], hw = N[(b - a)/2], x},
  half = If[beta > 0, hw Sinh[beta tau]/Sinh[beta], hw tau];
  x = Join[mid - Reverse[Rest[half]], {mid}, mid + Rest[half]];
  x[[1]] = N[a]; x[[-1]] = N[b];
  Developer`ToPackedArray[x]];

(* the axisymmetric profile theta(r), theta(0) = pi, theta(L) = 0, minimising the radial energy (Newton) *)
tfmProfile[mod_, L_, n_ : 3000] := Module[
  {A = N[mod["A"]], k = N[mod["DMI"]], K = N[mod["K"]], h = N[mod["Field"]], stiff, r, rm, dr, w, parts, th, lam = 10.^-3, en, Ld, Lm, Ldd, Lmm,
   g, d, o, gi, step, tn, it = 0, accepted, mat, R0, wm},
  stiff = Max[h + 2 K, 10.^-3];
  r = Join[{0.}, Exp[Subdivide[Log[Min[10.^-4, 10.^-3 N[L]]/(1 + stiff)], Log[N[L]], n - 1]]];
  r[[-1]] = N[L];
  rm = (Most[r] + Rest[r])/2; dr = Differences[r]; w = 2. Pi rm dr; wm = Mean[w];
  parts[t_] := Module[{dth = Differences[t]/dr, tm = (Most[t] + Rest[t])/2, s, c},
    s = Sin[tm]; c = Cos[tm];
    {Total[w (A (dth^2 + s^2/rm^2) + k (dth + s c/rm) + h (1 - c) + K s^2)], w (2 A dth + k),
     w (2 A s c/rm^2 + k (c^2 - s^2)/rm + h s + 2 K s c), w 2 A, w (2 A (c^2 - s^2)/rm^2 - 4 k s c/rm + h c + 2 K (c^2 - s^2))}];
  R0 = (Abs[k] + 10.^-3)/(stiff Max[Log[stiff + 2.], 1.]);
  th = 2 ArcTan[r, R0] Exp[-Sqrt[stiff/(2 A)] r/2];
  th[[1]] = N[Pi]; th[[-1]] = 0.;
  While[it++ < 500,
    {en, Ld, Lm, Ldd, Lmm} = parts[th];
    g = Join[{0.}, Ld/dr + Lm/2] + Join[-Ld/dr + Lm/2, {0.}];
    d = Join[{0.}, Ldd/dr^2 + Lmm/4] + Join[Ldd/dr^2 + Lmm/4, {0.}];
    o = -Ldd/dr^2 + Lmm/4;
    gi = g[[2 ;; -2]];
    If[Max[Abs[gi]] < 10.^-12 (1 + Abs[en]), Break[]];
    accepted = False;
    While[lam < 10.^12 && ! accepted,
      mat = SparseArray[{Band[{1, 1}] -> d[[2 ;; -2]] + lam wm, Band[{1, 2}] -> o[[2 ;; -2]], Band[{2, 1}] -> o[[2 ;; -2]]}, {Length[gi], Length[gi]}];
      step = Quiet@LinearSolve[mat, -gi];
      If[VectorQ[step, NumericQ] && gi . step < 0,
        tn = th; tn[[2 ;; -2]] += step;
        If[parts[tn][[1]] <= en + 10.^-4 gi . step, th = tn; lam = Max[lam/5, 10.^-12]; accepted = True]];
      If[! accepted, lam *= 10]];
    If[! accepted, Break[]]];
  Interpolation[Transpose[{r, th}], InterpolationOrder -> 1]];

(* nodal vectors of the axisymmetric field (-sin(theta) y/r, sin(theta) x/r, cos(theta)) about the centre of the
   domain; e3 on the boundary and outside the inscribed circle *)
tfmRadialNodes[spec_, x1_, x2_, theta_] := Module[{dom = spec["Domain"], cx, cy, L, X, Y, r, th, s, inside, Mn},
  cx = N[Mean[dom[[1]]]]; cy = N[Mean[dom[[2]]]]; L = N[Min[(dom[[1, 2]] - dom[[1, 1]])/2, (dom[[2, 2]] - dom[[2, 1]])/2]];
  X = Transpose[ConstantArray[x1 - cx, Length[x2]]]; Y = ConstantArray[x2 - cy, Length[x1]];
  r = Sqrt[X^2 + Y^2];
  th = ArrayReshape[tfmEvalList[theta, Flatten[Clip[r, {0., L}]]], Dimensions[r]];
  s = Sin[th]/(r + UnitStep[-r]);
  Mn = Transpose[{-s Y, s X, Cos[th]}, {3, 1, 2}];
  inside = 1 - UnitStep[r - L];
  inside[[1]] = inside[[-1]] = 0 inside[[1]]; inside[[All, 1]] = 0; inside[[All, -1]] = 0;
  Mn = Mn inside + Transpose[{0. inside, 0. inside, 1. - inside}, {3, 1, 2}];
  Developer`ToPackedArray[N[Mn]]];

tfmEvalList[f_, list_] := Module[{v = Quiet[f[list]]}, If[VectorQ[v, Internal`RealValuedNumericQ] && Length[v] == Length[list], N[v], N[f /@ list]]];

tfmAutoGrading[spec_, theta_, n_] := Module[{L, rr, core, f},
  L = N[Min[(spec["Domain"][[1, 2]] - spec["Domain"][[1, 1]])/2, (spec["Domain"][[2, 2]] - spec["Domain"][[2, 1]])/2]];
  rr = Subdivide[0., L, 4000];
  core = Max[rr[[Replace[Flatten[Position[UnitStep[Pi/2 - tfmEvalList[theta, rr]], 1, 1, 1]], {} -> {2}][[1]]]], L/1000];
  (* smallest cell about core/15 *)
  f[beta_?NumericQ] := L beta/(n Sinh[beta]) - core/15;
  If[f[0.01] <= 0, 0., Quiet@Check[beta /. FindRoot[f[beta], {beta, 2., 0.01, 30.}], 4.]]];

tfmCompetitor[spec_, cand_] := Module[{n = spec["Cells"]/2, theta, beta, x1, x2, Mn, f},
  theta = Which[
    cand === Automatic,
      If[spec["Degree"] =!= -1, Throw[failure["BadCandidate", "The automatic competitor has degree -1; for other degrees give a field Function[{x, y}, {m1, m2, m3}]."], "CertiField"]];
      tfmProfile[tfmModelOf[spec], Min[(spec["Domain"][[1, 2]] - spec["Domain"][[1, 1]])/2, (spec["Domain"][[2, 2]] - spec["Domain"][[2, 1]])/2]],
    True, None];
  f = If[theta === None && ! MatchQ[cand, _Function | _InterpolatingFunction], Throw[failure["BadCandidate", "The candidate must be Automatic, a profile Function[r, theta] or a field Function[{x, y}, {m1, m2, m3}]."], "CertiField"], cand];
  If[theta === None && NumericQ[Quiet[f[0.5]]], theta = f];
  beta = Replace[spec["Grading"], Automatic :> If[theta =!= None, tfmAutoGrading[spec, theta, n], 3.]];
  x1 = tfmGradedAxis[spec["Domain"][[1, 1]], spec["Domain"][[1, 2]], n, beta];
  x2 = tfmGradedAxis[spec["Domain"][[2, 1]], spec["Domain"][[2, 2]], n, beta];
  Mn = If[theta =!= None, tfmRadialNodes[spec, x1, x2, theta],
    Module[{vals = Outer[N[f[#1, #2]] &, x1, x2]},
      If[! ArrayQ[vals, 3, Internal`RealValuedNumericQ] || Dimensions[vals] =!= {Length[x1], Length[x2], 3},
        Throw[failure["BadCandidate", "The field must return three real numbers at every node."], "CertiField"]];
      vals[[1]] = vals[[-1]] = ConstantArray[{0., 0., 1.}, Length[x2]]; vals[[All, 1]] = vals[[All, -1]] = ConstantArray[{0., 0., 1.}, Length[x1]];
      Developer`ToPackedArray[N[vals]]]];
  {x1, x2, Mn, beta}];

(* ================================================================================================ *)
(* driver, description, verification                                                                *)
(* ================================================================================================ *)

certifyThinFilmMagnet[params_Association, cand_, opts_List] := Module[{spec = tfmSpec[params], t0 = AbsoluteTime[], dual, lb, x1, x2, Mn, beta, deg, ub},
  dual = tfmChooseDual[spec, Lookup[params, "DualField", Automatic]];
  lb = tfmLowerExact[dual, spec["Degree"], spec["Area"]];
  {x1, x2, Mn, beta} = tfmCompetitor[spec, cand];
  deg = tfmDegree[x1, x2, Mn];
  If[deg[[1]] =!= spec["Degree"],
    Throw[failure["Inconclusive", "The degree of the competitor could not be certified to be " <> ToString[spec["Degree"]] <> " (enclosure of 4\[Pi]Q: " <> ToString[N[deg[[2]]]] <> ")."], "CertiField"]];
  ub = tfmUpperBound[tfmModelOf[spec], x1, x2, Mn];
  <|"LowerBound" -> lb, "UpperBound" -> ub, "Seconds" -> AbsoluteTime[] - t0, "Dual" -> dual,
    "BoundCoefficients" -> tfmBoundCoefficients[dual, spec["Degree"], spec["Area"]], "Grid" -> {x1, x2}, "Grading" -> beta, "Nodal" -> Mn,
    "DegreeEnclosure" -> deg[[2]], "Spec" -> spec|>];

tfmDescription[params_, res_] := Module[{sp = res["Spec"], d = Max[displayDigits[res], 8], cd = res["BoundCoefficients"]},
  <|"Units" -> "energy units of the model",
    "Theorem" -> Row[{"Every m in H^1 of the rectangle ", exactString[sp["Domain"]], " with values in S^2, m = (0, 0, 1) on the boundary and degree ", sp["Degree"],
      " has E[m] \[GreaterEqual] ", lowerString[res["LowerBound"], d], " (exactly: E[m] \[GreaterEqual] ", exactString[cd[[1]]], " \[Pi]",
      If[cd[[2]] == 0, "", Row[{" + (", exactString[cd[[2]]], ")"}]], "); the competitor, of certified degree ", sp["Degree"], ", has E \[LessEqual] ", upperString[res["UpperBound"], d], "."}],
    "Energy" -> Row[{"E[m] = \[Integral] ", exactString[sp["A"]], " |\[Del]m|^2 + ", exactString[sp["DMI"]], " m\[CenterDot](\[Del]\[Times]m) + ", exactString[sp["K"]],
      " (1 - m3^2) + ", exactString[sp["Field"]], " (1 - m3) dx"}],
    "AdmissibleClass" -> Row[{"all m in H^1 of the rectangle ", exactString[sp["Domain"]], " with |m| = 1, m = (0, 0, 1) on the boundary (trace) and degree ",
      sp["Degree"], " (the degree is the integer (1/4\[Pi]) \[Integral] m\[CenterDot](\[PartialD]1 m \[Times] \[PartialD]2 m))"}],
    "LowerBoundMethod" -> Row[{"null Lagrangians of degree and divergence type with a dual field (\[Phi], u, \[CapitalGamma]) of piecewise polynomials in m3 on ",
      Length[res["Dual", "Mesh"]] - 1, " intervals (", res["Dual", "Source"], "), checked in exact rational arithmetic"}],
    "UpperBoundMethod" -> Row[{"normalised bilinear interpolant on a graded ", Length[res["Grid"][[1]]] - 1, "\[Times]", Length[res["Grid"][[2]]] - 1,
      " grid; energy and degree enclosed with outward-rounded arithmetic"}],
    "TrustedBase" -> "IEEE-754 machine arithmetic (+, \[Minus], \[Times], \[Divide]), Wolfram Language exact rational arithmetic, and docs/MAGNETISATION.md (Lemmas 1\[Dash]3, Theorem 4, Proposition 5; proofs checked by one author only)",
    "NotClaimed" -> "that the competitor is a minimiser, or that the axisymmetric skyrmion minimises the energy in its degree class; the stray field is not included",
    "References" -> "docs/MAGNETISATION.md; Belavin\[Dash]Polyakov (1975); D\[ODoubleDot]ring\[Dash]Melcher (2017); Li\[Dash]Melcher, arXiv:1711.07717."|>];

verifyThinFilmMagnet[a_Association] := Catch[Module[{c = a["Certificate"], spec, dual, lb, ub, deg, x1, x2},
  spec = tfmSpec[a["Parameters"]];
  dual = c["Dual"];
  If[dual["Model"] =!= tfmModelOf[spec] || ! tfmCheckDual[dual], Throw[False, "tfmv"]];
  lb = tfmLowerExact[dual, spec["Degree"], spec["Area"]];
  {x1, x2} = c["Grid"];
  If[{x1[[1]], x1[[-1]], x2[[1]], x2[[-1]]} != N[Flatten[spec["Domain"]]], Throw[False, "tfmv"]];
  deg = tfmDegree[x1, x2, c["Nodal"]];
  If[deg[[1]] =!= spec["Degree"], Throw[False, "tfmv"]];
  ub = tfmUpperBound[tfmModelOf[spec], x1, x2, c["Nodal"]];
  lb >= exactNumber[a["LowerBound"]] && ub <= exactNumber[a["UpperBound"]]], "tfmv"];

End[];
