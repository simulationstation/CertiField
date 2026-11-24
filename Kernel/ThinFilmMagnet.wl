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

tfmSpec[params_Association] := Module[{p = exactify[params], A, k, K, h, Q, dom, cells, beta, t},
  A = Lookup[p, "A", 1/2]; k = Lookup[p, "DMI", 0]; K = Lookup[p, "K", 0]; h = Lookup[p, "Field", 0]; Q = Lookup[p, "Degree", -1];
  dom = Lookup[p, "Domain", Missing[]]; cells = Lookup[p, "Cells", 200]; beta = Lookup[p, "Grading", Automatic]; t = Lookup[p, "Thickness", 0];
  If[! (tfmExactQ[t] && t >= 0), Throw[failure["BadModel", "\"Thickness\" must be a rational number >= 0 (0: no stray field)."], "CertiField"]];
  If[! (tfmExactQ[A] && A > 0 && tfmExactQ[k] && tfmExactQ[K] && tfmExactQ[h] && h >= 0),
    Throw[failure["BadModel", "\"A\" > 0, \"DMI\", \"K\" and \"Field\" >= 0 must be rational numbers."], "CertiField"]];
  If[! IntegerQ[Q], Throw[failure["BadModel", "\"Degree\" must be an integer."], "CertiField"]];
  If[! (MatchQ[dom, {{_, _}, {_, _}}] && And @@ (tfmDyadicQ /@ Flatten[dom]) && dom[[1, 1]] < dom[[1, 2]] && dom[[2, 1]] < dom[[2, 2]]),
    Throw[failure["BadModel", "\"Domain\" must be a rectangle {{a1, b1}, {a2, b2}} whose corners are dyadic rationals (exact machine numbers), e.g. {{-4, 4}, {-4, 4}}."], "CertiField"]];
  If[! (IntegerQ[cells] && cells >= 4 && EvenQ[cells]), Throw[failure["BadModel", "\"Cells\" must be an even integer >= 4 (cells per direction)."], "CertiField"]];
  If[! (beta === Automatic || (NumericQ[beta] && beta >= 0)), Throw[failure["BadModel", "\"Grading\" must be Automatic or a number >= 0."], "CertiField"]];
  <|"A" -> A, "DMI" -> k, "K" -> K, "Field" -> h, "Degree" -> Q, "Domain" -> dom, "Cells" -> cells, "Grading" -> beta, "Thickness" -> t,
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
  Do[If[! tfmCellOK[tfmCellPolys[dual, i], A, mesh[[i]], mesh[[i + 1]]], Throw[False, "tfm"]], {i, n}];
  True], "tfm"];

tfmCellOK[c_, A_, t0_, t1_] := Which[
  c["u"] === 0, c["P"] === 0 && tfmNonneg[4 A c["VG"] - c["R"], t0, t1],
  c["w"] === 0, c["R"] === 0 && tfmNonneg[4 A c["VG"] - c["P"], t0, t1],
  True, tfmNonneg[c["u"], t0, t1] && tfmNonneg[c["w"], t0, t1] && tfmNonneg[-(c["P"] c["w"] + c["R"] c["u"] - c["VG"] c["u"] c["w"]), t0, t1]];

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
    AppendTo[out, mk[{k}, {4 A}, 0]];
    AppendTo[out, mk[{k}, {4 A}, k^2/A + g0]]];
  out];

(* certificates computed by the Python reference and shipped with the paclet (Data/Magnetisation) *)
tfmImportDual[file_String] := Module[{d = Quiet@Import[file, "RawJSON"], mod},
  If[! AssociationQ[d] || d["format"] =!= "CertiField-magnet-dual-1", Throw[failure["MissingData", "Not a dual certificate: " <> file], "CertiField"]];
  mod = <|"A" -> ndRational[d["model", "A"]], "DMI" -> ndRational[d["model", "kappa"]], "K" -> ndRational[d["model", "K"]], "Field" -> ndRational[d["model", "h"]]|>;
  <|"Model" -> mod, "Mesh" -> ndRational /@ d["mesh"], "Phi" -> Map[ndRational, d["phi"], {2}], "U" -> Map[ndRational, d["u"], {2}],
    "Gamma" -> ndRational[d["Gamma"]], "Source" -> FileNameTake[file]|>];

(* the cache: certificates for this model, or for the model with -kappa with phi -> -phi (P and R are unchanged) *)
tfmShippedDuals[mod_] := Module[{files = FileNames["*.json", FileNameJoin[{$CertiFieldDirectory, "Data", "Magnetisation"}]], all},
  all = tfmImportDual /@ files;
  Join[Select[all, tfmSameModel[#["Model"], mod, 1] &], If[mod["DMI"] != 0, tfmFlipDual[#, mod] & /@ Select[all, tfmSameModel[#["Model"], mod, -1] &], {}]]];

tfmSameModel[m1_, m2_, sign_] := m1["A"] == m2["A"] && m1["DMI"] == sign m2["DMI"] && m1["K"] == m2["K"] && m1["Field"] == m2["Field"];

tfmFlipDual[dual_, mod_] := Join[dual, <|"Model" -> mod, "Phi" -> -dual["Phi"], "Source" -> dual["Source"] <> " (\[Phi] \[RightArrow] -\[Phi])"|>];

(* the best certificate that passes the exact check *)
tfmChooseDual[spec_, given_] := Module[{mod = tfmModelOf[spec], cands, valid, vals},
  cands = Which[
    given === Automatic, Module[{cache = tfmShippedDuals[mod]},
      Join[cache, If[cache === {} && spec["Degree"] < 0 && mod["DMI"] != 0, tfmSearchDual[mod], {}], tfmSimpleDuals[mod, spec["Degree"]]]],
    StringQ[given], {tfmImportDual[given]},
    AssociationQ[given], {Join[<|"Gamma" -> 0, "Source" -> "given"|>, given, <|"Model" -> mod|>]},
    True, Throw[failure["BadModel", "\"DualField\" must be Automatic, a file written by the Python reference, or an association with \"Mesh\", \"Phi\", \"U\" (and \"Gamma\")."], "CertiField"]];
  If[given =!= Automatic && ! tfmSameModel[cands[[1]]["Model"], mod, 1], Throw[failure["BadModel", "The dual certificate belongs to a different model."], "CertiField"]];
  If[given =!= Automatic, cands = {Append[cands[[1]], "Model" -> mod]}];
  valid = Select[cands, tfmCheckDual];
  If[valid === {}, Throw[failure["Inconclusive", "No dual certificate passed the exact check."], "CertiField"]];
  vals = N[tfmLowerExact[#, spec["Degree"], spec["Area"]]] & /@ valid;
  valid[[First[Ordering[vals, -1]]]]];

(* ================================================================================================ *)
(* the dual search (heuristic, port of Reference/Magnetisation/dual_opt.py; the result is checked    *)
(* exactly like any other certificate)                                                               *)
(* ================================================================================================ *)

(* Gauss-Legendre rule on [0, 1] (Golub-Welsch) *)
tfmGaussRule[n_Integer] := tfmGaussRule[n] = Module[{b = Table[j/Sqrt[4. j^2 - 1], {j, n - 1}], vals, vecs, ord},
  {vals, vecs} = Eigensystem[Normal[SparseArray[{Band[{1, 2}] -> b, Band[{2, 1}] -> b}, {n, n}]]];
  ord = Ordering[vals];
  {(vals[[ord]] + 1)/2, vecs[[ord, 1]]^2}];

(* knots t = 1 - s, s = 2, 2 (3/4), 2 (3/4)^2, ... > sLast, and t = 1; phi = kappa on the last cell *)
tfmDualMesh[sLast_] := Module[{s = {2.}}, While[Last[s] 0.75 > sLast, AppendTo[s, Last[s] 0.75]]; Append[1. - s, 1.]];

(* 1/2 Integral u_-(t) dt and its gradient; unknowns: values and scaled slopes (1 - t) phi' at the knots but the last two *)
tfmDualFG[z_, d_] := Module[{n = d["n"], hk = d["hk"], sk = d["sk"], A = d["A"], k = d["k"], K = d["K"], h = d["h"], eps = d["eps"],
   H0, H1, H2, H3, dH0, dH1, dH2, dH3, val, slope, t, phi, dphi, a1, b1, P, R, V, B, Dd, good, sq, u, dPf, dPp, dRf, dRp, dBf, dBp, dDf, dDp,
   duf, dup, W, Gf, Gp, dval, dslope},
  {H0, H1, H2, H3} = d["H"]; {dH0, dH1, dH2, dH3} = d["dH"];
  val = Join[z[[1 ;; n - 1]], {k, k}];
  slope = Join[z[[n ;; 2 n - 2]]/sk[[1 ;; n - 1]], {0., 0.}];
  t = d["t"];
  phi = Outer[Times, val[[1 ;; n]], H0] + Outer[Times, hk slope[[1 ;; n]], H1] + Outer[Times, val[[2 ;; n + 1]], H2] + Outer[Times, hk slope[[2 ;; n + 1]], H3];
  dphi = (Outer[Times, val[[1 ;; n]], dH0] + Outer[Times, hk slope[[1 ;; n]], dH1] + Outer[Times, val[[2 ;; n + 1]], dH2] + Outer[Times, hk slope[[2 ;; n + 1]], dH3])/hk;
  a1 = (phi + k) + (1 + t) dphi; b1 = -(phi - k) + (1 - t) dphi;
  P = ((1 - t) a1)^2/4; R = ((1 + t) b1)^2/4;
  V = (1 - eps) (K (1 - t^2) + h (1 - t));
  B = 4 A V + P - R; Dd = B^2 - 16 A V P;
  good = UnitStep[Dd] (1 - UnitStep[-B]);
  sq = Sqrt[Ramp[Dd]];
  dPf = (1 - t)^2 a1/2; dPp = dPf (1 + t); dRf = -(1 + t)^2 b1/2; dRp = (1 + t)^2 b1 (1 - t)/2;
  dBf = dPf - dRf; dBp = dPp - dRp; dDf = 2 B dBf - 16 A V dPf; dDp = 2 B dBp - 16 A V dPp;
  (* where there is no root: a penalty, as in the reference *)
  u = good (B - sq)/(2 V) + (1 - good) (Abs[B]/(2 V) + 50. Ramp[-Dd]/(4 A V)^2 4 A + 10. Ramp[-B]/V);
  duf = good (dBf - dDf/(2 (sq + 10.^-300)))/(2 V) + (1 - good) (Sign[B] dBf/(2 V) - 200. A (1 - UnitStep[Dd]) dDf/(4 A V)^2 - 10. UnitStep[-B] dBf/V);
  dup = good (dBp - dDp/(2 (sq + 10.^-300)))/(2 V) + (1 - good) (Sign[B] dBp/(2 V) - 200. A (1 - UnitStep[Dd]) dDp/(4 A V)^2 - 10. UnitStep[-B] dBp/V);
  W = d["W"];
  Gf = W duf; Gp = W dup;
  dval = Join[Gf . H0 + (Gp . dH0)/hk, {0.}] + Join[{0.}, Gf . H2 + (Gp . dH2)/hk];
  dslope = Join[hk (Gf . H1) + Gp . dH1, {0.}] + Join[{0.}, hk (Gf . H3) + Gp . dH3];
  {Total[W u, 2], Join[dval[[1 ;; n - 1]], dslope[[1 ;; n - 1]]/sk[[1 ;; n - 1]]]}];

tfmDualData[mod_, sLast_, eps_] := Module[{tk = tfmDualMesh[sLast], gx, gw, hk, n},
  {gx, gw} = tfmGaussRule[10]; n = Length[tk] - 1; hk = Differences[tk];
  <|"n" -> n, "tk" -> tk, "hk" -> hk, "sk" -> 1. - tk, "A" -> N[mod["A"]], "k" -> N[mod["DMI"]], "K" -> N[mod["K"]], "h" -> N[mod["Field"]], "eps" -> eps,
    "t" -> Most[tk] + Outer[Times, hk, gx], "W" -> 0.5 Outer[Times, hk, gw],
    "H" -> {2 gx^3 - 3 gx^2 + 1, gx^3 - 2 gx^2 + gx, -2 gx^3 + 3 gx^2, gx^3 - gx^2},
    "dH" -> {6 gx^2 - 6 gx, 3 gx^2 - 4 gx + 1, -6 gx^2 + 6 gx, 3 gx^2 - 2 gx}|>];

(* mesh continuation: the last cell (where phi = kappa) is subdivided stage by stage, each stage starting from the last *)
tfmOptimiseDual[mod_, eps_ : 2.*^-3] := Module[{z = None, d, z0, m, res},
  Do[
    d = tfmDualData[mod, sLast, eps];
    z0 = Join[ConstantArray[d["k"], d["n"] - 1], ConstantArray[0., d["n"] - 1]];
    If[z =!= None, m = Length[z]/2; z0[[1 ;; m]] = z[[1 ;; m]]; z0[[d["n"] ;; d["n"] + m - 1]] = z[[m + 1 ;; 2 m]]];
    res = ndLBFGS[tfmDualFG[#, d] &, z0, 3000];
    z = res[[1]],
    {sLast, {2.*^-3, 2.*^-4, 2.*^-5}}];
  {d, z}];

tfmHermiteCell[t0_, t1_, f0_, d0_, f1_, d1_] := Module[{hh = t1 - t0, s},
  s = (tfmT - t0)/hh;
  CoefficientList[Expand[f0 (2 s^3 - 3 s^2 + 1) + d0 hh (s^3 - 2 s^2 + s) + f1 (-2 s^3 + 3 s^2) + d1 hh (s^3 - s^2)], tfmT, 4]];

(* on every cell a polynomial u of degree 4 slightly above the smaller root, checked exactly cell by cell *)
tfmFitU[mod_, mesh_, phis_] := Module[{n = Length[mesh] - 1, A = mod["A"], k = mod["DMI"], K = mod["K"], h = mod["Field"], us = {}, bad = 0,
   t0, t1, nodes, ph, dph, a1, b1, P, R, V, B, Dd, um, up, ok, target, c, upoly, v0, cheb, vand},
  cheb = 0.5 (1 - Cos[Pi (Range[0, 4] + 0.5)/5]);
  vand = Outer[Power, cheb, Range[0, 4]];
  Do[
    t0 = mesh[[i]]; t1 = mesh[[i + 1]]; ok = False;
    If[i == n,
      v0 = h + K (1 + If[K >= 0, t0, 1]);
      If[v0 > 0 && 2 k^2/v0 <= 4 A, upoly = {k^2/v0, -k^2/v0}; ok = True],
      nodes = N[t0] + N[t1 - t0] cheb;
      ph = polyAt[N[phis[[i]]], nodes]; dph = polyAt[N[polyD[phis[[i]]]], nodes];
      a1 = (ph + N[k]) + (1 + nodes) dph; b1 = -(ph - N[k]) + (1 - nodes) dph;
      P = ((1 - nodes) a1)^2/4; R = ((1 + nodes) b1)^2/4; V = N[K] (1 - nodes^2) + N[h] (1 - nodes);
      B = 4 N[A] V + P - R; Dd = B^2 - 16 N[A] V P;
      If[Min[Dd] > 0 && Min[B] > 0 && Min[V] > 0,
        um = (B - Sqrt[Dd])/(2 V); up = (B + Sqrt[Dd])/(2 V);
        Do[
          target = If[mg > 0, um (1 + mg) + 10.^-3 mg Max[um], um + (-mg) (up - um)];
          c = Rationalize[LinearSolve[vand, target], 0];
          upoly = CoefficientList[Expand[c . ((tfmT - t0)/(t1 - t0))^Range[0, 4]], tfmT, 5];
          If[tfmCellOK[tfmCellPolys[<|"Model" -> mod, "Phi" -> {phis[[i]]}, "U" -> {upoly}, "Gamma" -> 0|>, 1], A, t0, t1], ok = True; Break[]],
          {mg, {10.^-6, 10.^-5, 10.^-4, 10.^-3, 10.^-2, -0.02, -0.1, -0.3, -0.5}}]]];
    If[! ok, bad++; upoly = {4 A}];
    AppendTo[us, upoly],
    {i, n}];
  {us, bad}];

(* the optimiser's field, rounded to rationals; if some cell fails, damped towards phi = kappa *)
tfmBuildDual[mod_, d_, z_] := Module[{n = d["n"], k = mod["DMI"], kf = d["k"], sk = d["sk"], mesh, zv, zs, vals, slopes, phis, us, bad, best = None},
  mesh = Rationalize[d["tk"], 0]; mesh[[1]] = -1; mesh[[-1]] = 1;
  Do[
    zv = kf + (1 - th) (z[[1 ;; n - 1]] - kf);
    zs = (1 - th) z[[n ;; 2 n - 2]]/sk[[1 ;; n - 1]];
    vals = Join[Rationalize[zv, 0], {k, k}]; slopes = Join[Rationalize[zs, 0], {0, 0}];
    phis = Table[tfmHermiteCell[mesh[[i]], mesh[[i + 1]], vals[[i]], slopes[[i]], vals[[i + 1]], slopes[[i + 1]]], {i, n}];
    {us, bad} = tfmFitU[mod, mesh, phis];
    If[best === None || bad < best[[2]], best = {<|"Model" -> mod, "Mesh" -> mesh, "Phi" -> phis, "U" -> us, "Gamma" -> 0, "Source" -> "search"|>, bad}];
    If[bad == 0, Break[]],
    {th, {0, 10^-4, 10^-3, 3 10^-3, 10^-2, 3 10^-2, 10^-1}}];
  If[best[[2]] == 0, {best[[1]]}, {}]];

(* for degree < 0 and kappa != 0, with V > 0 on (-1, 1); the search runs for |kappa| and phi -> -phi gives the field
   for -kappa.  Remembered for the session. *)
tfmSearchDual[mod_] := tfmSearchDual[mod] = Module[{m = Append[mod, "DMI" -> Abs[mod["DMI"]]], d, z, out},
  If[! (mod["Field"] > 0 && mod["Field"] + 2 mod["K"] > 0), Return[{}, Module]];
  {d, z} = Quiet[tfmOptimiseDual[m]];
  If[! VectorQ[z, Internal`RealValuedNumericQ], Return[{}, Module]];
  out = tfmBuildDual[m, d, z];
  If[mod["DMI"] < 0, tfmFlipDual[#, mod] & /@ out, out]];

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
tfmUpperBound[mod_, x1_, x2_, Mn_, extra_ : 0] := ndQuiet@Module[
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
  If[mod["Field"] != 0 || extra =!= 0,
    m3 = #[[3]] & /@ {B["00"], B["10"], B["01"], B["11"]};
    m0 = ivK[ivSum[m3], 1/4];
    m1 = ivSum[Flatten[Table[ivK[ndMul[{B["00"], B["01"], B["10"], B["11"]}[[1 + 2 i + j, 3]], c[[k, l]]], $tfmG12[[i + 1, k]] $tfmG12[[j + 1, l]]], {i, 0, 1}, {j, 0, 1}, {k, 3}, {l, 3}], 3]];
    neg = Ramp[-Fold[vMin, #[[1]] & /@ m3]];
    gz = vMax[ndDiv[ivSqr[ndSub[{1., 1.}, a]], ivK[a, 2]][[2]], ndDiv[ivSqr[ndSub[{1., 1.}, b]], ivK[b, 2]][[2]]];
    ze = ivSum[{{1., 1.}, ivK[m0, -3/2], ivK[m1, 1/2], ndMul[ivPt[neg], ivPt[gz]]}];
    total = ndAdd[total, ndMul[ndMul[ze, area], ndAdd[ivC[mod["Field"]], {0. alo, 0. alo + extra}]]]];
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

(* nodal vectors of an axisymmetric field about the centre of the domain, e3 on the boundary and outside the inscribed
   circle: degree -1: (-sin(theta) y/r, sin(theta) x/r, cos(theta)) for kappa >= 0 and its image under
   (m1, m2) -> (-m1, -m2) for kappa < 0; degree +1: (sin(theta) x/r, -sin(theta) y/r, cos(theta)) *)
tfmRadialNodes[spec_, x1_, x2_, theta_, Ltex_] := Module[{dom = spec["Domain"], cx, cy, L = N[Ltex], X, Y, r, th, s, inside, Mn, m1, m2},
  cx = N[Mean[dom[[1]]]]; cy = N[Mean[dom[[2]]]];
  X = Transpose[ConstantArray[x1 - cx, Length[x2]]]; Y = ConstantArray[x2 - cy, Length[x1]];
  r = Sqrt[X^2 + Y^2];
  th = ArrayReshape[tfmEvalList[theta, Flatten[Clip[r, {0., L}]]], Dimensions[r]];
  s = Sin[th]/(r + UnitStep[-r]);
  {m1, m2} = Which[spec["Degree"] > 0, {s X, -s Y}, spec["DMI"] < 0, {s Y, -s X}, True, {-s Y, s X}];
  Mn = Transpose[{m1, m2, Cos[th]}, {3, 1, 2}];
  inside = 1 - UnitStep[r - L];
  inside[[1]] = inside[[-1]] = 0 inside[[1]]; inside[[All, 1]] = 0; inside[[All, -1]] = 0;
  Mn = Mn inside + Transpose[{0. inside, 0. inside, 1. - inside}, {3, 1, 2}];
  Developer`ToPackedArray[N[Mn]]];

tfmEvalList[f_, list_] := Module[{v = Quiet[f[list]]}, If[VectorQ[v, Internal`RealValuedNumericQ] && Length[v] == Length[list], N[v], N[f /@ list]]];

tfmAutoGrading[spec_, theta_, n_, Ltex_] := Module[{L, rr, core, f},
  L = N[Min[(spec["Domain"][[1, 2]] - spec["Domain"][[1, 1]])/2, (spec["Domain"][[2, 2]] - spec["Domain"][[2, 1]])/2]];
  rr = Subdivide[0., N[Ltex], 4000];
  core = Max[rr[[Replace[Flatten[Position[UnitStep[Pi/2 - tfmEvalList[theta, rr]], 1, 1, 1]], {} -> {2}][[1]]]], L/1000];
  (* smallest cell about core/40 *)
  f[beta_?NumericQ] := L beta/(n Sinh[beta]) - core/40;
  If[f[0.01] <= 0, 0., Quiet@Check[beta /. FindRoot[f[beta], {beta, 2., 0.01, 30.}], 4.]]];

(* a Belavin-Polyakov profile of small radius, cut off at the inscribed circle (energy just above 8 pi A) *)
tfmSmallBP[spec_, L_] := Module[{stiff = N[spec["Field"] + 2 Abs[spec["K"]]], R},
  R = Min[N[L]/40, If[stiff > 0, 0.05/Sqrt[stiff], Infinity]];
  With[{RR = R, LL = N[L]}, Function[r, 2 ArcTan[r, RR (1 - (r/LL)^2)^2]]]];

tfmCompetitor[spec_, cand_, frac_ : 1] := Module[{n = spec["Cells"]/2, theta, beta, x1, x2, Mn, f, L},
  L = frac Min[(spec["Domain"][[1, 2]] - spec["Domain"][[1, 1]])/2, (spec["Domain"][[2, 2]] - spec["Domain"][[2, 1]])/2];
  theta = Which[
    cand === Automatic,
      If[! MemberQ[{-1, 1}, spec["Degree"]], Throw[failure["BadCandidate", "The automatic competitor has degree -1 or +1; for other degrees give a field Function[{x, y}, {m1, m2, m3}]."], "CertiField"]];
      (* degree -1 with a Dzyaloshinskii-Moriya term: the radial minimiser; otherwise a small Belavin-Polyakov map *)
      If[spec["Degree"] == -1 && spec["DMI"] != 0, tfmProfile[Append[tfmModelOf[spec], "DMI" -> Abs[spec["DMI"]]], L], tfmSmallBP[spec, L]],
    True, None];
  f = If[theta === None && ! MatchQ[cand, _Function | _InterpolatingFunction], Throw[failure["BadCandidate", "The candidate must be Automatic, a profile Function[r, theta] or a field Function[{x, y}, {m1, m2, m3}]."], "CertiField"], cand];
  If[theta === None && NumericQ[Quiet[f[0.5]]], theta = f];
  beta = Replace[spec["Grading"], Automatic :> If[theta =!= None, tfmAutoGrading[spec, theta, n, L], 3.]];
  x1 = tfmGradedAxis[spec["Domain"][[1, 1]], spec["Domain"][[1, 2]], n, beta];
  x2 = tfmGradedAxis[spec["Domain"][[2, 1]], spec["Domain"][[2, 2]], n, beta];
  Mn = If[theta =!= None, tfmRadialNodes[spec, x1, x2, theta, L],
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
  If[spec["Thickness"] > 0, Return[tfmCertifyStray[spec, cand], Module]];
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

(* ---- the stray field (docs/MAGNETISATION.md, Section 10): E[m] - E[e3] for the film of thickness t ----
   lower: t E2D[A - t lam/8, kappa, K - 1/2 + t/(8 lam), h - t/(4 lam)]  (any lam > 0),
   upper: t E2D[A, kappa, K - 1/2, h + 1 - clow][m_h] + t^2/4 Sqrt[2 Int |grad m_h|^2 Int (1 - m3^2)],
   clow = 1 - t/(2 rho) <= the mean of the 2D Poisson kernels over the thickness, rho the distance of the texture to the edge *)

tfmSqrtUp[q_] := Module[{s = fHi[Sqrt[N[q, 30]]]}, While[exactNumber[s]^2 < q, s = nextUp[s]]; exactNumber[s]];

(* per cell an upper bound for 1 - c on the cell: min(1, t/(2 rho)), rho the distance of the cell to the edge *)
tfmEdgeFactor[spec_, t_, x1_, x2_] := Module[{dom = spec["Domain"], dx, dy, rho},
  dx = MapThread[Min[#1 - dom[[1, 1]], dom[[1, 2]] - #2] &, {exactNumber /@ Most[x1], exactNumber /@ Rest[x1]}];
  dy = MapThread[Min[#1 - dom[[2, 1]], dom[[2, 2]] - #2] &, {exactNumber /@ Most[x2], exactNumber /@ Rest[x2]}];
  rho = Outer[Min, dx, dy];
  Developer`ToPackedArray[Map[If[# <= 0, 1., fHi[Min[1, t/(2 #)]]] &, rho, {2}]]];

(* RIGOROUS: an exact rational >= Integral (div m')^2, m' = (m1, m2), by interval evaluation on every cell *)
tfmDivBound[x1_, x2_, Mn_] := ndQuiet@Module[{B = tfmCells[x1, x2, Mn], c, alo, bhi, a, b, hull, M, dxi, deta, d1, d2, n1, n2, num, isd, dv},
  c = tfmDCoefficients[B, {1, 2, 3}];
  {alo, bhi} = tfmRange[Flatten[c, 1]];
  If[Min[alo] <= 0, Throw[failure["Inconclusive", "The interpolant of the candidate may vanish."], "CertiField"]];
  a = ivPt[alo]; b = ivPt[bhi];
  hull[ivs_List] := {Fold[vMin, #[[1]] & /@ ivs], Fold[vMax, #[[2]] & /@ ivs]};
  M = Table[hull[{B["00"][[j]], B["10"][[j]], B["01"][[j]], B["11"][[j]]}], {j, 3}];
  dxi = Table[hull[{ndSub[B["10"][[j]], B["00"][[j]]], ndSub[B["11"][[j]], B["01"][[j]]]}], {j, 3}];
  deta = Table[hull[{ndSub[B["01"][[j]], B["00"][[j]]], ndSub[B["11"][[j]], B["10"][[j]]]}], {j, 3}];
  d1 = ndDiv[#, B["h1"]] & /@ dxi; d2 = ndDiv[#, B["h2"]] & /@ deta;
  n1 = tfmDot[M, d1]; n2 = tfmDot[M, d2];
  num = ndSub[ndAdd[d1[[1]], d2[[2]]], ndDiv[ndAdd[ndMul[M[[1]], n1], ndMul[M[[2]], n2]], {alo, bhi}]];
  isd = {ndDiv[{2., 2.}, ndAdd[b, {1., 1.}]][[1]], ndDiv[ndAdd[a, {1., 1.}], ivK[a, 2]][[2]]};    (* 2/(1 + D) <= 1/Sqrt[D] <= (1 + D)/(2 D) *)
  dv = ndMul[num, isd];
  exactTotal[ndMul[ivSqr[dv], ndMul[B["h1"], B["h2"]]][[2]]]];

tfmCLow[spec_, t_, x1_, x2_, Mn_] := Module[{dom = spec["Domain"], nz, cm, dx, dy, rows, cols, rho},
  nz = 1 - Map[Boole[# === {0., 0., 1.}] &, Mn, {2}];
  cm = Unitize[nz[[1 ;; -2, 1 ;; -2]] + nz[[2 ;;, 1 ;; -2]] + nz[[1 ;; -2, 2 ;;]] + nz[[2 ;;, 2 ;;]]];
  rows = Flatten[Position[Total[cm, {2}], _?Positive]]; cols = Flatten[Position[Total[cm], _?Positive]];
  If[rows === {}, Return[1, Module]];
  dx = MapThread[Min[#1 - dom[[1, 1]], dom[[1, 2]] - #2] &, {exactNumber /@ Most[x1], exactNumber /@ Rest[x1]}];
  dy = MapThread[Min[#1 - dom[[2, 1]], dom[[2, 2]] - #2] &, {exactNumber /@ Most[x2], exactNumber /@ Rest[x2]}];
  rho = Min[dx[[rows]], dy[[cols]]];
  If[rho <= 0, 0, Max[0, 1 - t/(2 rho)]]];

(* heuristic: lam balancing t^2 lam/8 Int |grad m|^2 against t^2/(8 lam) Int (1 - m3)^2 on the competitor *)
tfmLambda[spec_, X_, x1_, x2_, Mn_] := Module[{t = spec["Thickness"], m3, w2, area, yw, lam},
  m3 = Mn[[All, All, 3]]/Sqrt[Total[Mn^2, {3}]];
  w2 = (1 - m3)^2;
  area = Outer[Times, Differences[x1], Differences[x2]];
  yw = Total[area (w2[[1 ;; -2, 1 ;; -2]] + w2[[2 ;;, 1 ;; -2]] + w2[[1 ;; -2, 2 ;;]] + w2[[2 ;;, 2 ;;]])/4, 2];
  lam = If[N[X] > 0 && yw > 0, Sqrt[yw/N[X]], 1.];
  lam = Clip[lam, {10.^-3, 2. spec["A"]/t}];
  Rationalize[lam, lam 10^-3]];

tfmShifted[spec_, lam_] := With[{t = spec["Thickness"]},
  Join[spec, <|"A" -> spec["A"] - t lam/8, "K" -> spec["K"] - 1/2 + t/(8 lam), "Field" -> spec["Field"] - t/(4 lam), "Thickness" -> 0|>]];

tfmStrayUpper[spec_, x1_, x2_, Mn_] := Module[{t = spec["Thickness"], X, Y, dv, clow, ubLoc},
  X = tfmUpperBound[<|"A" -> 1, "DMI" -> 0, "K" -> 0, "Field" -> 0|>, x1, x2, Mn];
  Y = tfmUpperBound[<|"A" -> 0, "DMI" -> 0, "K" -> 1, "Field" -> 0|>, x1, x2, Mn];
  dv = Min[tfmDivBound[x1, x2, Mn], 2 X];
  clow = tfmCLow[spec, t, x1, x2, Mn];
  ubLoc = tfmUpperBound[<|"A" -> spec["A"], "DMI" -> spec["DMI"], "K" -> spec["K"] - 1/2, "Field" -> spec["Field"]|>, x1, x2, Mn, tfmEdgeFactor[spec, t, x1, x2]];
  {t ubLoc + t^2/4 tfmSqrtUp[dv Y], X, Y, clow}];

tfmCertifyStray[spec_, cand_] := Quiet[Module[{t = spec["Thickness"], t0 = AbsoluteTime[], x1, x2, Mn, beta, deg, ub, X, Y, clow, lam, shifted, dual, lb},
  {x1, x2, Mn, beta} = tfmCompetitor[Join[spec, <|"K" -> spec["K"] - 1/2|>], cand, 3/4];
  deg = tfmDegree[x1, x2, Mn];
  If[deg[[1]] =!= spec["Degree"],
    Throw[failure["Inconclusive", "The degree of the competitor could not be certified to be " <> ToString[spec["Degree"]] <> "."], "CertiField"]];
  {ub, X, Y, clow} = tfmStrayUpper[spec, x1, x2, Mn];
  lam = tfmLambda[spec, X, x1, x2, Mn];
  shifted = tfmShifted[spec, lam];
  dual = tfmChooseDual[shifted, Automatic];
  lb = t tfmLowerExact[dual, spec["Degree"], spec["Area"]];
  <|"LowerBound" -> lb, "UpperBound" -> ub, "Seconds" -> AbsoluteTime[] - t0, "Dual" -> dual, "Lambda" -> lam, "CLow" -> clow,
    "BoundCoefficients" -> t tfmBoundCoefficients[dual, spec["Degree"], spec["Area"]], "Grid" -> {x1, x2}, "Grading" -> beta, "Nodal" -> Mn,
    "DegreeEnclosure" -> deg[[2]], "Spec" -> spec|>], {General::munfl}];

tfmDescription[params_, res_] := Module[{sp = res["Spec"], d = Max[displayDigits[res], 8], stray, en},
  stray = sp["Thickness"] > 0;
  en = If[stray, "E[m] - E[e3]", "E[m]"];
  <|"Units" -> If[stray, "\[Mu]0 Ms^2 times length^3 (lengths in the units of the domain)", "energy units of the model"],
    "Theorem" -> Row[{"Every m in H^1 of the rectangle ", exactString[sp["Domain"]], " with values in S^2, m = (0, 0, 1) on the boundary and degree ", sp["Degree"],
      " has ", en, " \[GreaterEqual] ", lowerString[res["LowerBound"], d], "; the competitor, of certified degree ", sp["Degree"], ", has ", en, " \[LessEqual] ", upperString[res["UpperBound"], d],
      If[stray, Row[{" (stray field of the film of thickness ", exactString[sp["Thickness"]], " included; e3 is the uniformly magnetised film)"}], ""], "."}],
    "Energy" -> Row[{"E[m] = ", If[stray, Row[{exactString[sp["Thickness"]], " "}], ""], "\[Integral] ", exactString[sp["A"]], " |\[Del]m|^2 + ", exactString[sp["DMI"]], " m\[CenterDot](\[Del]\[Times]m) + ", exactString[sp["K"]],
      " (1 - m3^2) + ", exactString[sp["Field"]], " (1 - m3) dx",
      If[stray, " + (1/2) \[Integral] |\[Del]U|^2 over \[DoubleStruckCapitalR]^3, \[CapitalDelta]U = div(m \[Chi]), \[Chi] the film \[CapitalOmega] \[Times] (0, t), m independent of x3", ""]}],
    "AdmissibleClass" -> Row[{"all m in H^1 of the rectangle ", exactString[sp["Domain"]], " with |m| = 1, m = (0, 0, 1) on the boundary (trace) and degree ",
      sp["Degree"], " (the degree is the integer (1/4\[Pi]) \[Integral] m\[CenterDot](\[PartialD]1 m \[Times] \[PartialD]2 m))"}],
    "LowerBoundMethod" -> Row[{"null Lagrangians of degree and divergence type with a dual field (\[Phi], u, \[CapitalGamma]) of piecewise polynomials in m3 on ",
      Length[res["Dual", "Mesh"]] - 1, " intervals (", res["Dual", "Source"], "), checked in exact rational arithmetic"}],
    "UpperBoundMethod" -> Row[{"normalised bilinear interpolant on a graded ", Length[res["Grid"][[1]]] - 1, "\[Times]", Length[res["Grid"][[2]]] - 1,
      " grid; energy and degree enclosed with outward-rounded arithmetic"}],
    "TrustedBase" -> "IEEE-754 machine arithmetic (+, \[Minus], \[Times], \[Divide]), Wolfram Language exact rational arithmetic, and docs/MAGNETISATION.md (Lemmas 1\[Dash]3, Theorem 4, Proposition 5; proofs checked by one author only)",
    "NotClaimed" -> If[stray, "that the competitor is a minimiser; magnetisations that vary through the thickness are not in the admissible class",
      "that the competitor is a minimiser, or that the axisymmetric skyrmion minimises the energy in its degree class; the stray field is not included"],
    "References" -> "docs/MAGNETISATION.md; Belavin\[Dash]Polyakov (1975); D\[ODoubleDot]ring\[Dash]Melcher (2017); Li\[Dash]Melcher, arXiv:1711.07717."|>];

verifyThinFilmMagnet[a_Association] := Catch[Module[{c = a["Certificate"], spec, dual, lb, ub, deg, x1, x2, t},
  spec = tfmSpec[a["Parameters"]]; t = spec["Thickness"];
  dual = c["Dual"];
  If[! tfmSameModel[dual["Model"], tfmModelOf[If[t > 0, tfmShifted[spec, c["Lambda"]], spec]], 1] || ! tfmCheckDual[dual], Throw[False, "tfmv"]];
  lb = If[t > 0, t, 1] tfmLowerExact[dual, spec["Degree"], spec["Area"]];
  {x1, x2} = c["Grid"];
  If[{x1[[1]], x1[[-1]], x2[[1]], x2[[-1]]} != N[Flatten[spec["Domain"]]], Throw[False, "tfmv"]];
  deg = tfmDegree[x1, x2, c["Nodal"]];
  If[deg[[1]] =!= spec["Degree"], Throw[False, "tfmv"]];
  ub = If[t > 0, tfmStrayUpper[spec, x1, x2, c["Nodal"]][[1]], tfmUpperBound[tfmModelOf[spec], x1, x2, c["Nodal"]]];
  lb >= exactNumber[a["LowerBound"]] && ub <= exactNumber[a["UpperBound"]]], "tfmv"];

End[];
