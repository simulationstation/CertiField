(* ::Package:: *)

(* ::Title:: *)
(*CertiField: certified affine minorants of h_kappa*)

(* ::Text:: *)
(*h_kappa(D) = 2 Sqrt[(1 + 2D)(D^2 + 2D + kappa)].  A triple (alpha, beta, kappa) >= 0 is a minorant if*)
(*alpha + beta D <= h_kappa(D) for all D >= 0.  Since both sides are nonnegative, this holds iff*)
(**)
(*    P(D) = 8 D^3 + (20 - beta^2) D^2 + (8 + 8 kappa - 2 alpha beta) D + 4 kappa - alpha^2 >= 0   on [0, oo).*)
(**)
(*The table is built from tangent lines of h_kappa (heuristic: placement only).  Every entry is then certified here*)
(*by an outward-rounded evaluation of P at its interior local minimum.  Row kappa = 0 is (0, 6, 0), certified by the*)
(*identity 4 q_0(D) - 36 D^2 = 8 D (D - 1)^2.*)
(*Reference: note_levelset.tex, "Minorants" (M1)-(M3).*)

Begin["CertiField`Private`"];

(* ---- h_kappa and its tangent geometry (machine precision; used only to choose table entries) ---- *)
hq[k_, d_] := (1 + 2 d) (d d + 2 d + k);
hh[k_, d_] := 2 Sqrt[hq[k, d]];
hDer[k_, d_] := (2 (d d + 2 d + k) + (1 + 2 d) (2 d + 2))/Sqrt[hq[k, d]];
hIntercept[k_, d_] := hh[k, d] - d hDer[k, d];

$kConvexLo = (3 - Sqrt[5.])/2;    (* h_kappa is convex on [0, oo) iff kappa in [kConvexLo, kConvexHi] *)
$kConvexHi = (3 + Sqrt[5.])/2;

bisectRoot[f_, a0_, b0_] := Module[{a = a0, b = b0, fa = f[a0], m, fm},
  Do[
    m = 0.5 (a + b);
    If[m == a || m == b, Break[]];
    fm = f[m];
    If[(fm > 0) === (fa > 0), a = m; fa = fm, b = m],
    {300}];
  0.5 (a + b)];

dInflection[k_] := Module[{n0 = -4 (k k - 3 k + 1), nf, hi = 1.},
  If[n0 >= 0, 0.,
    nf = Function[d, 12 d^4 + 40 d^3 + 24 (1 + k) d d + 24 k d + n0];
    While[nf[hi] < 0, hi *= 2];
    bisectRoot[nf, 0., hi]]];

(* tangent point of the line through (0, 2 Sqrt[kappa]) *)
dTangent[k_] := Which[
  $kConvexLo <= k <= $kConvexHi, 0.,
  k <= 0, 1.,
  True, Module[{di = dInflection[k], g, hi},
    g = Function[d, hIntercept[k, d] - 2 Sqrt[k]];
    hi = Max[2 di, 1.*^-6];
    While[g[hi] > 0, hi *= 2];
    bisectRoot[g, di, hi]]];

(* point where the tangent passes through the origin *)
dStar[k_] := Module[{lo, hi},
  lo = Max[If[$kConvexLo <= k <= $kConvexHi, 1.*^-12, dTangent[k]], 1.*^-12];
  hi = Max[2 lo, 1.];
  While[hIntercept[k, hi] > 0, hi *= 2];
  If[hIntercept[k, lo] <= 0, lo, bisectRoot[Function[d, hIntercept[k, d]], lo, hi]]];

(* ---- vectorised certification: 1 where (alpha, beta, kappa) is a certified minorant, else 0 ---- *)
certifyMinorantMask[al_, be_, k_] := Module[
  {bb, ab, c2, c1, c0, ok0, monotone, disc, noCrit, sq, d2, d2neg, X, p, pPos},
  bb = {rDn[be be], rUp[be be]};
  ab = {rDn[al be], rUp[al be]};
  c2 = ivSub[{20., 20.}, bb];
  c1 = {rDn[rDn[8. + 8. k] - 2. ab[[2]]], rUp[rUp[8. + 8. k] - 2. ab[[1]]]};   (* 8k and 2ab are exact scalings *)
  c0 = {rDn[4. k - rUp[al al]], rUp[4. k - rDn[al al]]};
  ok0 = geMask[al, 0.] geMask[be, 0.] geMask[k, 0.] gtMask[c0[[1]], 0.];
  monotone = geMask[c2[[1]], 0.] geMask[c1[[1]], 0.];
  disc = ivSub[ivMul[c2, c2], ivMul[{24., 24.}, c1]];
  noCrit = ltMask[disc[[2]], 0.];
  sq = {sqrtLo[disc[[1]]], sqrtHi[disc[[2]]]};
  d2 = ivDivPos[ivAdd[{-c2[[2]], -c2[[1]]}, sq], {24., 24.}];
  d2neg = leMask[d2[[2]], 0.];
  X = {Ramp[d2[[1]]], Ramp[d2[[2]]]};
  p = ivAdd[ivMul[ivAdd[ivMul[ivAdd[ivMul[{8., 8.}, X], c2], X], c1], X], c0];
  pPos = gtMask[p[[1]], 0.];
  ok0 (monotone + (1 - monotone) (noCrit + (1 - noCrit) (d2neg + (1 - d2neg) pPos)))];

(* identity behind row kappa = 0, checked symbolically once *)
$rowZeroIdentity := (Expand[4 (1 + 2 #) (#^2 + 2 #) - 36 #^2 - 8 # (# - 1)^2] === 0) &[Unique["d"]];

minorantRow[k_, dg_] := Module[{dt, ds, lo, hi, de, b, a, acc, alpha, beta, s, aa, bb, ok, fresh, fb},
  dt = dTangent[k]; ds = dStar[k];
  lo = dt (1 + 1.*^-9); hi = ds (1 - 1.*^-9);
  de = If[hi > lo, Clip[dg, {lo, hi}], ConstantArray[lo, Length[dg]]];
  b = hDer[k, de]; a = hh[k, de] - de b;
  acc = ConstantArray[0, Length[dg]];
  alpha = ConstantArray[0., Length[dg]]; beta = ConstantArray[0., Length[dg]];
  Do[
    aa = Ramp[a] (1 - s); bb = b (1 - s);
    ok = certifyMinorantMask[aa, bb, k];
    fresh = ok (1 - acc);
    alpha = vWhere[fresh, aa, alpha]; beta = vWhere[fresh, bb, beta];
    acc = acc + fresh,
    {s, {1.*^-9, 1.*^-7, 1.*^-5, 1.*^-3}}];
  (* trivial minorant 2 Sqrt[kappa] (1 - 10^-9) where nothing else certified *)
  fb = 2 Sqrt[k] (1 - 1.*^-9);
  If[certifyMinorantMask[{fb}, {0.}, k] =!= {1}, Throw[$Failed, "CertiFieldTable"]];
  alpha = vWhere[acc, alpha, fb]; beta = vWhere[acc, beta, 0.];
  {alpha, beta, Length[dg] - Total[acc]}];

(* Build the minorant table and certify every entry.  Parameters as in levelset_cert.build_table. *)
buildMinorantTable[nK_Integer, nD_Integer, kmin_?NumericQ, kmax_: 400., dmin_: 1.*^-4, dmax_: 100.] := Module[
  {KN, DG, rows, AL, BE, fallbacks},
  If[! TrueQ[$rowZeroIdentity], Throw[$Failed, "CertiFieldTable"]];
  KN = Developer`ToPackedArray[N@Join[{0.}, Exp[Subdivide[Log[N[kmin]], Log[N[kmax]], nK - 1]]]];
  DG = Developer`ToPackedArray[N@Join[{0.}, Exp[Subdivide[Log[N[dmin]], Log[N[dmax]], nD - 1]]]];
  rows = Table[minorantRow[KN[[n]], DG], {n, 2, Length[KN]}];
  AL = Developer`ToPackedArray[Join[{ConstantArray[0., Length[DG]]}, rows[[All, 1]]]];
  BE = Developer`ToPackedArray[Join[{ConstantArray[6., Length[DG]]}, rows[[All, 2]]]];
  fallbacks = Total[rows[[All, 3]]];
  makeMinorantTable[KN, DG, AL, BE, <|"Source" -> "built and certified in Wolfram Language",
    "Parameters" -> {nK, nD, kmin, kmax, dmin, dmax}, "Fallbacks" -> fallbacks|>]];

(* Wrap raw arrays as a table, certifying every entry (this is what makes externally supplied tables safe to use) *)
makeMinorantTable[KN_, DG_, AL_, BE_, info_Association] := Module[{nKt = Length[KN], nDt = Length[DG], HH, NU, okRows},
  If[! (VectorQ[KN, Developer`MachineRealQ] && OrderedQ[KN] && KN[[1]] == 0. && Dimensions[AL] === {nKt, nDt} &&
      Dimensions[BE] === {nKt, nDt}), Throw[$Failed, "CertiFieldTable"]];
  (* row 0 must be exactly (0, 6) at kappa = 0 *)
  If[! (Union[AL[[1]]] === {0.} && Union[BE[[1]]] === {6.}), Throw[$Failed, "CertiFieldTable"]];
  okRows = Table[Total[certifyMinorantMask[AL[[n]], BE[[n]], KN[[n]]]], {n, 2, nKt}];
  If[Total[okRows] =!= (nKt - 1) nDt, Throw[$Failed, "CertiFieldTable"]];
  HH = AL + BE (ConstantArray[DG, nKt]);
  NU = Differences[HH]/Differences[KN];
  <|"KN" -> KN, "DG" -> DG, "ALf" -> Flatten[AL], "BEf" -> Flatten[BE], "NUf" -> Flatten[NU],
    "nK" -> nKt, "nD" -> nDt, "ld0" -> Log[DG[[2]]], "lstep" -> (Log[DG[[-1]]] - Log[DG[[2]]])/(nDt - 2),
    "Certified" -> True, "Entries" -> nKt nDt, "Info" -> info|>];

iMin[a_, b_] := a - Ramp[a - b];     (* exact elementwise minimum of integer arrays *)

(* ---- table queries (vectorised); indices d are 0-based columns, as in levelset_cert2.Table ---- *)
tabDIndex[tab_, D_] := Clip[Round[(Log[vMax[D, 1.*^-300]] - tab["ld0"])/tab["lstep"]] + 1, {0, tab["nD"] - 1}];

tabRowIndex[tab_, kap_] := Clip[searchRight[tab["KN"], kap] - 1, {0, tab["nK"] - 2}];

(* heuristic (floating point) convex combination: used only to choose multipliers *)
tabMinor[tab_, kap_, d_] := Module[{n = tabRowIndex[tab, kap], k0, k1, th, i0, i1, nD = tab["nD"]},
  k0 = tab["KN"][[n + 1]]; k1 = tab["KN"][[n + 2]];
  th = Clip[(kap - k0)/(k1 - k0), {0., 1.}];
  i0 = n nD + d + 1; i1 = (n + 1) nD + d + 1;
  {(1 - th) tab["ALf"][[i0]] + th tab["ALf"][[i1]], (1 - th) tab["BEf"][[i0]] + th tab["BEf"][[i1]], (1 - th) k0 + th k1}];

(* rigorous: enclosures of the convex combination of two neighbouring certified rows (M2) *)
tabMinorIv[tab_, kap_, d_] := Module[{n = tabRowIndex[tab, kap], k0, k1, th, i0, i1, nD = tab["nD"], comb, al, be, kh},
  k0 = tab["KN"][[n + 1]]; k1 = tab["KN"][[n + 2]];
  th = ivPt[Clip[(kap - k0)/(k1 - k0), {0., 1.}]];
  i0 = n nD + d + 1; i1 = (n + 1) nD + d + 1;
  comb[a0_, a1_] := ivAdd[ivPt[a0], ivMul[ivSub[ivPt[a1], ivPt[a0]], th]];
  al = comb[tab["ALf"][[i0]], tab["ALf"][[i1]]];
  be = comb[tab["BEf"][[i0]], tab["BEf"][[i1]]];
  kh = comb[k0, k1];
  {Ramp[al[[1]]], Ramp[be[[1]]], {Ramp[kh[[1]]], kh[[2]]}}];

(* heuristic inverse of d h**/d kappa along a column (11 bisection steps, as in levelset_cert2.Table.kinv) *)
tabKInv[tab_, nuT_, d_] := Module[{nint = tab["nK"] - 1, lo, hi, mid, cond, nD = tab["nD"]},
  lo = ConstantArray[0, Length[nuT]]; hi = ConstantArray[nint, Length[nuT]];
  Do[
    mid = Quotient[lo + hi, 2];
    cond = leMask[tab["NUf"][[Clip[mid, {0, nint - 1}] nD + d + 1]], nuT];
    hi = vWhere[cond, mid, hi];
    lo = vWhere[cond, lo, iMin[mid + 1, hi]],
    {11}];
  tab["KN"][[Clip[lo, {0, tab["nK"] - 1}] + 1]]];

End[];
