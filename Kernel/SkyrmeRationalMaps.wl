(* ::Package:: *)

(* ::Title:: *)
(*CertiField: certified upper bounds for the Skyrme energy in degree B (rational-map ansatz)*)

(* ::Text:: *)
(*Statement and proofs: docs/SKYRME_RATIONAL_MAPS.md.*)
(**)
(*The competitor is U = Exp[I f(r) n_R . sigma] with R = p/q a rational map of degree B and f a radial profile.*)
(*Its energy, in the units of the paclet, is*)
(*    E = 4 Pi Integrate[ r^2 f'^2 + 2 B (f'^2 + 1) Sin[f]^2 + II Sin[f]^4/r^2 + 2 mu^2 r^2 (1 - Cos[f]), {r, 0, Infinity} ],*)
(*    II = (1/(4 Pi)) Integrate[ J^4 dA over the sphere ],   J = (1 + |z|^2) |R'(z)| / (1 + |R(z)|^2).*)
(*Rigorous in this file:*)
(*  rationalMapIntegral   an enclosure of II: midpoint rule on dyadic cells with the remainder bounded through*)
(*                        enclosures of the second derivatives (interval arithmetic on second-order Taylor data,*)
(*                        outward-rounded machine arithmetic of RigorousArithmetic.wl), exact sums;*)
(*  rationalMapUpperBound the radial integral for a piecewise cubic profile with rational data, as in*)
(*                        SkyrmeHedgehog.wl, with B and an upper bound of II in the two places where the hedgehog has 1.*)
(*Heuristic: the choice of the map (rationalMapOfDegree) and of the profile (rationalMapProfile, rationalMapTrialData).*)

Begin["CertiField`Private`"];

(* ------------------------------------------------------------------------------------------------ *)
(* interval arithmetic on second-order Taylor data in two variables                                  *)
(* ------------------------------------------------------------------------------------------------ *)
(* An interval is {lo, hi} with lo, hi packed arrays of the same length.  A jet is a list of six intervals
   {f, f_1, f_2, f_11, f_12, f_22}: enclosures of a function and of its partial derivatives over a box.  The functions
   below act on intervals (Length 2) and on jets (Length 6) alike, so the same code evaluates a function at the
   centres of the cells and encloses its second derivatives over the cells. *)

rmNegIv[{lo_, hi_}] := {-hi, -lo};
rmTwice[{lo_, hi_}] := {2. lo, 2. hi};                      (* exact *)
rmJetQ[a_] := Length[a] == 6;

rmAdd[a_, b_] := If[rmJetQ[a], MapThread[ivAdd, {a, b}], ivAdd[a, b]];
rmSub[a_, b_] := If[rmJetQ[a], MapThread[ivSub, {a, b}], ivSub[a, b]];
rmNeg[a_] := If[rmJetQ[a], rmNegIv /@ a, rmNegIv[a]];
(* multiplication by a scalar interval c = {clo, chi} *)
rmScale[c_, a_] := If[rmJetQ[a], ivMul[c, #] & /@ a, ivMul[c, a]];
(* the product rule up to second order *)
rmMul[a_, b_] := If[! rmJetQ[a], ivMul[a, b],
  {ivMul[a[[1]], b[[1]]],
   ivAdd[ivMul[a[[2]], b[[1]]], ivMul[a[[1]], b[[2]]]],
   ivAdd[ivMul[a[[3]], b[[1]]], ivMul[a[[1]], b[[3]]]],
   ivAdd[ivAdd[ivMul[a[[4]], b[[1]]], ivMul[a[[1]], b[[4]]]], rmTwice[ivMul[a[[2]], b[[2]]]]],
   ivAdd[ivAdd[ivMul[a[[5]], b[[1]]], ivMul[a[[1]], b[[5]]]], ivAdd[ivMul[a[[2]], b[[3]]], ivMul[a[[3]], b[[2]]]]],
   ivAdd[ivAdd[ivMul[a[[6]], b[[1]]], ivMul[a[[1]], b[[6]]]], rmTwice[ivMul[a[[3]], b[[3]]]]]}];

(* a constant with the shape of the element `like` *)
rmConst[c_, like_] := If[rmJetQ[like],
  With[{n = Length[like[[1, 1]]]}, With[{z = ConstantArray[0., n]},
    {{ConstantArray[c[[1]], n], ConstantArray[c[[2]], n]}, {z, z}, {z, z}, {z, z}, {z, z}, {z, z}}]],
  With[{n = Length[like[[1]]]}, {ConstantArray[c[[1]], n], ConstantArray[c[[2]], n]}]];

(* 1/s for s with positive lower end.  Returns {1/s, ok}: ok is the 0/1 mask of the entries where the enclosure of s
   is positive; elsewhere s is replaced by 1 and the caller must discard the entry. *)
rmRecip[s_] := Module[{v = If[rmJetQ[s], s[[1]], s], ok, safe, r, rr, rrr},
  ok = gtMask[v[[1]], 0.];
  safe = {vWhere[ok, v[[1]], 1.], vWhere[ok, v[[2]], 1.]};
  r = ivDivPos[{1., 1.}, safe];
  If[! rmJetQ[s], Return[{r, ok}, Module]];
  rr = ivMul[r, r]; rrr = ivMul[rr, r];
  {{r,
    rmNegIv[ivMul[s[[2]], rr]],
    rmNegIv[ivMul[s[[3]], rr]],
    ivSub[rmTwice[ivMul[ivMul[s[[2]], s[[2]]], rrr]], ivMul[s[[4]], rr]],
    ivSub[rmTwice[ivMul[ivMul[s[[2]], s[[3]]], rrr]], ivMul[s[[5]], rr]],
    ivSub[rmTwice[ivMul[ivMul[s[[3]], s[[3]]], rrr]], ivMul[s[[6]], rr]]}, ok}];

(* complex numbers are pairs {re, im} of elements; a complex scalar coefficient is {cr, ci} with cr, ci scalar
   intervals or the exact symbol 0 *)
rmCAdd[{a_, b_}, {c_, d_}] := {rmAdd[a, c], rmAdd[b, d]};
rmCMul[{a_, b_}, {c_, d_}] := {rmSub[rmMul[a, c], rmMul[b, d]], rmAdd[rmMul[a, d], rmMul[b, c]]};
rmCScale[{cr_, ci_}, {a_, b_}] := Which[
  ci === 0, {rmScale[cr, a], rmScale[cr, b]},
  cr === 0, {rmNeg[rmScale[ci, b]], rmScale[ci, a]},
  True, {rmSub[rmScale[cr, a], rmScale[ci, b]], rmAdd[rmScale[cr, b], rmScale[ci, a]]}];
rmCConst[{cr_, ci_}, like_] := {rmConst[If[cr === 0, {0., 0.}, cr], like], rmConst[If[ci === 0, {0., 0.}, ci], like]};

(* a Gaussian rational as a complex scalar coefficient, None if it is zero *)
rmCoefIv[c_] := If[c == 0, None, {If[Re[c] == 0, 0, {fLo[Re[c]], fHi[Re[c]]}], If[Im[c] == 0, 0, {fLo[Im[c]], fHi[Im[c]]}]}];

(* Sum_k coef_k a^k c^(deg - k), with cpow[[j + 1]] = c^j *)
rmHorner[coef_List, a_, cpow_List, like_] := Module[{deg = Length[coef] - 1, acc},
  acc = If[coef[[deg + 1]] === None, rmCConst[{0, 0}, like], rmCConst[coef[[deg + 1]], like]];
  Do[acc = rmCMul[acc, a];
    If[coef[[k + 1]] =!= None, acc = rmCAdd[acc, rmCScale[coef[[k + 1]], cpow[[deg - k + 1]]]]],
    {k, deg - 1, 0, -1}];
  acc];

(* G(rho, t) = 2 rho (1 + rho^2)^2 (1 + t^2)^3 |W|^4 / (|P|^2 + |Q|^2)^4 for the elements rho, t (both intervals or
   both jets), where P = Sum_k a_k A^k C^(B - k), A = rho (1 + I t), C = 1 - I t, and likewise Q and W.
   Returns {G, ok}. *)
rmEvalG[{pc_, qc_, wc_}, rho_, t_] := Module[{one, a, c, cpow, pp, qq, ww, s, m, a1, a2, pref, num, rec, ok, rec2},
  one = rmConst[{1., 1.}, rho];
  a = {rho, rmMul[rho, t]};
  c = {one, rmNeg[t]};
  cpow = FoldList[rmCMul[#1, c] &, {one, rmConst[{0., 0.}, rho]}, Range[Max[Length[pc], Length[wc]] - 1]];
  pp = rmHorner[pc, a, cpow, rho]; qq = rmHorner[qc, a, cpow, rho]; ww = rmHorner[wc, a, cpow, rho];
  s = rmAdd[rmAdd[rmMul[pp[[1]], pp[[1]]], rmMul[pp[[2]], pp[[2]]]], rmAdd[rmMul[qq[[1]], qq[[1]]], rmMul[qq[[2]], qq[[2]]]]];
  m = rmAdd[rmMul[ww[[1]], ww[[1]]], rmMul[ww[[2]], ww[[2]]]];
  a1 = rmAdd[one, rmMul[rho, rho]]; a2 = rmAdd[one, rmMul[t, t]];
  pref = rmMul[rmScale[{2., 2.}, rho], rmMul[rmMul[a1, a1], rmMul[a2, rmMul[a2, a2]]]];
  num = rmMul[pref, rmMul[m, m]];
  {rec, ok} = rmRecip[s];
  rec2 = rmMul[rec, rec];
  {rmMul[num, rmMul[rec2, rec2]], ok}];

(* ------------------------------------------------------------------------------------------------ *)
(* the angular integral                                                                              *)
(* ------------------------------------------------------------------------------------------------ *)

rmExactTotal[x_] := Total[SetPrecision[x, Infinity]];

(* Enclosure {lo, hi} (exact rationals) of the integral of G over [0, 1] x [-1, 1] for one piece: adaptive midpoint
   rule.  On a cell with centre m and sides h1, h2 the integral lies in
       h1 h2 (G(m) + h1^2/24 [G_11] + h2^2/24 [G_22]),     [.] the enclosure over the cell
   (docs/SKYRME_RATIONAL_MAPS.md, Lemma 3).  A cell is accepted when the width of its enclosure is below its share of
   the tolerance, and bisected in both directions otherwise.  Corners and centres are dyadic, hence exact. *)
rmPieceBound[piece_, tol_, maxLevel_] := Module[
  {n1 = 8, n2 = 16, rlo, rhi, tlo, thi, accLo = 0, accHi = 0, est = None, cells = 0, n, h1, h2, rm, tm, g0, ok0, gj, ok2,
   one, zero, c1, c2, cell, ok, w, acc, rej, level = 0},
  {rlo, tlo} = Transpose[Tuples[{N[Range[0, n1 - 1]/n1], N[-1 + 2 Range[0, n2 - 1]/n2]}]];
  rhi = rlo + 1./n1; thi = tlo + 2./n2;
  While[Length[rlo] > 0,
    n = Length[rlo]; cells += n;
    h1 = rhi - rlo; h2 = thi - tlo; rm = 0.5 (rlo + rhi); tm = 0.5 (tlo + thi);
    {g0, ok0} = rmEvalG[piece, {rm, rm}, {tm, tm}];
    one = ConstantArray[1., n]; zero = ConstantArray[0., n];
    {gj, ok2} = rmEvalG[piece,
      {{rlo, rhi}, {one, one}, {zero, zero}, {zero, zero}, {zero, zero}, {zero, zero}},
      {{tlo, thi}, {zero, zero}, {one, one}, {zero, zero}, {zero, zero}, {zero, zero}}];
    c1 = ivDivPos[{h1 h1, h1 h1}, {24., 24.}]; c2 = ivDivPos[{h2 h2, h2 h2}, {24., 24.}];
    cell = ivMul[{h1 h2, h1 h2}, ivAdd[g0, ivAdd[ivMul[c1, gj[[4]]], ivMul[c2, gj[[6]]]]]];
    ok = ok0 ok2;
    If[! (VectorQ[cell[[1]], Developer`MachineRealQ] && VectorQ[cell[[2]], Developer`MachineRealQ] && Max[Abs[cell[[2]]]] < $big),
      Throw[failure["Inconclusive", "The enclosure of the angular integral overflowed."], "CertiField"]];
    If[est === None, est = Max[Total[ok0 Abs[ivMul[{h1 h2, h1 h2}, g0][[2]]]], $MinMachineNumber]];
    w = cell[[2]] - cell[[1]];
    acc = Unitize[ok If[level >= maxLevel, one, leMask[w, 0.5 tol est h1 h2]]];
    If[level >= maxLevel && Total[ok] < n,
      Throw[failure["Inconclusive", "The denominator of the angular integrand could not be shown to be positive."], "CertiField"]];
    accLo += rmExactTotal[Pick[cell[[1]], acc, 1]]; accHi += rmExactTotal[Pick[cell[[2]], acc, 1]];
    rej = Flatten[Position[acc, 0]];
    If[Length[rej] == 0, Break[]];
    If[cells + 4 Length[rej] > 4000000,
      Throw[failure["Inconclusive", "The angular integral needs more than four million cells."], "CertiField"]];
    {rlo, rhi, tlo, thi, rm, tm} = #[[rej]] & /@ {rlo, rhi, tlo, thi, rm, tm};
    {rlo, rhi, tlo, thi} = {Join[rlo, rm, rlo, rm], Join[rm, rhi, rm, rhi], Join[tlo, tlo, tm, tm], Join[tm, tm, thi, thi]};
    level++];
  {accLo, accHi, cells}];

(* p, q: polynomials in z with Gaussian rational coefficients, without a common zero, of degree b as a map.
   The four pieces of the sphere: (p, q), z -> -z, and the same for the map composed with the inversion. *)
rationalMapQ[p_, q_, z_Symbol, b_Integer] := b >= 1 && PolynomialQ[p, z] && PolynomialQ[q, z] &&
  AllTrue[Join[CoefficientList[p, z], CoefficientList[q, z]], MatchQ[#, _Integer | _Rational | Complex[_Integer | _Rational, _Integer | _Rational]] &] &&
  Max[Exponent[p, z], Exponent[q, z]] == b && ! PossibleZeroQ[q] && Resultant[p, q, z] =!= 0;

rmPieces[p_, q_, z_, b_] := Module[{a = PadRight[CoefficientList[p, z], b + 1], c = PadRight[CoefficientList[q, z], b + 1],
   sgn = (-1)^Range[0, b], mk},
  mk[aa_, cc_] := Module[{pp = aa . z^Range[0, b], qq = cc . z^Range[0, b]},
    {rmCoefIv /@ aa, rmCoefIv /@ cc, rmCoefIv /@ PadRight[CoefficientList[Expand[D[pp, z] qq - pp D[qq, z]], z], Max[2 b - 1, 1]]}];
  {mk[a, c], mk[sgn a, sgn c], mk[Reverse[a], Reverse[c]], mk[sgn Reverse[a], sgn Reverse[c]]}];

Options[rationalMapIntegral] = {"Tolerance" -> 10^-5, "MaxLevel" -> 9};
(* enclosure of the angular integral: <|"Lower" -> rational, "Upper" -> rational, "Cells" -> number of cells|> *)
rationalMapIntegral[p_, q_, z_Symbol, b_Integer, OptionsPattern[]] := Module[{res, lo, hi},
  If[! rationalMapQ[p, q, z, b],
    Throw[failure["BadModel", "The rational map must be given by two polynomials with Gaussian rational coefficients and no common zero, of the stated degree."], "CertiField"]];
  res = rmPieceBound[#, N[OptionValue["Tolerance"]], OptionValue["MaxLevel"]] & /@ rmPieces[p, q, z, b];
  lo = Total[res[[All, 1]]]; hi = Total[res[[All, 2]]];
  <|"Lower" -> If[lo >= 0, lo/$piHi, lo/$piLo], "Upper" -> If[hi >= 0, hi/$piLo, hi/$piHi], "Cells" -> Total[res[[All, 3]]]|>];

(* ------------------------------------------------------------------------------------------------ *)
(* the radial integral (after hedgehogUpperBound)                                                    *)
(* ------------------------------------------------------------------------------------------------ *)

rationalMapDensityP[p_, dp_, r_, mu2_, b_, ii_] := r^2 dp^2 + 2 b Sin[p]^2 (1 + dp^2) + ii Sin[p]^4/r^2 + 2 mu2 r^2 (1 + Cos[p]);

(* Rigorous upper bound, in units of 12 Pi^2, for the energy of the competitor with degree b, a rational number
   ii >= II (an upper bound of the angular integral) and the profile f = Pi - p of the trial data. *)
Options[rationalMapUpperBound] = {WorkingPrecision -> 30, "Tolerance" -> 10^-10};
rationalMapUpperBound[trial_Association, b_Integer, ii_, muq_, OptionsPattern[]] := Module[
  {r = trial["r"], p = trial["p"], dp = trial["dp"], n, x, mu2 = muq^2, poly, q, maj, first, middle, last, fmax, df,
   prec = OptionValue[WorkingPrecision], tol = OptionValue["Tolerance"], total, rn1, rr, t0 = AbsoluteTime[], piIv},
  n = Length[r] - 1;
  If[! (b >= 1 && MatchQ[ii, _Integer | _Rational] && ii >= 0 && MatchQ[muq, _Integer | _Rational] && muq >= 0 && n >= 2 &&
      r[[1]] === 0 && p[[1]] === 0 && OrderedQ[r] && Length[p] == n && Length[dp] == n &&
      AllTrue[Join[r, p, dp], MatchQ[#, _Integer | _Rational] &]),
    Throw[failure["BadCandidate", "The trial profile must have exact rational data with p(0) = 0."], "CertiField"]];
  (* first cell: polynomial majorant (|sin p| <= |p|, 1 + cos p <= 2), integrated exactly *)
  poly = hermiteCubic[r[[{1, 2}]], p[[{1, 2}]], dp[[{1, 2}]], x];
  q = Cancel[poly/x];
  If[! PolynomialQ[q, x], Throw[failure["Internal", "p(0) != 0"], "CertiField"]];
  maj = x^2 D[poly, x]^2 + 2 b poly^2 (1 + D[poly, x]^2) + ii q^4 x^2 + 4 mu2 x^2;
  first = Integrate[maj, {x, 0, r[[2]]}];
  (* middle cells *)
  middle = Total[Table[
    poly = hermiteCubic[r[[{k, k + 1}]], p[[{k, k + 1}]], dp[[{k, k + 1}]], x];
    rigorousIntegral[rationalMapDensityP[poly, D[poly, x], x, mu2, b, ii], {x, r[[k]], r[[k + 1]]},
      WorkingPrecision -> prec, "Cells" -> 1, "Tolerance" -> tol (r[[k + 1]] - r[[k]])/Last[r]],
    {k, 2, n - 1}]];
  (* last cell: f = (Pi - p_{n-1}) (R - x)/(R - r_{n-1}); |sin f| <= |f|, 1 - cos f <= f^2/2, 1/x^2 <= 1/r_{n-1}^2 *)
  rn1 = r[[n]]; rr = r[[n + 1]];
  fmax = Max[Abs[$piHi - p[[n]]], Abs[$piLo - p[[n]]]];
  df = fmax/(rr - rn1);
  maj = With[{fx = fmax (rr - x)/(rr - rn1)}, x^2 df^2 + 2 b fx^2 (1 + df^2) + ii fx^4/rn1^2 + mu2 x^2 fx^2];
  last = Integrate[maj, {x, rn1, rr}];
  total = toIv[first + last, prec] + middle;
  piIv = Interval[{$piLo, $piHi}];
  <|"UpperBound" -> Max[total/(3 piIv)], "Enclosure" -> total/(3 piIv), "Cells" -> n,
    "Seconds" -> AbsoluteTime[] - t0, "Radius" -> rr|>];

(* ------------------------------------------------------------------------------------------------ *)
(* heuristics: the maps of the literature and the profile                                            *)
(* ------------------------------------------------------------------------------------------------ *)

(* Houghton, Manton and Sutcliffe, Nucl. Phys. B 510 (1998) 507, Section 3: equations (3.6) with a = I, (3.9) with
   c = 1, (3.11) with a = 3.07, b = 3.94, (3.15) with a = 0.16, (3.17) with a = -1/7 and (3.16) with a = 0.136.
   s3 is a rational number close to Sqrt[3]: the maps of degree 3 and 4 are then not exactly the symmetric ones, which
   changes only the value of the integral.  In (3.17) the paper prints a = -3; that value gives the same Wronskian
   up to a factor but an integral of 216.5, and a = -1/7 gives the icosahedral map with 60.87.  In (3.16) the paper
   gives a = 0.14; 0.136 gives a slightly smaller integral. *)
rationalMapOfDegree[b_Integer, z_] := With[{s3 = 1732051/1000000},
  Switch[b,
    1, {z, 1},
    2, {z^2, 1},
    3, {s3 I z^2 - 1, z (z^2 - s3 I)},
    4, {z^4 + 2 s3 I z^2 + 1, z^4 - 2 s3 I z^2 + 1},
    5, {z (z^4 + 197/50 z^2 + 307/100), 307/100 z^4 - 197/50 z^2 + 1},
    6, {z^4 + 4/25 I, z^2 (4/25 I z^4 + 1)},
    7, {7 z^5 + 1, z^2 (7 - z^5)},
    8, {z^6 - 17/125, z^2 (17/125 z^6 + 1)},
    _, $Failed]];

(* nodes of the trial profile: steps 1/25 up to the radius core, then growing by 4 % up to 1/2 *)
rationalMapNodes[core_, rmax_] := Module[{nodes = {0}, h = 1/25},
  While[Last[nodes] < rmax, If[Last[nodes] >= core, h = Min[h 104/100, 1/2]]; AppendTo[nodes, Min[Last[nodes] + h, rmax]]];
  nodes];

(* the minimiser of the midpoint discretisation of the radial energy on the nodes, by a damped Newton iteration
   (floating point; any values give a valid trial profile) *)
rationalMapProfile[b_, ii_, mu_, nodes_] := Module[
  {r = N[nodes], n, h, rm, f, s, r0, bb = N[b], iv = N[ii], m2 = N[mu]^2, energy, fm, d, sf, cf, s2, c2, ef, ed, eff, efd, edd,
   grad, diag, off, hess, delta, e0, e1, step, lam = 0., ok},
  n = Length[r] - 1; h = Differences[r]; rm = (Most[r] + Rest[r])/2.;
  s = (Sqrt[1. + 8. bb] - 1.)/2.; r0 = 0.8 iv^0.25;
  f = Join[{N[Pi]}, 2. ArcTan[(r0/r[[2 ;; -2]])^(s + 1.)], {0.}];
  energy[ff_] := Module[{mid = (Most[ff] + Rest[ff])/2., dd = Differences[ff]/h},
    Total[h (rm^2 dd^2 + 2. bb (dd^2 + 1.) Sin[mid]^2 + iv Sin[mid]^4/rm^2 + 2. m2 rm^2 (1. - Cos[mid]))]];
  e0 = energy[f];
  Do[
    fm = (Most[f] + Rest[f])/2.; d = Differences[f]/h;
    sf = Sin[fm]; cf = Cos[fm]; s2 = Sin[2. fm]; c2 = Cos[2. fm];
    ef = 2. bb (d^2 + 1.) s2 + 2. iv sf^2 s2/rm^2 + 2. m2 rm^2 sf;
    ed = 2. rm^2 d + 4. bb d sf^2;
    eff = 4. bb (d^2 + 1.) c2 + iv (2. s2^2 + 4. sf^2 c2)/rm^2 + 2. m2 rm^2 cf;
    efd = 4. bb d s2;
    edd = 2. rm^2 + 4. bb sf^2;
    (* interior nodes 2 .. n: the right end of cell k - 1 and the left end of cell k *)
    grad = Most[h (ef/2. + ed/h)] + Rest[h (ef/2. - ed/h)];
    diag = Most[h (eff/4. + efd/h + edd/h^2)] + Rest[h (eff/4. - efd/h + edd/h^2)];
    off = (h (eff/4. - edd/h^2))[[2 ;; -2]];
    If[Max[Abs[grad]] < 1.*^-11 (1. + Abs[e0]), Break[]];
    ok = False;
    Do[
      hess = SparseArray[{Band[{1, 1}] -> diag + lam, Band[{1, 2}] -> off, Band[{2, 1}] -> off}, {n - 1, n - 1}];
      delta = Quiet@LinearSolve[hess, -grad];
      If[VectorQ[delta, Developer`MachineRealQ] && delta . grad < 0.,
        step = 1.;
        Do[e1 = energy[Join[{N[Pi]}, f[[2 ;; -2]] + step delta, {0.}]];
          If[TrueQ[e1 < e0], ok = True; Break[]]; step = step/2., {30}]];
      If[ok, Break[]];
      lam = If[lam == 0., 1.*^-3 Max[Abs[diag]], 10. lam], {12}];
    If[! ok, Break[]];
    f = Join[{N[Pi]}, f[[2 ;; -2]] + step delta, {0.}]; e0 = e1; lam = lam/10.,
    {80}];
  f];

(* rational trial data {r, p, dp} for rationalMapUpperBound: p = Pi - f at the nodes, slopes from the three-point
   difference formula *)
rationalMapTrialData[b_Integer, ii_, mu_] := Module[{core, rmax, nodes, f, p, n, h, dp, rat},
  core = Max[3, Ceiling[2 N[ii]^(1/4)]];
  rmax = If[mu == 0, 60, Min[60, Max[Ceiling[15/mu], core + 5]]];
  nodes = rationalMapNodes[core, rmax];
  f = rationalMapProfile[b, ii, mu, nodes];
  p = N[Pi] - f; n = Length[nodes] - 1; h = Differences[N[nodes]];
  dp = Join[{Max[0., (4. p[[2]] - p[[3]] - 3. p[[1]])/(2. h[[1]])]},
    Table[(h[[k - 1]]/(h[[k]] (h[[k - 1]] + h[[k]]))) (p[[k + 1]] - p[[k]]) + (h[[k]]/(h[[k - 1]] (h[[k - 1]] + h[[k]]))) (p[[k]] - p[[k - 1]]), {k, 2, n}]];
  rat[v_] := Rationalize[SetPrecision[v, 22], 10^-18];
  <|"Degree" -> b, "Mu" -> mu, "r" -> nodes, "p" -> Join[{0}, rat /@ p[[2 ;; n]]], "dp" -> rat /@ dp|>];

(* ------------------------------------------------------------------------------------------------ *)
(* the model SkyrmeRationalMap                                                                       *)
(* ------------------------------------------------------------------------------------------------ *)

CertiFieldModelInfo["SkyrmeRationalMap"] := <|
  "Name" -> "SkyrmeRationalMap",
  "Energy" -> "E(U) = \[Integral] |dU|^2 + |\[CapitalLambda]^2 dU|^2 + 2\[Mu]^2 (1 - cos f) over \[DoubleStruckCapitalR]^3, f = dist(U, 1); reported as E/(12\[Pi]^2)",
  "Parameters" -> <|"Degree" -> "B, a positive integer (default 2); maps are built in for B = 1, ..., 8",
    "PionMass" -> "a rational \[Mu] with 0 \[LessEqual] \[Mu] < 2 (default 0)",
    "RationalMap" -> "optional: a pure function of z giving a rational map of degree B with Gaussian rational coefficients; it replaces the built-in map"|>,
  "AdmissibleClass" -> "finite-energy maps U: \[DoubleStruckCapitalR]^3 \[RightArrow] S^3 of degree B with U - 1 \[Element] L^6",
  "LowerBound" -> "the topological bound E \[GreaterEqual] 12\[Pi]^2 |B| (not tight)",
  "UpperBound" -> "energy of the field of the rational-map ansatz, U = exp(i f(r) n_R\[CenterDot]\[Sigma]): certified enclosure of the angular integral of the map and rigorous quadrature of the radial integral",
  "Reference" -> "docs/SKYRME_RATIONAL_MAPS.md",
  "Status" -> "Supported"|>;

(* {B, mu, p, q}: p, q polynomials in \[FormalZ] *)
rationalMapSpec[params_Association] := Module[{b = Lookup[params, "Degree", 2], mu = Lookup[params, "PionMass", 0],
   map = Lookup[params, "RationalMap", Automatic], z = \[FormalZ], r, pq},
  mu = If[InexactNumberQ[mu], Rationalize[mu, 10^-12], mu];
  If[! (MatchQ[mu, _Integer | _Rational] && 0 <= mu < 2),
    Throw[failure["BadModel", "The pion mass must be a rational number with 0 \[LessEqual] \[Mu] < 2."], "CertiField"]];
  If[! (IntegerQ[b] && b >= 1), Throw[failure["BadModel", "\"Degree\" must be a positive integer."], "CertiField"]];
  pq = If[map === Automatic,
    If[b > 8, Throw[failure["BadModel", "Maps are built in for the degrees 1 to 8; give \"RationalMap\" for another degree."], "CertiField"]];
    rationalMapOfDegree[b, z],
    r = Quiet@Together[map[z]];
    {Numerator[r], Denominator[r]}];
  If[! (ListQ[pq] && Length[pq] == 2 && TrueQ[rationalMapQ[pq[[1]], pq[[2]], z, b]]),
    Throw[failure["BadModel", "\"RationalMap\" must be a pure function of z giving a quotient of two polynomials with Gaussian rational coefficients and no common zero, of degree \"Degree\" as a map."], "CertiField"]];
  {b, mu, pq[[1]], pq[[2]]}];

rationalMapIntegralCached[p_, q_, z_, b_] := rationalMapIntegralCached[p, q, z, b] = rationalMapIntegral[p, q, z, b];

(* The stored upper bound is the computed one weakened by one part in 10^9. *)
certifySkyrmeRationalMap[params_Association, cand_, opts_List] := Module[
  {spec = rationalMapSpec[params], verbose = TrueQ[OptionValue[CertifiedVariationalBound, opts, "Verbose"]], z = \[FormalZ],
   b, mu, p, q, int, ii, trial, ub, t0 = AbsoluteTime[]},
  {b, mu, p, q} = spec;
  If[verbose, Print["  enclosing the angular integral of the map"]];
  int = rationalMapIntegralCached[p, q, z, b];
  ii = upperDecimal[int["Upper"], 12];
  If[verbose, Print["  angular integral in [", N[int["Lower"], 12], ", ", N[int["Upper"], 12], "], ", int["Cells"], " cells"]];
  trial = Which[
    cand === Automatic, rationalMapTrialData[b, ii, mu],
    AssociationQ[cand] && KeyExistsQ[cand, "r"], cand,
    True, Throw[failure["BadCandidate", "The candidate must be trial data <|\"r\" -> ..., \"p\" -> ..., \"dp\" -> ...|> (rational numbers, p = \[Pi] - f)."], "CertiField"]];
  trial = KeyTake[trial, {"r", "p", "dp"}];
  ub = rationalMapUpperBound[trial, b, ii, mu];
  <|"LowerBound" -> b, "UpperBound" -> SetPrecision[ub["UpperBound"], Infinity] (1 + 1/10^9), "Seconds" -> AbsoluteTime[] - t0,
    "Degree" -> b, "Mu" -> mu, "Map" -> {CoefficientList[p, z], CoefficientList[q, z]},
    "AngularIntegral" -> {lowerDecimal[int["Lower"], 12], ii}, "AngularCells" -> int["Cells"], "Trial" -> trial,
    "Radial" -> KeyDrop[ub, "Enclosure"], "ComputedUpperBound" -> ub["UpperBound"]|>];

skyrmeRationalMapDescription[res_] := With[{b = res["Degree"]},
  <|"Units" -> "12\[Pi]^2 (Faddeev\[Dash]Bogomolny units), Battye\[Dash]Sutcliffe normalisation",
    "Theorem" -> Row[{"Every finite-energy map U: \[DoubleStruckCapitalR]^3 \[RightArrow] S^3 of degree ", b, " with pion mass \[Mu] = ", exactString[res["Mu"]],
      " has E(U)/(12\[Pi]^2) \[GreaterEqual] ", b, " (the topological bound).  The field of the rational-map ansatz with the stored map and profile has degree ", b,
      " and E/(12\[Pi]^2) \[LessEqual] ", upperString[res["UpperBound"], 8], ", that is ", upperString[res["UpperBound"]/b, 8], " per unit of degree."}],
    "AdmissibleClass" -> "all finite-energy maps \[DoubleStruckCapitalR]^3 \[RightArrow] S^3 of the given degree (no symmetry assumed)",
    "LowerBoundMethod" -> "the topological bound E \[GreaterEqual] 12\[Pi]^2 |B| (inequality of the arithmetic and geometric means); it is not attained, so the interval is wide",
    "UpperBoundMethod" -> Row[{"rational-map ansatz: the angular integral of the map lies in [", lowerString[res["AngularIntegral"][[1]], 10], ", ",
      upperString[res["AngularIntegral"][[2]], 10], "] (midpoint rule with enclosures of the second derivatives, ", res["AngularCells"],
      " cells); radial integral for a profile with ", res["Radial", "Cells"], " cubic cells, rigorous quadrature"}],
    "TrustedBase" -> "IEEE-754 machine arithmetic (+, \[Minus], \[Times], \[Divide]) with outward rounding, exact rational arithmetic, Wolfram Language Interval arithmetic, and docs/SKYRME_RATIONAL_MAPS.md",
    "NotClaimed" -> "that the field of the ansatz is a minimiser or a critical point of the energy, or that the infimum is attained; the lower end of the interval is only the topological bound",
    "References" -> "Houghton, Manton, Sutcliffe, Nucl. Phys. B 510 (1998) 507 (the ansatz and the maps); docs/SKYRME_RATIONAL_MAPS.md"|>];

(* recompute both integrals from the stored map and profile; the model's parameters must describe the stored map *)
verifySkyrmeRationalMap[a_Association] := Module[{c = a["Certificate"], z = \[FormalZ], spec, b, mu, p, q, int, ii, ub},
  spec = rationalMapSpec[a["Parameters"]];
  {b, mu, p, q} = spec;
  If[! (c["Degree"] === b && c["Mu"] === mu && c["Map"] === {CoefficientList[p, z], CoefficientList[q, z]}), Return[False, Module]];
  int = rationalMapIntegral[p, q, z, b];
  ii = upperDecimal[int["Upper"], 12];
  ub = rationalMapUpperBound[c["Trial"], b, ii, mu]["UpperBound"];
  NumericQ[ub] && a["LowerBound"] <= b && ub <= a["UpperBound"]];

(* ---- certificates on disk (JSON with exact numbers as strings), for the cross-check in Reference/SkyrmeRationalMaps ---- *)
rmString[x_] := ToString[x, InputForm];
rmParse[s_String] := If[StringMatchQ[s, RegularExpression["-?[0-9]+(/[0-9]+)?"]], ToExpression[s],
  Throw[failure["BadCertificate", "Not a rational number: " <> s], "CertiField"]];

exportRationalMapCertificate[file_String, cert_CertifiedBound] := Module[{c = cert["Certificate"], pair},
  pair[v_] := {rmString[Re[v]], rmString[Im[v]]};
  Export[file, <|"format" -> "CertiField-rational-map-certificate-1", "degree" -> c["Degree"], "mu" -> rmString[c["Mu"]],
    "map" -> <|"p" -> (pair /@ c["Map"][[1]]), "q" -> (pair /@ c["Map"][[2]])|>,
    "angular_integral" -> <|"lower" -> rmString[c["AngularIntegral"][[1]]], "upper" -> rmString[c["AngularIntegral"][[2]]]|>,
    "trial" -> (rmString /@ # & /@ c["Trial"]),
    "upper_bound" -> rmString[upperDecimal[cert["UpperBound"], 12]],
    "upper_bound_per_degree" -> ToString[N[upperDecimal[cert["UpperBound"], 12]/c["Degree"], 8]]|>, "RawJSON"]];

(* True if the bounds recomputed from the stored map and profile are at least as strong as the stored claim *)
verifyRationalMapCertificateFile[file_String] := Catch[Module[{j = Import[file, "RawJSON"], z = \[FormalZ], b, mu, p, q, trial, int, ii, ub, cpoly},
  cpoly[list_] := (rmParse[#[[1]]] + I rmParse[#[[2]]] & /@ list) . z^Range[0, Length[list] - 1];
  b = j["degree"]; mu = rmParse[j["mu"]];
  p = cpoly[j["map", "p"]]; q = cpoly[j["map", "q"]];
  If[! (IntegerQ[b] && b >= 1 && mu >= 0 && rationalMapQ[p, q, z, b]), Throw[False, "CertiField"]];
  trial = Association[# -> (rmParse /@ j["trial", #]) & /@ {"r", "p", "dp"}];
  int = rationalMapIntegral[p, q, z, b];
  ii = upperDecimal[int["Upper"], 12];
  ub = rationalMapUpperBound[trial, b, ii, mu]["UpperBound"];
  TrueQ[NumericQ[ub] && ub <= rmParse[j["upper_bound"]] && ii <= rmParse[j["angular_integral", "upper"]]]], "CertiField"] === True;

End[];
