(* ::Package:: *)

(* ::Title:: *)
(*CertiField: scalar fields on boxes, extensions (docs/SCALAR_ND_2.md)*)

(* ::Text:: *)
(*E[u] = Integrate[k |grad u|^2 + W(x, u)] over a box; u = g on the Dirichlet faces, free on the natural faces,*)
(*periodic in the periodic directions, optionally with a prescribed mean.*)
(*    W(x, t) = P(t) [+ c1 (t - lo) Log(t - lo) + c2 (hi - t) Log(hi - t)] [+ D_c(t)] - f(x) t ,*)
(*W = +Infinity outside the value range [lo, hi] when one is given; D_c is a table per cell, piecewise linear in t;*)
(*f is the piecewise multilinear interpolant of nodal values, or constant on each cell.*)
(*Reached through the model ScalarFieldND whenever one of the keys in $nd2Keys is present.*)
(**)
(*Rigorous: nd2UpperBound, nd2LowerBound, nd2PoissonLower, nd2LnPt, nd2XLnX, nd2LnQ, nd2SqrtUp, nd2PiLower.  The*)
(*candidate solver, the equilibrated flux of nd2MixedField and the optimiser are heuristic: they affect tightness only.*)

Begin["CertiField`Private`"];

$nd2Keys = {"Source", "Mean", "PeriodicDirections", "ValueRange", "DataTerm"};
nd2ParamsQ[params_Association] := AnyTrue[$nd2Keys, KeyExistsQ[params, #] &];
nd2ParamsQ[_] := False;

slab[n_, i_, s_] := Sequence @@ ReplacePart[ConstantArray[All, n], i -> s];

(* ------------------------------------------------------------------------------------------------ *)
(* logarithms, square roots and pi from exact rationals and the four operations (Lemma 12)           *)
(* ------------------------------------------------------------------------------------------------ *)

nd2Ln2Q[] := nd2Ln2Q[] = Module[{z = 1/3, m = 40, s},
  s = 2 Sum[z^(2 j + 1)/(2 j + 1), {j, 0, m}];
  {s, s + 2 z^(2 m + 3)/((2 m + 3) (1 - z^2))}];

(* exact rational bounds {lo, hi} for Log[q], q > 0 rational: q = r 2^e with r in [2/3, 4/3], 2 ArcTanh series *)
nd2LnQ[q_] := nd2LnQ[q] = Module[{e, r, z, s, rem, l2 = nd2Ln2Q[], m = 24},
  e = Floor[Log[2, q]]; r = q/2^e;
  If[r > 4/3, r = r/2; e = e + 1];
  z = (r - 1)/(r + 1);
  s = 2 Sum[z^(2 j + 1)/(2 j + 1), {j, 0, m}];
  rem = 2 Abs[z]^(2 m + 3)/((2 m + 3) (1 - z^2));
  {s - rem + If[e >= 0, e l2[[1]], e l2[[2]]], s + rem + If[e >= 0, e l2[[2]], e l2[[1]]]}];

nd2XLnXQ[s_] := If[s == 0, {0, 0}, s nd2LnQ[s]];

$nd2LnRem := $nd2LnRem = fHi[2 (9/50)^27/(27 (1 - (9/50)^2))];
$nd2Ln2Iv := $nd2Ln2Iv = {fLo[nd2Ln2Q[][[1]]], fHi[nd2Ln2Q[][[2]]]};

(* RIGOROUS enclosure {lo, hi} of Log at a packed array of positive machine numbers.  x = mu 2^e exactly, with mu in
   [0.707, 1.415] (the exponent estimate of Log2 is corrected by exact comparisons and checked), Log[x] = e Log[2] +
   2 ArcTanh[z], z = (mu - 1)/(mu + 1) with |z| < 9/50; 13 terms of the series and the bound for the remainder.
   Below 10^-300 the enclosure [-745.2, -690.7] is used. *)
nd2LnPt[x0_] := ndQuiet@Module[{x, tiny, e, mu, big, small, z, y, acc, s, r = $nd2LnRem},
  tiny = ltMask[x0, 1.*^-300];
  x = vWhere[tiny, 1., x0];
  e = Floor[Log2[x]];
  mu = x Developer`ToPackedArray[N[2^-e]];
  big = geMask[mu, 2.]; mu = vWhere[big, mu/2., mu]; e = e + big;
  small = ltMask[mu, 1.]; mu = vWhere[small, 2. mu, mu]; e = e - small;
  big = gtMask[mu, 1.4142135623730951]; mu = vWhere[big, mu/2., mu]; e = e + big;
  If[! (Min[mu] >= 0.7 && Max[mu] <= 1.42), Throw[failure["Internal", "Logarithm: argument reduction."], "CertiField"]];
  z = ndDiv[ivPt[mu - 1.], ndAdd[ivPt[mu], {1., 1.}]];
  y = ivSqr[z];
  acc = ivC[1/25];
  Do[acc = ndAdd[ndMul[acc, y], ivC[1/(2 j + 1)]], {j, 11, 0, -1}];
  s = ivK[ndMul[z, acc], 2];
  s = {ndDn[s[[1]] - r], ndUp[s[[2]] + r]};
  s = ndAdd[ndMul[ivPt[N[e]], $nd2Ln2Iv], s];
  {vWhere[tiny, -745.2, s[[1]]], vWhere[tiny, -690.7, s[[2]]]}];

(* RIGOROUS enclosure of s Log[s] over interval arrays s >= 0 (value 0 at 0).  Where the lower end is 0: on [0, b]
   the function is at most max(0, b Log b), and at least b Log b if b <= 1/4 (it decreases on [0, 1/E]), else -1. *)
nd2XLnX[{lo0_, hi0_}] := Module[{lo = Ramp[lo0], hi = Ramp[hi0], pos, hz, l1, h1, hh, full, b, zlo, zhi},
  pos = gtMask[lo, 0.]; hz = gtMask[hi, 0.];
  l1 = vWhere[pos, lo, 1.]; hh = vWhere[pos, vMax[hi, l1], 1.];
  full = ndMul[{l1, hh}, {nd2LnPt[l1][[1]], nd2LnPt[hh][[2]]}];
  h1 = vWhere[hz, hi, 1.];
  b = ndMul[ivPt[h1], nd2LnPt[h1]];
  zlo = vWhere[leMask[h1, 0.25], vMin[b[[1]], 0.], -1. + 0. h1];
  zhi = vMax[b[[2]], 0.];
  {vWhere[pos, full[[1]], vWhere[hz, zlo, 0. h1]], vWhere[pos, full[[2]], vWhere[hz, zhi, 0. h1]]}];

(* a machine number s >= Sqrt[x] for a machine number x >= 0, verified: s s rounded down is at least x *)
nd2SqrtUp[x_] := If[! (x > 0.), 0., Module[{s = rUp[Sqrt[x]]}, While[rDn[s s] < x, s = rUp[s]]; s]];

(* an exact rational lower bound for Pi: Machin's formula with alternating series cut after a negative term *)
nd2PiLower[] := nd2PiLower[] = Module[{at},
  at[q_, last_] := Sum[(-1)^j/((2 j + 1) q^(2 j + 1)), {j, 0, last}];
  16 at[5, 25] - 4 at[239, 24]];

(* machine x Log[x] with 0 at 0 (heuristics only) *)
nd2XLogMach[s_] := With[{p = Ramp[s]}, p Log[p + (1 - Unitize[p])]];

(* ------------------------------------------------------------------------------------------------ *)
(* the model                                                                                         *)
(* ------------------------------------------------------------------------------------------------ *)

nd2Spec[params_Association] := Module[
  {dom, n, per, nat, vr, V, s, l1, l2, vexpr, logs, e2, c1 = 0, c2 = 0, P, base, cells, verts, tab, tm = 0, src, lin = None,
   linNodal = False, av = None, mean, g, kind, gf},
  dom = exactify[Lookup[params, "Domain", Missing[]]];
  If[! (MatchQ[dom, {{_, _} ..}] && 2 <= Length[dom] <= 3),
    Throw[failure["BadModel", "\"Domain\" must be {{a1, b1}, {a2, b2}} or {{a1, b1}, {a2, b2}, {a3, b3}} with rational a < b."], "CertiField"]];
  n = Length[dom];
  per = Lookup[params, "PeriodicDirections", {}];
  If[! (VectorQ[per, IntegerQ] && DuplicateFreeQ[per] && AllTrue[per, 1 <= # <= n &]),
    Throw[failure["BadModel", "\"PeriodicDirections\" must be a list of distinct directions between 1 and the dimension."], "CertiField"]];
  nat = Lookup[params, "NaturalFaces", {}];
  If[! ListQ[nat], Throw[failure["BadModel", "\"NaturalFaces\" must be a list of {direction, \"Lower\" | \"Upper\"}."], "CertiField"]];
  If[AnyTrue[nat, MatchQ[#, {i_, _} /; MemberQ[per, i]] &],
    Throw[failure["BadModel", "A periodic direction has no natural faces."], "CertiField"]];
  vr = Lookup[params, "ValueRange", None];
  If[vr =!= None,
    vr = exactify[vr];
    If[! (MatchQ[vr, {_Integer | _Rational, _Integer | _Rational}] && vr[[1]] < vr[[2]]),
      Throw[failure["BadModel", "\"ValueRange\" must be {lo, hi} with rational lo < hi."], "CertiField"]]];
  V = Lookup[params, "Potential", Function[u, 0]];
  If[! MatchQ[V, _Function | _Symbol], Throw[failure["BadModel", "\"Potential\" must be a pure function of one variable."], "CertiField"]];
  vexpr = Expand[exactify[V[s]]];
  If[PolynomialQ[vexpr, s],
    P = vexpr,
    If[vr === None, Throw[failure["BadModel", "A potential with logarithms needs \"ValueRange\" -> {lo, hi}."], "CertiField"]];
    logs = Union[Cases[vexpr, Log[_], Infinity]];
    e2 = Expand[vexpr /. Join[
      Cases[logs, Log[a_] /; Expand[a - (s - vr[[1]])] === 0 :> (Log[a] -> l1)],
      Cases[logs, Log[a_] /; Expand[a - (vr[[2]] - s)] === 0 :> (Log[a] -> l2)]]];
    c1 = Coefficient[Coefficient[e2, l1], s, 1]; c2 = -Coefficient[Coefficient[e2, l2], s, 1];
    P = Expand[e2 /. {l1 -> 0, l2 -> 0}];
    If[! (PolynomialQ[e2, {s, l1, l2}] && Exponent[e2, l1] <= 1 && Exponent[e2, l2] <= 1 &&
        MatchQ[c1, _Integer | _Rational] && MatchQ[c2, _Integer | _Rational] && c1 >= 0 && c2 >= 0 &&
        Expand[Coefficient[e2, l1] - c1 (s - vr[[1]])] === 0 && Expand[Coefficient[e2, l2] - c2 (vr[[2]] - s)] === 0),
      Throw[failure["BadModel", "With \"ValueRange\" -> {lo, hi} the potential must be a polynomial plus c1 (u - lo) Log[u - lo] + c2 (hi - u) Log[hi - u] with rational c1, c2 >= 0."], "CertiField"]]];
  P = CoefficientList[P, s];
  If[P === {}, P = {0}];
  If[! VectorQ[P, MatchQ[#, _Integer | _Rational] &],
    Throw[failure["BadModel", "The potential must have rational coefficients and no free parameters."], "CertiField"]];
  (* grid, Dirichlet faces and data from the base model; periodic faces are passed to it as natural faces *)
  base = scalarFieldNDSpec[Join[KeyDrop[params, Join[$nd2Keys, {"Field"}]],
    <|"Potential" -> Function @@ {\[FormalU], P . \[FormalU]^Range[0, Length[P] - 1]},
      "NaturalFaces" -> Join[nat, Flatten[Table[{{i, "Lower"}, {i, "Upper"}}, {i, per}], 1]]|>]];
  cells = base["Cells"];
  verts = Tuples[{0, 1}, n];
  tab = Lookup[params, "DataTerm", None];
  If[tab =!= None,
    If[vr === None, Throw[failure["BadModel", "\"DataTerm\" needs \"ValueRange\" -> {lo, hi}."], "CertiField"]];
    If[c1 != 0 || c2 != 0, Throw[failure["NotImplemented", "A data term together with logarithmic terms is not supported in this version."], "CertiField"]];
    tm = If[ArrayQ[tab] && ArrayDepth[tab] == n + 1 && Most[Dimensions[tab]] === cells, Last[Dimensions[tab]] - 1, 0];
    tab = If[tm >= 1, ndRealArray[tab, Append[cells, tm + 1]], $Failed];
    If[tab === $Failed,
      Throw[failure["BadModel", "\"DataTerm\" must be a real array of dimensions Cells x (M + 1), M >= 1: in each cell the values at lo + j (hi - lo)/M, j = 0, ..., M."], "CertiField"]];
    tab = ArrayReshape[tab, {Times @@ cells, tm + 1}]];
  src = Lookup[params, "Source", None];
  If[src =!= None,
    lin = Which[
      NumericQ[src], ndRealArray[ConstantArray[N[src], cells], cells],
      ArrayQ[src, n] && Dimensions[src] === cells + 1, linNodal = True; ndRealArray[src, cells + 1],
      ArrayQ[src, n] && Dimensions[src] === cells, ndRealArray[src, cells],
      True, linNodal = True; ndRealArray[Quiet@Apply[src, Outer[List, Sequence @@ base["Axes"]], {n}], cells + 1]];
    If[lin === $Failed,
      Throw[failure["BadModel", "\"Source\" must be a real number, a real function of the coordinates, or an array of nodal (Cells + 1) or cell (Cells) values."], "CertiField"]];
    (* the coefficient a of the term a(x) t of the potential is -f *)
    av = Association[# -> If[linNodal, -vertexPart[lin, n, #], -lin] & /@ verts]];
  mean = Lookup[params, "Mean", None];
  If[mean =!= None,
    mean = exactify[mean];
    If[! MatchQ[mean, _Integer | _Rational], Throw[failure["BadModel", "\"Mean\" must be a rational number."], "CertiField"]];
    If[vr =!= None && ! (vr[[1]] < mean < vr[[2]]), Throw[failure["BadModel", "\"Mean\" must lie inside the value range."], "CertiField"]]];
  g = base["G"];
  Do[
    If[Max[Abs[sliceAt[g, n, i, 1] - sliceAt[g, n, i, -1]]] > 10.^-12 (1. + Max[Abs[g]]),
      Throw[failure["BadModel", "The boundary data must be periodic in the periodic directions."], "CertiField"]];
    g[[slab[n, i, -1]]] = g[[slab[n, i, 1]]],
    {i, per}];
  If[vr =!= None,
    gf = SetPrecision[Pick[Flatten[g], Flatten[base["Fixed"]], 1], Infinity];
    If[gf =!= {} && ! (Min[gf] >= vr[[1]] && Max[gf] <= vr[[2]]),
      Throw[failure["BadModel", "The boundary data leave the value range."], "CertiField"]]];
  kind = If[vr === None, "Line", "Range"];
  If[kind === "Line" && Length[P] <= 2 && mean === None && ! MemberQ[Flatten[base["Dirichlet"]], True],
    Throw[failure["BadModel", "With a potential that is at most linear in u, the energy is bounded below only with a Dirichlet face or a prescribed mean."], "CertiField"]];
  Join[base, <|"G" -> g, "P" -> P, "VR" -> vr, "C1" -> c1, "C2" -> c2, "TableMat" -> tab, "TM" -> tm, "Lin" -> lin,
    "LinNodal" -> linNodal, "Av" -> av, "Mean" -> mean, "Periodic" -> Table[MemberQ[per, i], {i, n}], "Kind" -> kind,
    "Linear" -> (kind === "Line" && Length[P] <= 2), "Volume" -> Times @@ (dom[[All, 2]] - dom[[All, 1]]), "ND2" -> True|>]];

nd2Tie[spec_, u0_] := Module[{u = u0, n = spec["n"]},
  Do[If[spec["Periodic"][[i]], u[[slab[n, i, -1]]] = u[[slab[n, i, 1]]]], {i, n}];
  u];
nd2TiedQ[spec_, u_] := AllTrue[Range[spec["n"]], ! spec["Periodic"][[#]] || sliceAt[u, spec["n"], #, 1] === sliceAt[u, spec["n"], #, -1] &];

(* ------------------------------------------------------------------------------------------------ *)
(* RIGOROUS: the potential                                                                          *)
(* ------------------------------------------------------------------------------------------------ *)

(* the smooth part P(t) + c1 (t - lo) Log(t - lo) + c2 (hi - t) Log(hi - t) over interval arrays t inside the range *)
nd2SmoothIv[spec_, t_] := Module[{v = ivHornerQ[spec["P"], t], vr = spec["VR"]},
  If[vr =!= None && spec["C1"] != 0, v = ndAdd[v, ivK[nd2XLnX[ndSub[t, ivC[vr[[1]]]]], spec["C1"]]]];
  If[vr =!= None && spec["C2"] != 0, v = ndAdd[v, ivK[nd2XLnX[ndSub[ivC[vr[[2]]], t]], spec["C2"]]]];
  v];

(* c/s over interval arrays s, c > 0 rational: [c/s_hi, c/s_lo], and $big above where s_lo <= 0 *)
nd2InvPos[{a_, b_}, c_] := Module[{pos = gtMask[a, 0.], a1, b1, q},
  a1 = vWhere[pos, a, 1.]; b1 = vWhere[pos, vMax[b, a1], 1.];
  q = ivK[ndDiv[{1. + 0. a1, 1. + 0. a1}, {a1, b1}], c];
  {vWhere[pos, q[[1]], 0. a1], vWhere[pos, q[[2]], $big + 0. a1]}];

(* W'' of the smooth part over interval arrays t *)
nd2SmoothW2Iv[spec_, t_] := Module[{v = ivHornerQ[spec["W2"], t], vr = spec["VR"]},
  If[vr =!= None && spec["C1"] != 0, v = ndAdd[v, nd2InvPos[ndSub[t, ivC[vr[[1]]]], spec["C1"]]]];
  If[vr =!= None && spec["C2"] != 0, v = ndAdd[v, nd2InvPos[ndSub[ivC[vr[[2]]], t], spec["C2"]]]];
  v];

(* the potential at an exact t, without the data term: a scalar interval *)
nd2SmoothAtQ[spec_, t_] := Module[{lo = polyAt[spec["P"], t], hi, vr = spec["VR"], b},
  hi = lo;
  If[vr =!= None && spec["C1"] != 0, b = spec["C1"] nd2XLnXQ[t - vr[[1]]]; lo += b[[1]]; hi += b[[2]]];
  If[vr =!= None && spec["C2"] != 0, b = spec["C2"] nd2XLnXQ[vr[[2]] - t]; lo += b[[1]]; hi += b[[2]]];
  {fLo[lo], fHi[hi]}];

(* the data term at an exact t in the range, all cells: an interval array with the dimensions of the grid *)
nd2TableAtQ[spec_, t_] := Module[{vr = spec["VR"], m = spec["TM"], tb = spec["TableMat"], z, j, th},
  z = (t - vr[[1]]) m/(vr[[2]] - vr[[1]]);
  j = Clip[Floor[z], {0, m - 1}]; th = z - j;
  ArrayReshape[#, spec["Cells"]] & /@ ndAdd[ivPt[tb[[All, j + 1]]], ivK[ndSub[ivPt[tb[[All, j + 2]]], ivPt[tb[[All, j + 1]]]], th]]];

nd2WAtQ[spec_, t_] := If[spec["TableMat"] === None, nd2SmoothAtQ[spec, t], ndAdd[nd2SmoothAtQ[spec, t], nd2TableAtQ[spec, t]]];

(* the data term over thin interval arrays {a, b} (one interval per cell, flattened): the hull of its values at a, b
   and at the table node between them; at a point t the piece is found from an enclosure of (t - lo)/dt, so the
   true piece is one of the (at most two) pieces evaluated *)
nd2TableIv[spec_, {a_, b_}] := Module[{tb = spec["TableMat"], m = spec["TM"], vr = spec["VR"], nc, dt, at, pa, pb, node, ni},
  nc = Length[a]; dt = (vr[[2]] - vr[[1]])/m;
  at[t_] := Module[{z = ivK[ndSub[ivPt[t], ivC[vr[[1]]]], 1/dt], j1, j2, piece},
    j1 = Clip[Floor[z[[1]]], {0, m - 1}]; j2 = Clip[Floor[z[[2]]], {0, m - 1}];
    piece[j_] := Module[{t0 = Extract[tb, Transpose[{Range[nc], j + 1}]], t1 = Extract[tb, Transpose[{Range[nc], j + 2}]]},
      ndAdd[ivPt[t0], ndMul[ndSub[z, ivPt[N[j]]], ndSub[ivPt[t1], ivPt[t0]]]]];
    With[{p1 = piece[j1], p2 = piece[j2]}, {{vMin[p1[[1]], p2[[1]]], vMax[p1[[2]], p2[[2]]]}, j1, j2}]];
  pa = at[a]; pb = at[b];
  ni = gtMask[pb[[3]], pa[[2]]];
  node = Extract[tb, Transpose[{Range[nc], pb[[3]] + 1}]];
  {vWhere[ni, vMin[vMin[pa[[1, 1]], pb[[1, 1]]], node], vMin[pa[[1, 1]], pb[[1, 1]]]],
   vWhere[ni, vMax[vMax[pa[[1, 2]], pb[[1, 2]]], node], vMax[pa[[1, 2]], pb[[1, 2]]]]}];

(* Proposition 13.3: the integral of D_c(u_h) over each cell divided by |c|, enclosed: D_c = W+ - W- with W- the sum
   of the concave kinks; between D_c(centre) - e and the vertex mean of D_c + e, e = mean_v W-(U_v) - W-(centre) *)
nd2TableCells[spec_, uc_, mean_, verts_] := Module[{n = spec["n"], tb = spec["TableMat"], m = spec["TM"], vr = spec["VR"], dt, flat, vv, ce, em = 0., tj, jump, cj, br},
  dt = (vr[[2]] - vr[[1]])/m;
  flat[iv_] := Flatten /@ iv;
  vv = ivK[ivSum[nd2TableIv[spec, flat[uc[#]]] & /@ verts], 1/2^n];
  ce = nd2TableIv[spec, flat[mean]];
  Do[
    tj = ivC[vr[[1]] + j dt];
    jump = ivK[ndAdd[ndSub[ivPt[tb[[All, j + 2]]], ivK[ivPt[tb[[All, j + 1]]], 2]], ivPt[tb[[All, j]]]], 1/dt];
    cj = Ramp[-jump[[1]]];
    br = ndSub[ivK[ivSum[ivPos[ndSub[flat[uc[#]], tj]] & /@ verts], 1/2^n], ivPos[ndSub[flat[mean], tj]]];
    em = ndUp[em + ndMul[ivPt[cj], ivPt[Ramp[br[[2]]]]][[2]]],
    {j, 1, m - 1}];
  ArrayReshape[#, spec["Cells"]] & /@ {ndDn[ce[[1]] - em], ndUp[vv[[2]] + em]}];

(* ------------------------------------------------------------------------------------------------ *)
(* RIGOROUS: the competitor                                                                          *)
(* ------------------------------------------------------------------------------------------------ *)

(* the nodal values of the competitor as intervals.  With a prescribed mean the competitor is u_h + delta phi_h, phi_h
   the Q1 function that is 1 at the nodes off the Dirichlet faces, delta exact (Proposition 13.4) *)
nd2NodalEnclosure[spec_, u_] := Module[{n = spec["n"], verts, phi, iu, iphi, delta},
  If[spec["Mean"] === None, Return[{u, u}]];
  verts = Tuples[{0, 1}, n];
  phi = 1 - spec["Fixed"];
  iu = spec["CellVolume"]/2^n exactTotal[Flatten[vertexPart[u, n, #] & /@ verts]];
  iphi = spec["CellVolume"]/2^n Total[Flatten[vertexPart[phi, n, #] & /@ verts]];
  If[iphi == 0, Throw[failure["BadModel", "The mean cannot be prescribed: every node lies on a Dirichlet face."], "CertiField"]];
  delta = (spec["Mean"] spec["Volume"] - iu)/iphi;
  {vWhere[phi, ndDn[u + fLo[delta]], u], vWhere[phi, ndUp[u + fHi[delta]], u]}];

(* RIGOROUS enclosure {lo, hi} (exact rationals) of the energy of the competitor (Propositions 7 and 13) *)
nd2UpperBound[spec_, u0_] := ndQuiet@Module[
  {n = spec["n"], vr = spec["VR"], u, verts, lohi, uc, invh, h2, umin, umax, corners, d, acc, kin = None, dmax2, ri, w2, fin, w2pm,
   remLo, remHi, mean, mid, cvx, vtx, potLo, potHi, tab, lin, cellLo, cellHi},
  verts = Tuples[{0, 1}, n];
  u = ndRealArray[u0, spec["Cells"] + 1];
  If[u === $Failed, Throw[failure["BadCandidate", "The nodal values must be a real array of dimensions Cells + 1."], "CertiField"]];
  If[! (Max[Abs[spec["Fixed"] (u - spec["G"])]] == 0.), Throw[failure["Internal", "The Dirichlet nodal values are not imposed exactly."], "CertiField"]];
  If[! nd2TiedQ[spec, u], Throw[failure["Internal", "The nodal values differ on identified periodic nodes."], "CertiField"]];
  lohi = nd2NodalEnclosure[spec, u];
  If[vr =!= None && ! (exactNumber[Min[lohi[[1]]]] >= vr[[1]] && exactNumber[Max[lohi[[2]]]] <= vr[[2]]),
    Throw[failure["Inconclusive", "The candidate leaves the value range of the potential, where its energy is infinite."], "CertiField"]];
  uc = Association[# -> {vertexPart[lohi[[1]], n, #], vertexPart[lohi[[2]], n, #]} & /@ verts];
  invh = ivC[1/#] & /@ spec["h"];
  h2 = ivC[#^2] & /@ spec["h"];
  umin = Fold[vMin, uc[#][[1]] & /@ verts];
  umax = Fold[vMax, uc[#][[2]] & /@ verts];
  w2 = If[vr === None, ivHornerQ[spec["W2"], {umin, umax}], nd2SmoothW2Iv[spec, {umin, umax}]];
  fin = ltMask[w2[[2]], $big/2];
  w2pm = {vWhere[fin, vMin[w2[[1]], 0.], 0. umin], vWhere[fin, vMax[w2[[2]], 0.], 0. umin]};
  corners = Tuples[{0, 1}, n - 1];
  remLo = 0. umin; remHi = 0. umin;
  Do[
    d = Association[# -> ndMul[ndSub[uc[Insert[#, 1, i]], uc[Insert[#, 0, i]]], invh[[i]]] & /@ corners];
    acc = ivSum[Flatten[Table[
      ivK[If[e === e2, ivSqr[d[e]], ndMul[d[e], d[e2]]], Times @@ MapThread[If[#1 == #2, 1/3, 1/6] &, {e, e2}]],
      {e, corners}, {e2, corners}], 1]];
    kin = If[kin === None, acc, ndAdd[kin, acc]];
    dmax2 = Fold[vMax, ivSqr[d[#]][[2]] & /@ corners];
    ri = ivK[ndMul[ndMul[w2pm, {0. dmax2, dmax2}], h2[[i]]], 1/24];
    remLo = ndDn[remLo + ri[[1]]]; remHi = ndUp[remHi + ri[[2]]],
    {i, n}];
  kin = ivK[kin, spec["k"]];
  mean = ivK[ivSum[uc /@ verts], 1/2^n];
  If[vr === None,
    mid = ivHornerQ[spec["P"], mean];
    potLo = ndDn[mid[[1]] + remLo]; potHi = ndUp[mid[[2]] + remHi],
    mid = nd2SmoothIv[spec, mean];
    cvx = geMask[w2[[1]], 0.];
    vtx = ivK[ivSum[nd2SmoothIv[spec, uc[#]] & /@ verts], 1/2^n];
    potHi = vMin[vWhere[fin, ndUp[mid[[2]] + remHi], $big + 0. umin], vWhere[cvx, vtx[[2]], $big + 0. umin]];
    potLo = vMax[vWhere[fin, ndDn[mid[[1]] + remLo], -$big + 0. umin], vWhere[cvx, mid[[1]], -$big + 0. umin]];
    If[Max[potHi] >= $big/2,
      Throw[failure["Inconclusive", "The candidate comes so close to an end of the value range that its energy cannot be enclosed there."], "CertiField"]];
    If[spec["TableMat"] =!= None,
      tab = nd2TableCells[spec, uc, mean, verts];
      potLo = ndDn[potLo + tab[[1]]]; potHi = ndUp[potHi + tab[[2]]]]];
  If[spec["Av"] =!= None,
    lin = ivSum[Flatten[Table[
      ivK[ndMul[ivPt[spec["Av"][v]], uc[v2]], Times @@ MapThread[If[#1 == #2, 1/3, 1/6] &, {v, v2}]], {v, verts}, {v2, verts}], 1]];
    potLo = ndDn[potLo + lin[[1]]]; potHi = ndUp[potHi + lin[[2]]]];
  cellLo = ndDn[kin[[1]] + potLo]; cellHi = ndUp[kin[[2]] + potHi];
  spec["CellVolume"] {exactTotal[cellLo], exactTotal[cellHi]}];

(* ------------------------------------------------------------------------------------------------ *)
(* RIGOROUS: the lower bound                                                                         *)
(* ------------------------------------------------------------------------------------------------ *)

(* the field with packed arrays and a multiplier "Lambda" (a machine number, 0 without a prescribed mean) *)
nd2PackField[spec_, f_] := Module[{g = ndPackField[spec, f], lam},
  If[g === $Failed, Return[$Failed]];
  lam = N[Lookup[f, "Lambda", 0.]];
  If[! MachineNumberQ[lam], Return[$Failed]];
  Append[g, "Lambda" -> lam]];

(* zero normal component on the natural faces, equal values on identified periodic faces *)
nd2RespectsFaces[spec_, f_] := And @@ Flatten[Table[
  If[spec["Periodic"][[i]],
    {sliceAt[f["S0"][[i]], spec["n"], i, 1] === sliceAt[f["S0"][[i]], spec["n"], i, -1] &&
      (f["K"] == 0 || sliceAt[f["C"][[i]], spec["n"], i, 1] === sliceAt[f["C"][[i]], spec["n"], i, -1])},
    Table[With[{sp = If[side == 1, 1, -1]},
      spec["Dirichlet"][[i, side]] || (Max[Abs[sliceAt[f["S0"][[i]], spec["n"], i, sp]]] == 0. &&
        (f["K"] == 0 || Max[Abs[sliceAt[f["C"][[i]], spec["n"], i, sp]]] == 0.))], {side, 2}]],
  {i, spec["n"]}]];

(* an upper bound h2 for the regular part of W'' on [a, b] and the corrections ca, cb of Lemma 11 at a singular end *)
nd2SubBound[spec_, a_, b_] := Module[{vr = spec["VR"], h = Last[exactRange[spec["W2"], a, b]], ca = 0, cb = 0},
  If[vr =!= None && spec["C1"] != 0, If[a > vr[[1]], h += spec["C1"]/(a - vr[[1]]), ca = spec["C1"] (b - a)]];
  If[vr =!= None && spec["C2"] != 0, If[b < vr[[2]], h += spec["C2"]/(vr[[2]] - b), cb = spec["C2"] (b - a)]];
  {fHi[h], fHi[ca], fHi[cb]}];

Options[nd2LowerBound] = {"Subintervals" -> 16};
nd2LowerBound[spec_, f0_, OptionsPattern[]] := ndQuiet[nd2LowerBoundRaw[spec, f0, OptionValue["Subintervals"]]];
nd2LowerBoundRaw[spec_, f0_, qs_] := Module[
  {n = spec["n"], vr = spec["VR"], f, kk, t0, t1, lam, verts, invh, inv4k, s0, d0, a0, tilt, tl, th, flux, mup, dl, w, xi, cc, c, d01, d12, dd2,
   prev, cur, tt, theta, n0, n1, n2, w0, w1, dv, aa, base, ta, h2, ca, cb, corr, dend, intM, args},
  f = nd2PackField[spec, f0];
  If[f === $Failed || ! nd2RespectsFaces[spec, f],
    Throw[failure["BadField", "The field does not fit the grid, has a normal component on a natural face, or differs on identified periodic faces."], "CertiField"]];
  If[! (IntegerQ[qs] && qs >= 1), Throw[failure["BadOption", "\"Subintervals\" must be a positive integer."], "CertiField"]];
  kk = f["K"]; {t0, t1} = f["TRange"]; lam = f["Lambda"];
  If[spec["Mean"] === None && lam != 0., Throw[failure["BadField", "A multiplier needs a prescribed mean."], "CertiField"]];
  If[spec["Linear"],
    If[kk != 0, Throw[failure["BadField", "A potential that is at most linear in u needs a field with K = 0."], "CertiField"]];
    Return[nd2PoissonLower[spec, f]]];
  If[kk == 0, Throw[failure["NotImplemented", "With these model keys the lower bound needs a calibration field (K > 0) unless the potential is at most linear."], "CertiField"]];
  If[vr =!= None && {t0, t1} =!= vr, Throw[failure["BadField", "With a value range the t-range of the field must be the value range."], "CertiField"]];
  If[vr =!= None && kk qs < 2, Throw[failure["BadOption", "Need at least two subintervals in the value range."], "CertiField"]];
  If[! f["Affine"] && Max[Abs[Flatten[f["S0"]]]] != 0., Throw[failure["BadField", "A field without affine part must have S0 = 0."], "CertiField"]];
  If[vr === None && ! superlinearQ[spec["P"]],
    Throw[failure["Inconclusive", "The potential must have even degree and a positive leading coefficient, or a value range."], "CertiField"]];
  If[spec["TableMat"] =!= None && ! Divisible[kk qs, spec["TM"]],
    Throw[failure["BadOption", "The table nodes of the data term must be subinterval ends: M must divide K times \"Subintervals\"."], "CertiField"]];
  verts = Tuples[{0, 1}, n];
  invh = ivC[1/#] & /@ spec["h"];
  inv4k = ivC[1/(4 spec["k"])];
  s0 = ivPt /@ f["S0"];
  d0 = ivSum[Table[ndMul[ndSub[ivHi[s0[[i]], n, i], ivLo[s0[[i]], n, i]], invh[[i]]], {i, n}]];
  a0 = Table[With[{q = ndMul[ivSqr[s0[[i]]], inv4k]}, {ivLo[q, n, i], ivHi[q, n, i]}], {i, n}];
  (* the factor of t at the vertex v of every cell: lam - a_v (Proposition 10, Corollary 9) *)
  tilt = Association[# -> If[spec["Av"] === None, {lam + 0. d0[[1]], lam + 0. d0[[1]]}, ndSub[{lam, lam}, ivPt[spec["Av"][#]]]] & /@ verts];
  If[vr === None,
    (* Lemma 4 at the shifted arguments d0 + lam - a_v, all vertices in one table *)
    args = Developer`ToPackedArray /@ Transpose[ndAdd[d0, tilt[#]] & /@ verts];
    tl = conjCells[spec["P"], args, t0, -1]; th = conjCells[spec["P"], args, t1, +1]];
  flux = ndFluxLower[spec, f];
  dl = (t1 - t0)/kk;
  w = dl/qs;
  xi = t0 + (Range[0, kk + 1] - 1/2) dl;
  cc = f["C"];
  mup = Association[# -> None & /@ verts];
  Do[
    c = Table[ndAdd[ivPt[lastPart[cc[[i]], n, kap + r - 1]], ivK[s0[[i]], xi[[kap + r - 1]]]], {i, n}, {r, 3}];
    d01 = Table[ndSub[c[[i, 2]], c[[i, 1]]], {i, n}];
    d12 = Table[ndSub[c[[i, 3]], c[[i, 2]]], {i, n}];
    dd2 = ivSum[Table[
      With[{s2 = ivK[ndSub[d12[[i]], d01[[i]]], 1/dl^2]}, ndMul[ndSub[ivHi[s2, n, i], ivLo[s2, n, i]], invh[[i]]]], {i, n}]];
    prev = None;
    Do[
      theta = j/qs;
      tt = t0 + (kap - 1 + theta) dl;
      {n0, n1, n2} = {(1 - theta)^2/2, (-2 theta^2 + 2 theta + 1)/2, theta^2/2};
      {w0, w1} = {1 - theta, theta};
      dv = ivSum[Table[
        With[{s = ivSum[{ivK[c[[i, 1]], n0], ivK[c[[i, 2]], n1], ivK[c[[i, 3]], n2]}]},
          ndMul[ndSub[ivHi[s, n, i], ivLo[s, n, i]], invh[[i]]]], {i, n}]];
      aa = Table[
        With[{q = ndMul[ivSqr[ivK[ndAdd[ivK[d01[[i]], w0], ivK[d12[[i]], w1]], 1/dl]], inv4k]}, {ivLo[q, n, i], ivHi[q, n, i]}], {i, n}];
      base = ndSub[dv, nd2WAtQ[spec, tt]];
      cur = Association[# -> ivSum[Join[{base, ndMul[tilt[#], ivC[tt]]}, Table[aa[[i, #[[i]] + 1]], {i, n}]]][[2]] & /@ verts];
      If[prev =!= None,
        ta = t0 + (kap - 1 + (j - 1)/qs) dl;
        {h2, ca, cb} = nd2SubBound[spec, ta, tt];
        corr = ivK[ivPos[ndSub[ivPt[h2 + 0. dd2[[1]]], dd2]], w^2/8][[2]];
        Do[
          With[{cand = ndUp[vMax[ndUp[prev[v] + ca], ndUp[cur[v] + cb]] + corr]},
            mup[v] = If[mup[v] === None, cand, vMax[mup[v], cand]]], {v, verts}]];
      prev = cur,
      {j, 0, qs}],
    {kap, kk}];
  If[vr === None,
    dend = Table[
      ivSum[Table[
        With[{ev = ivK[ndAdd[ivPt[lastPart[cc[[i]], n, m0]], ivPt[lastPart[cc[[i]], n, m0 + 1]]], 1/2]},
          ndMul[ndSub[ivHi[ev, n, i], ivLo[ev, n, i]], invh[[i]]]], {i, n}]],
      {m0, {1, kk + 1}}];
    Do[
      With[{v = verts[[vi]], av = ivSum[Table[a0[[i, verts[[vi, i]] + 1]], {i, n}]]},
        mup[v] = vMax[mup[v], ivSum[{dend[[1]], av, ivPt[tl[[vi]]]}][[2]]];
        mup[v] = vMax[mup[v], ivSum[{dend[[2]], av, ivPt[th[[vi]]]}][[2]]]],
      {vi, Length[verts]}]];
  intM = spec["CellVolume"]/2^n exactTotal[ndChk1 /@ Values[mup]];
  flux - intM + If[spec["Mean"] === None, 0, exactNumber[lam] spec["Mean"] spec["Volume"]]];

(* RIGOROUS: Proposition 16, a potential at most linear in u, W = w0 + (w1 + a(x)) t, for any Raviart-Thomas field *)
nd2PoissonLower[spec_, f_] := Module[
  {n = spec["n"], cv = spec["CellVolume"], verts, invh, s0, d0, sig2, abar, rho, r1sq, srho, ssig, sr, ipi, flux, lam = f["Lambda"],
   w0, w1, c2, c0, bigB, pen, faces, li, side, gmax, sumRho, cst},
  verts = Tuples[{0, 1}, n];
  {w0, w1} = PadRight[spec["P"], 2];
  invh = ivC[1/#] & /@ spec["h"];
  s0 = ivPt /@ f["S0"];
  d0 = ivSum[Table[ndMul[ndSub[ivHi[s0[[i]], n, i], ivLo[s0[[i]], n, i]], invh[[i]]], {i, n}]];
  sig2 = ivSum[Table[With[{lo = ivLo[s0[[i]], n, i], hi = ivHi[s0[[i]], n, i]},
    ivK[ivSum[{ivSqr[lo], ndMul[lo, hi], ivSqr[hi]}], 1/3]], {i, n}]];
  abar = If[spec["Av"] === None, {0. d0[[1]], 0. d0[[1]]}, ivK[ivSum[ivPt[spec["Av"][#]] & /@ verts], 1/2^n]];
  rho = ndSub[ndSub[d0, abar], ndSub[ivC[w1], {lam, lam}]];
  r1sq = If[spec["Av"] === None || ! spec["LinNodal"], 0,
    With[{q = ndSub[ivSum[Flatten[Table[ivK[ndMul[ivPt[spec["Av"][v]], ivPt[spec["Av"][v2]]], Times @@ MapThread[If[#1 == #2, 1/3, 1/6] &, {v, v2}]],
        {v, verts}, {v2, verts}], 1]], ivSqr[abar]]},
      cv exactTotal[Ramp[q[[2]]]]]];
  srho = nd2SqrtUp[fHi[cv exactTotal[ivSqr[rho][[2]]]]];
  ssig = nd2SqrtUp[fHi[cv exactTotal[sig2[[2]]]]];
  sr = nd2SqrtUp[fHi[r1sq]];
  ipi = fHi[1/nd2PiLower[]];
  flux = ndFluxLower[spec, f];
  cst = w0 spec["Volume"];
  If[spec["Mean"] =!= None,
    c2 = rUp[fHi[Max[spec["Domain"][[All, 2]] - spec["Domain"][[All, 1]]]] ipi];
    sumRho = cv {exactTotal[rho[[1]]], exactTotal[rho[[2]]]};
    c0 = Max[spec["Mean"] sumRho];
    cst += exactNumber[lam] spec["Mean"] spec["Volume"],
    faces = SortBy[Select[Flatten[Table[{spec["Domain"][[i, 2]] - spec["Domain"][[i, 1]], i, side}, {i, n}, {side, 2}], 1], spec["Dirichlet"][[#[[2]], #[[3]]]] &], First];
    {li, side} = Rest[First[faces]];
    c2 = nd2SqrtUp[fHi[(spec["Domain"][[li, 2]] - spec["Domain"][[li, 1]])^2/2]];
    gmax = Max[Abs[sliceAt[spec["G"], n, li, If[side == 1, 1, -1]]]];
    c0 = exactNumber[rUp[rUp[srho nd2SqrtUp[fHi[spec["Volume"]]]] gmax]]];
  bigB = rUp[rUp[ssig + rUp[rUp[fHi[Max[spec["h"]]] ipi] sr]] + rUp[c2 srho]];
  pen = rUp[rUp[bigB bigB] fHi[1/(4 spec["k"])]];
  cst + flux - exactNumber[pen] - c0];

(* ------------------------------------------------------------------------------------------------ *)
(* heuristic: the equilibrated flux of a potential at most linear in u (mixed Raviart-Thomas solve)   *)
(* ------------------------------------------------------------------------------------------------ *)

(* maximise  flux - |sigma|^2/(4k) + lam m |Omega|  over Raviart-Thomas fields with div sigma = A - lam on every cell,
   A = w1 + (cell mean of a): the bound of Proposition 16 with rho = 0, a linear system *)
nd2MixedField[spec_] := Module[
  {n = spec["n"], cells = spec["Cells"], hf = N[spec["h"]], cv = N[spec["CellVolume"]], k = N[spec["k"]], fs, nf, nc, lo, hi, map, nu = 0,
   sel, plo, phi, mm, bb, fl, avec, kkt, rhs, sol, col, lam = 0., w1},
  nc = Times @@ cells;
  fs = Table[faceShape[spec, i], {i, n}]; nf = Times @@@ fs;
  lo = Table[1 + Flatten[Outer[Plus, Sequence @@ Table[Range[0, cells[[j]] - 1] (Times @@ fs[[i, j + 1 ;;]]), {j, n}]]], {i, n}];
  hi = Table[lo[[i]] + Times @@ fs[[i, i + 1 ;;]], {i, n}];
  map = Table[
    Module[{keep = ConstantArray[1, fs[[i]]], idx, first, last, all = ArrayReshape[Range[nf[[i]]], fs[[i]]]},
      If[spec["Periodic"][[i]],
        keep[[slab[n, i, -1]]] = 0,
        If[! spec["Dirichlet"][[i, 1]], keep[[slab[n, i, 1]]] = 0];
        If[! spec["Dirichlet"][[i, 2]], keep[[slab[n, i, -1]]] = 0]];
      keep = Flatten[keep];
      idx = (Accumulate[keep] + nu) keep;
      nu = Max[nu, Max[idx]];
      If[spec["Periodic"][[i]],
        first = Flatten[all[[slab[n, i, 1]]]]; last = Flatten[all[[slab[n, i, -1]]]];
        idx[[last]] = idx[[first]]];
      idx],
    {i, n}];
  sel[i_, which_] := Module[{uu = map[[i]][[which[[i]]]], keep},
    keep = Flatten[Position[uu, _?Positive]];
    If[keep === {}, SparseArray[{}, {nc, nu}], SparseArray[Transpose[{keep, uu[[keep]]}] -> 1., {nc, nu}]]];
  plo = Table[sel[i, lo], {i, n}]; phi = Table[sel[i, hi], {i, n}];
  mm = cv Sum[(Transpose[plo[[i]]] . plo[[i]] + Transpose[phi[[i]]] . phi[[i]])/3. + (Transpose[plo[[i]]] . phi[[i]] + Transpose[phi[[i]]] . plo[[i]])/6., {i, n}];
  bb = Sum[(phi[[i]] - plo[[i]])/hf[[i]], {i, n}];
  fl = ConstantArray[0., nu];
  Do[
    If[spec["Dirichlet"][[i, side]],
      With[{ids = Flatten[ArrayReshape[map[[i]], fs[[i]]][[slab[n, i, If[side == 1, 1, -1]]]]], gb = Flatten[faceMeans[spec, i, side]]},
        fl[[ids]] += If[side == 1, -1., 1.] (cv/hf[[i]]) gb]],
    {i, n}, {side, 2}];
  w1 = N[PadRight[spec["P"], 2][[2]]];
  avec = w1 + If[spec["Av"] === None, ConstantArray[0., nc], Flatten[Mean[Values[spec["Av"]]]]];
  If[spec["Mean"] === None,
    kkt = ArrayFlatten[{{mm/(2 k), cv Transpose[bb]}, {cv bb, SparseArray[{}, {nc, nc}]}}];
    rhs = Join[fl, cv avec],
    col = SparseArray[Table[{j, 1} -> cv, {j, nc}], {nc, 1}];
    kkt = ArrayFlatten[{{mm/(2 k), cv Transpose[bb], SparseArray[{}, {nu, 1}]}, {cv bb, SparseArray[{}, {nc, nc}], col},
      {SparseArray[{}, {1, nu}], Transpose[col], SparseArray[{}, {1, 1}]}}];
    rhs = Join[fl, cv avec, {N[spec["Mean"] spec["Volume"]]}]];
  sol = Quiet@Check[LinearSolve[kkt, rhs], $Failed];
  If[! VectorQ[sol, MachineNumberQ], Throw[failure["Inconclusive", "The equilibrated flux could not be computed."], "CertiField"]];
  If[spec["Mean"] =!= None, lam = sol[[-1]]];
  <|"TRange" -> {0, 0}, "K" -> 0, "Affine" -> True, "C" -> ConstantArray[{}, n], "Lambda" -> lam,
    "S0" -> Table[Developer`ToPackedArray[ArrayReshape[With[{m = map[[i]]}, (Unitize[m] sol[[Clip[m, {1, nu}]]])], fs[[i]]]], {i, n}]|>];

(* ------------------------------------------------------------------------------------------------ *)
(* heuristic: the candidate                                                                          *)
(* ------------------------------------------------------------------------------------------------ *)

(* the potential, its derivative and a curvature for Newton steps, at machine points t with their cells cg *)
nd2WMach[spec_, t_, cg_] := Module[{pf = N[spec["P"]], w, w1, w2, vr = spec["VR"], lo, hi, s, tb = spec["TableMat"], m, dt, z, j, th, a0, a1, jn},
  w = polyAt[pf, t]; w1 = polyAt[polyD[pf], t]; w2 = polyAt[polyD[polyD[pf]], t];
  If[vr =!= None,
    {lo, hi} = N[vr];
    If[spec["C1"] != 0, s = Ramp[t - lo] + 10.^-300;
      w += N[spec["C1"]] s Log[s]; w1 += N[spec["C1"]] (Log[s] + 1.); w2 += N[spec["C1"]]/s];
    If[spec["C2"] != 0, s = Ramp[hi - t] + 10.^-300;
      w += N[spec["C2"]] s Log[s]; w1 -= N[spec["C2"]] (Log[s] + 1.); w2 += N[spec["C2"]]/s];
    If[tb =!= None,
      m = spec["TM"]; dt = (hi - lo)/m;
      z = Clip[(t - lo)/dt, {0., N[m]}]; j = Clip[Floor[z], {0, m - 1}]; th = z - j;
      a0 = Extract[tb, Transpose[{cg, j + 1}]]; a1 = Extract[tb, Transpose[{cg, j + 2}]];
      w += a0 + th (a1 - a0); w1 += (a1 - a0)/dt;
      If[m >= 2,
        jn = Clip[Round[z], {1, m - 1}];
        w2 += Ramp[(Extract[tb, Transpose[{cg, jn + 2}]] - 2 Extract[tb, Transpose[{cg, jn + 1}]] + Extract[tb, Transpose[{cg, jn}]])/dt^2]]]];
  {w, w1, w2}];

nd2WSmoothMach[spec_, t_] := Module[{w = polyAt[N[spec["P"]], t], vr = spec["VR"]},
  If[vr =!= None,
    If[spec["C1"] != 0, w += N[spec["C1"]] nd2XLogMach[t - N[vr[[1]]]]];
    If[spec["C2"] != 0, w += N[spec["C2"]] nd2XLogMach[N[vr[[2]]] - t]]];
  w];

(* the cell of every Gauss point (rows of the Kronecker product of the 1D Gauss matrices) *)
nd2GaussCells[cells_] := Module[{n = Length[cells], strides},
  strides = Reverse[FoldList[Times, 1, Reverse[Rest[cells]]]];
  1 + Flatten[Outer[Plus, Sequence @@ Table[(Ceiling[Range[2 cells[[i]]]/2] - 1) strides[[i]], {i, n}]]]];

(* nodes of a periodic direction: the last one is the first one.  Returns {map full -> reduced, fixed, data} *)
nd2NodeMap[spec_] := Module[{n = spec["n"], dims = spec["Cells"] + 1, per = spec["Periodic"], rdims, strides, red, nfull, nred, rr, fx, cnt},
  rdims = MapThread[If[#2, #1 - 1, #1] &, {dims, per}];
  strides = Reverse[FoldList[Times, 1, Reverse[Rest[rdims]]]];
  red = (1 + (MapThread[If[#2 && #1 == #3, 1, #1] &, {#, per, dims}] - 1) . strides) & /@ Tuples[Range /@ dims];
  nfull = Length[red]; nred = Times @@ rdims;
  rr = SparseArray[Transpose[{Range[nfull], red}] -> 1., {nfull, nred}];
  fx = N[Flatten[spec["Fixed"]]];
  cnt = Transpose[rr] . fx;
  {rr, Unitize[cnt], (Transpose[rr] . (fx Flatten[spec["G"]]))/(Max[1., #] & /@ cnt)}];

nd2KKTSolve[h_, w_, r_] := Module[{m = Length[r], sol},
  sol = Quiet@Check[LinearSolve[ArrayFlatten[{{h, SparseArray[Transpose[{w}]]}, {SparseArray[{w}], SparseArray[{}, {1, 1}]}}], Append[r, 0.]], $Failed];
  If[VectorQ[sol, MachineNumberQ], sol[[;; m]], $Failed]];

(* heuristic: a local minimiser of the Q1 energy with two-point Gauss quadrature (Newton steps with the positive part
   of W'', an Armijo search, the mean as a linear constraint, nodal values kept inside the value range) *)
nd2SolveCandidate[spec_, u0_] := Module[
  {n = spec["n"], cells = spec["Cells"], k = N[spec["k"]], wq = N[spec["CellVolume"]]/2^spec["n"], mats, pg, gm, stiff, cg, ag, rr, fixedR, gR,
   pgR, pgRT, stR, free, x, wm, en, grad, hess, d, e0, e1, t, xn, it = 0, slope, wmv, vr = spec["VR"], eps, lo, hi, tmax, mean = spec["Mean"], maxit, ug, wv},
  mats = Table[ndGauss1D[spec, i], {i, n}];
  pg = KroneckerProduct @@ mats[[All, 1]];
  gm = Table[KroneckerProduct @@ ReplacePart[mats[[All, 1]], i -> mats[[i, 2]]], {i, n}];
  stiff = 2 k wq Total[Transpose[#] . # & /@ gm];
  cg = nd2GaussCells[cells];
  ag = Which[spec["Av"] === None, 0., spec["LinNodal"], pg . Flatten[-spec["Lin"]], True, Flatten[-spec["Lin"]][[cg]]];
  {rr, fixedR, gR} = nd2NodeMap[spec];
  pgR = pg . rr; pgRT = Transpose[pgR]; stR = Transpose[rr] . stiff . rr;
  free = Flatten[Position[fixedR, 0]];
  x = (Transpose[rr] . N[Flatten[u0]])/(Transpose[rr] . ConstantArray[1., Length[Flatten[u0]]]);
  x = fixedR gR + (1 - fixedR) x;
  If[free === {}, Return[ArrayReshape[rr . x, cells + 1]]];
  wm = wq Total[pgR];
  eps = If[vr =!= None && (spec["C1"] != 0 || spec["C2"] != 0), 10.^-9 N[vr[[2]] - vr[[1]]], 0.];
  If[vr =!= None, {lo, hi} = N[vr] + {eps, -eps}; x[[free]] = Clip[x[[free]], {lo, hi}]];
  If[mean =!= None,
    Do[x[[free]] += (N[mean spec["Volume"]] - wm . x)/Total[wm[[free]]]; If[vr =!= None, x[[free]] = Clip[x[[free]], {lo, hi}]], {2}]];
  en[y_] := Module[{yg = pgR . y}, wq Total[nd2WMach[spec, yg, cg][[1]]] + 0.5 y . (stR . y) + If[ag === 0., 0., wq ag . yg]];
  e0 = en[x];
  maxit = If[spec["Linear"], 5, 300];
  While[it++ < maxit,
    ug = pgR . x;
    wv = nd2WMach[spec, ug, cg];
    grad = (wq (pgRT . (wv[[2]] + ag)) + stR . x)[[free]];
    hess = (stR + wq pgRT . (DiagonalMatrix[SparseArray[Ramp[wv[[3]]] + 10.^-10 (1. + Max[Abs[wv[[3]]]])]] . pgR))[[free, free]];
    d = If[mean === None, ndSolveSparse[hess, -grad], nd2KKTSolve[hess, wm[[free]], -grad]];
    If[d === $Failed, Break[]];
    slope = grad . d;
    If[! (slope < 0.) || Max[Abs[d]] <= 10.^-13 (1. + Max[Abs[x]]), Break[]];
    tmax = 1.;
    If[vr =!= None,
      tmax = Min[1., 0.99 Min[Join[{2.}, Pick[hi - x[[free]], Positive[d]]/Pick[d, Positive[d]], Pick[lo - x[[free]], Negative[d]]/Pick[d, Negative[d]]]]]];
    If[! (tmax > 0.), Break[]];
    t = tmax; xn = x; xn[[free]] = x[[free]] + t d; e1 = en[xn];
    While[! (e1 <= e0 + 10.^-4 t slope) && t > 10.^-10 tmax, t *= 0.5; xn = x; xn[[free]] = x[[free]] + t d; e1 = en[xn]];
    If[! (e1 < e0), Break[]];
    x = xn;
    If[e0 - e1 <= 10.^-15 Abs[e0], Break[]];
    e0 = e1];
  ArrayReshape[rr . x, cells + 1]];

nd2AutoStart[spec_, gfun_] := Module[{n = spec["n"], dims = spec["Cells"] + 1, u = $Failed, vr = spec["VR"], tb, m},
  Which[
    spec["Mean"] =!= None, u = ConstantArray[N[spec["Mean"]], dims],
    spec["TableMat"] =!= None,
      tb = spec["TableMat"]; m = spec["TM"];
      u = ndNodalFromCells[spec, ArrayReshape[N[vr[[1]]] + (First[Ordering[#, 1]] & /@ tb - 1) N[(vr[[2]] - vr[[1]])/m], spec["Cells"]]],
    gfun =!= None, u = ndRealArray[Quiet@Apply[gfun, Outer[List, Sequence @@ spec["Axes"]], {n}], dims]];
  If[u === $Failed, u = ConstantArray[If[vr === None, 0., N[Mean[vr]]], dims]];
  If[! spec["Linear"], u = u + 10.^-3 (1. + Max[Abs[u]]) BlockRandom[SeedRandom[20261124]; RandomReal[{-1, 1}, dims]]];
  u];

nd2Competitor[spec_, cand_, refine_, gfun_] := Module[{n = spec["n"], dims = spec["Cells"] + 1, u},
  u = If[cand === Automatic, nd2AutoStart[spec, gfun],
    ndRealArray[If[ArrayQ[cand, n] && Dimensions[cand] === dims, cand, Quiet@Apply[cand, Outer[List, Sequence @@ spec["Axes"]], {n}]], dims]];
  If[u === $Failed,
    Throw[failure["BadCandidate", "The candidate must be a real function of the coordinates or an array of nodal values of dimensions Cells + 1."], "CertiField"]];
  u = ndImposeDirichlet[spec, nd2Tie[spec, u]];
  If[cand === Automatic || TrueQ[refine], u = nd2SolveCandidate[spec, u]];
  ndImposeDirichlet[spec, nd2Tie[spec, u]]];

(* ------------------------------------------------------------------------------------------------ *)
(* heuristic: calibration mode with a linear term, a multiplier, periodic faces, a value range        *)
(* ------------------------------------------------------------------------------------------------ *)

(* float restricted conjugate sup {t s - P(t) : t >= t0 (side +1) or t <= t0 (side -1)} and its maximiser, tabulated *)
nd2TailFloat[pf_, t0_, side_, {s1_, s2_}] := Module[{g = Subdivide[s1, s2, 256], wr, args, vals},
  args = If[side > 0, conjArgmax[pf, #, t0] & /@ g,
    wr = pf (-1.)^Range[0, Length[pf] - 1]; -(conjArgmax[wr, -#, -t0] & /@ g)];
  vals = args g - polyAt[pf, args];
  Function[s, Module[{z = Clip[(s - s1)/(s2 - s1) 256., {0., 255.999999}], j, th},
    j = Floor[z]; th = z - j;
    {vals[[j + 1]] (1 - th) + vals[[j + 2]] th, args[[j + 1]] (1 - th) + args[[j + 2]] th}]]];

nd2CalSetup[spec_, trange_, kk_, qo_, lam0_] := Module[
  {n = spec["n"], vr = spec["VR"], st, nc, verts, fs, mask, ties, tq, q, w, vq, fine, p2, e1, e2, segmax, w2max, endAdd, tabM, tilt, tails = None, sa, pad, ns},
  st = ndCalSetup[Join[spec, <|"W" -> {0, 0, 1}, "W2" -> {2}|>], trange, kk, qo];
  nc = st["NC"]; verts = st["Verts"];
  fs = Table[faceShape[spec, i], {i, n}];
  mask = Table[Module[{m = ConstantArray[1., fs[[i]]]},
      If[! spec["Periodic"][[i]],
        If[! spec["Dirichlet"][[i, 1]], m[[slab[n, i, 1]]] = 0.];
        If[! spec["Dirichlet"][[i, 2]], m[[slab[n, i, -1]]] = 0.]];
      Flatten[m]], {i, n}];
  ties = Table[If[spec["Periodic"][[i]],
      With[{a = ArrayReshape[Range[Times @@ fs[[i]]], fs[[i]]]}, {Flatten[a[[slab[n, i, 1]]]], Flatten[a[[slab[n, i, -1]]]]}], None], {i, n}];
  tq = First[ndSplineTables[kk, trange[[1]], trange[[2]], qo]];
  q = Length[tq]; w = N[(trange[[2]] - trange[[1]])/kk/qo];
  ns = kk qo;
  vq = nd2WSmoothMach[spec, tq];
  fine = Subdivide[N[trange[[1]]], N[trange[[2]]], 8 ns];
  p2 = polyAt[N[spec["W2"]], fine];
  e1 = If[vr =!= None && spec["C1"] != 0, N[spec["C1"]]/(Ramp[fine - N[vr[[1]]]] + 10.^-300), 0. fine];
  e2 = If[vr =!= None && spec["C2"] != 0, N[spec["C2"]]/(Ramp[N[vr[[2]]] - fine] + 10.^-300), 0. fine];
  segmax = Table[Max[p2[[8 s + 1 ;; 8 s + 9]] + If[s == 0 && vr =!= None, 0., e1[[8 s + 1 ;; 8 s + 9]]] + If[s == ns - 1 && vr =!= None, 0., e2[[8 s + 1 ;; 8 s + 9]]]], {s, 0, ns - 1}];
  w2max = Flatten[Table[Max[segmax[[Select[{kap qo + j, kap qo + j + 1}, kap qo + 1 <= # <= (kap + 1) qo &]]]], {kap, 0, kk - 1}, {j, 0, qo}]];
  endAdd = ConstantArray[0., q];
  If[vr =!= None, endAdd[[1]] = N[spec["C1"]] w; endAdd[[q]] = N[spec["C2"]] w];
  tabM = 0.;
  If[spec["TableMat"] =!= None,
    With[{m = spec["TM"], z = Clip[(tq - N[vr[[1]]])/N[(vr[[2]] - vr[[1]])/spec["TM"]], {0., N[spec["TM"]]}]},
      Module[{j = Clip[Floor[z], {0, m - 1}], th},
        th = z - j;
        tabM = Normal[spec["TableMat"] . Transpose[SparseArray[Join[Transpose[{Range[q], j + 1}], Transpose[{Range[q], j + 2}]] -> Join[1 - th, th], {q, m + 1}]]]]]];
  tilt = Association[# -> If[spec["Av"] === None, ConstantArray[0., nc], -Flatten[spec["Av"][#]]] & /@ verts];
  If[vr === None,
    sa = Flatten[Values[tilt]] + lam0;
    pad = Max[1., 2 Abs[lam0], Max[sa] - Min[sa]];
    tails = {nd2TailFloat[N[spec["P"]], N[trange[[1]]], -1, {Min[sa] - pad, Max[sa] + pad}],
      nd2TailFloat[N[spec["P"]], N[trange[[2]]], +1, {Min[sa] - pad, Max[sa] + pad}]}];
  Join[st, <|"Mask" -> mask, "Ties" -> ties, "TQ" -> Developer`ToPackedArray[tq],
    "VqM" -> Developer`ToPackedArray[ConstantArray[vq - endAdd, nc] + tabM], "W2maxM" -> ConstantArray[Developer`ToPackedArray[N[w2max]], nc],
    "Tilt" -> tilt, "Tails" -> tails, "MeanVol" -> If[spec["Mean"] === None, None, N[spec["Mean"] spec["Volume"]]]|>]];

nd2CalUnpack[st_, x_] := Module[{n = st["n"], cm},
  cm = Table[ArrayReshape[x[[st["Off"][[i]] + 1 ;; st["Off"][[i + 1]]]], {st["NF"][[i]], st["K"] + 2}] st["Mask"][[i]], {i, n}];
  Do[If[st["Ties"][[i]] =!= None, cm[[i, st["Ties"][[i, 2]]]] = cm[[i, st["Ties"][[i, 1]]]]], {i, n}];
  {cm, If[st["MeanVol"] === None, 0., x[[-1]]]}];

(* minus the smoothed functional and its gradient (ndCalValueGradient with the term (lam - a_v) t at the vertices, the
   tails at shifted arguments or no tails, a multiplier, and tied periodic faces) *)
nd2CalValueGradient[st_, x_, tau_] := Module[
  {n = st["n"], kk = st["K"], q = st["Q"], hf = st["h"], k = st["k"], wv = st["wv"], nc = st["NC"], verts = st["Verts"], it = 1./tau, cw = st["cw"],
   tq = st["TQ"], tails = st["Tails"], cm, lam, dc, d2, ds, aa, pos, base, obj = 0., gtl, gth, glam = 0., exs, sv, ex, mx, sm, wgt, tl, dtl, th, dth,
   tlv, thv, etl, eth, ga, gd, g2, gdc, grads, gi, gw},
  {cm, lam} = nd2CalUnpack[st, x];
  dc = Sum[(cm[[i, st["Hi"][[i]]]] - cm[[i, st["Lo"][[i]]]])/hf[[i]], {i, n}];
  d2 = dc . st["B2T"];
  ds = Table[cm[[i]] . st["DBT"], {i, n}];
  aa = Table[With[{sq = ds[[i]]^2 (it/(4 k))}, {sq[[st["Lo"][[i]]]], sq[[st["Hi"][[i]]]]}], {i, n}];
  pos = Ramp[st["W2maxM"] - d2[[All, st["Kq"]]]];
  base = (dc . st["BT"] + cw pos - st["VqM"]) it;
  gtl = ConstantArray[0., nc]; gth = gtl;
  exs = Table[
    sv = lam + st["Tilt"][v];
    ex = base + Sum[aa[[i, v[[i]] + 1]], {i, n}] + Outer[Times, sv, tq] it;
    If[tails === None,
      mx = Max /@ ex;
      ex = Exp[ex - mx]; sm = Total[ex, {2}];
      obj -= wv tau (Total[mx] + Total[Log[sm/q]]);
      wgt = -wv/sm;
      ex = ex wgt;
      glam += Total[ex . tq],
      {tl, dtl} = tails[[1]][sv]; {th, dth} = tails[[2]][sv];
      tlv = ((dc[[All, 1]] + dc[[All, 2]])/2 + tl) it; thv = ((dc[[All, kk + 1]] + dc[[All, kk + 2]])/2 + th) it;
      mx = Max /@ ex; mx = mx + Ramp[tlv - mx]; mx = mx + Ramp[thv - mx];
      ex = Exp[ex - mx]; etl = Exp[tlv - mx]; eth = Exp[thv - mx];
      sm = Total[ex, {2}] + etl + eth;
      obj -= wv tau (Total[mx] + Total[Log[sm/(q + 2)]]);
      wgt = -wv/sm;
      ex = ex wgt;
      gtl += etl wgt; gth += eth wgt;
      glam += Total[ex . tq] + (etl wgt) . dtl + (eth wgt) . dth];
    ex,
    {v, verts}];
  ga = Table[Total[Pick[exs, verts[[All, i]], side]], {i, n}, {side, {0, 1}}];
  gd = ga[[1, 1]] + ga[[1, 2]];
  g2 = -cw Total[ArrayReshape[gd Unitize[pos], {nc, kk, st["Qo"] + 1}], {3}];
  gdc = gd . st["B"] + g2 . st["B2"];
  If[tails =!= None,
    gdc[[All, 1]] += gtl/2; gdc[[All, 2]] += gtl/2; gdc[[All, kk + 1]] += gth/2; gdc[[All, kk + 2]] += gth/2];
  grads = Table[
    gi = ConstantArray[0., {st["NF"][[i]], kk + 2}];
    gi[[st["Hi"][[i]]]] = gdc/hf[[i]]; gi[[st["Lo"][[i]]]] -= gdc/hf[[i]];
    gw = ConstantArray[0., {st["NF"][[i]], q}];
    gw[[st["Lo"][[i]]]] = ga[[i, 1]]; gw[[st["Hi"][[i]]]] += ga[[i, 2]];
    obj += Total[st["FWC"][[i]] cm[[i]], 2];
    gi = (gi + ((gw ds[[i]]) . st["DB"])/(2 k) + st["FWC"][[i]]) st["Mask"][[i]];
    If[st["Ties"][[i]] =!= None,
      gi[[st["Ties"][[i, 1]]]] += gi[[st["Ties"][[i, 2]]]]; gi[[st["Ties"][[i, 2]]]] = 0. gi[[st["Ties"][[i, 2]]]]];
    Flatten[gi],
    {i, n}];
  If[st["MeanVol"] =!= None, obj += lam st["MeanVol"]; glam += st["MeanVol"]];
  {-obj, -Join[Flatten[grads], If[st["MeanVol"] === None, {}, {glam}]]}];

nd2CalStage[spec_, st_, field_, tau_, iters_] := Module[{x, res, cm, lam},
  x = Join[Flatten[Flatten /@ field["C"]], If[st["MeanVol"] === None, {}, {N[field["Lambda"]]}]];
  res = Quiet[ndLBFGS[nd2CalValueGradient[st, #, tau] &, x, iters], {General::munfl}];
  {cm, lam} = nd2CalUnpack[st, res[[1]]];
  {Join[field, <|"C" -> Table[Developer`ToPackedArray[ArrayReshape[cm[[i]], Append[faceShape[spec, i], st["K"] + 2]]], {i, spec["n"]}],
     "Lambda" -> lam|>], -res[[2]], res[[3]]}];

nd2ProlongX[f_, specFine_] := Append[
  ndProlongX[f, Append[specFine, "Dirichlet" -> MapThread[If[#2, {True, True}, #1] &, {specFine["Dirichlet"], specFine["Periodic"]}]]],
  "Lambda" -> f["Lambda"]];

(* the model on a coarser grid (coarse levels of calibration mode): nodal arrays are restricted, cell arrays averaged *)
nd2ParamsOnGrid[params_, cells_List, cellsFine_List] := Module[{p, r = cellsFine/cells, n = Length[cells], src, tab, block},
  If[cells === cellsFine, Return[params]];
  p = ndParamsOnGrid[params, cells, cellsFine];
  block[a_] := BlockMap[Mean[Flatten[#]] &, a, r];
  src = Lookup[params, "Source", None];
  Which[
    ArrayQ[src, n] && Dimensions[src] === cellsFine + 1, p["Source"] = src[[Sequence @@ Table[Span[1, All, r[[j]]], {j, n}]]],
    ArrayQ[src, n] && Dimensions[src] === cellsFine, p["Source"] = block[src]];
  tab = Lookup[params, "DataTerm", None];
  If[ArrayQ[tab, n + 1] && Most[Dimensions[tab]] === cellsFine,
    p["DataTerm"] = Transpose[Table[block[tab[[Sequence @@ ConstantArray[All, n], j]]], {j, Last[Dimensions[tab]]}], Prepend[Range[n], n + 1]]];
  p];

nd2Calibrate[params_, levels_List, trange_, lam0_, qo_, qs_, scale_, sched_, monitor_] := Module[
  {cellsFine = levels[[-1, 1]], spec, prev = None, field = None, best, bestLB, lb, st, stages, rows = {}, t0 = AbsoluteTime[], res},
  Do[
    With[{cells = levels[[li, 1]], kk = levels[[li, 2]]},
      spec = nd2Spec[nd2ParamsOnGrid[params, cells, cellsFine]];
      If[field === None,
        field = Append[ndZeroField[spec, trange, kk], "Lambda" -> lam0],
        If[cells =!= prev["Cells"], field = nd2ProlongX[field, spec]];
        If[kk =!= field["K"], field = Append[ndRefineT[field], "Lambda" -> field["Lambda"]]]];
      best = field; bestLB = nd2LowerBound[spec, field, "Subintervals" -> qs];
      st = nd2CalSetup[spec, trange, kk, qo, field["Lambda"]];
      stages = Which[li == 1, sched["First"], li == Length[levels], sched["Last"], True, sched["Later"]];
      Do[
        res = nd2CalStage[spec, st, field, stage[[1]] scale, stage[[2]]];
        field = res[[1]];
        lb = nd2LowerBound[spec, field, "Subintervals" -> qs];
        If[lb > bestLB, bestLB = lb; best = field];
        If[monitor =!= None, monitor[<|"Level" -> li, "Levels" -> Length[levels], "Cells" -> cells, "K" -> kk, "Tau" -> stage[[1]] scale,
          "Iterations" -> res[[3]], "Smoothed" -> res[[2]], "Certified" -> N[lb], "Best" -> N[bestLB], "Seconds" -> AbsoluteTime[] - t0|>]],
        {stage, stages}];
      field = best; prev = spec;
      AppendTo[rows, <|"Cells" -> cells, "K" -> kk, "LowerBound" -> N[bestLB], "Seconds" -> AbsoluteTime[] - t0|>]],
    {li, Length[levels]}];
  {spec, field, rows, bestLB}];

(* the t-range without a value range: the hull of the Dirichlet data, the mean and the critical points of P(t) + a t for
   the extreme values of a, rounded outward to a binary fraction of the span (heuristic) *)
nd2AutoTRange[spec_] := Module[{t, as, roots, vals, lo, hi, span, qn},
  as = Union[{0}, If[spec["Lin"] === None, {}, -SetPrecision[{Min[spec["Lin"]], Max[spec["Lin"]]}, Infinity]]];
  roots = Flatten[Table[Quiet@Check[t /. NSolve[polyAt[spec["W1"], t] + a == 0, t, Reals], {}], {a, as}]];
  roots = Rationalize[Select[roots, NumericQ[#] && Im[#] == 0 &], 0];
  vals = Join[SetPrecision[Pick[Flatten[spec["G"]], Flatten[spec["Fixed"]], 1], Infinity], roots, If[spec["Mean"] === None, {}, {spec["Mean"]}]];
  If[vals === {}, Return[{-1, 1}]];
  lo = Min[vals]; hi = Max[vals]; span = N[hi - lo];
  If[! (span > 0.), span = 1.; lo = lo - 1/2; hi = hi + 1/2];
  qn = 2^Floor[Log2[span/32]];
  {Floor[lo/qn] qn, Ceiling[hi/qn] qn}];

nd2AutoScale[spec_, {t0_, t1_}] := Module[{ts = Subdivide[N[t0], N[t1], 256], w, tb = spec["TableMat"]},
  w = nd2WSmoothMach[spec, ts];
  Max[1., Max[w] - Min[w] + If[tb === None, 0., Max[Max[#] - Min[#] & /@ tb]]]];

(* the multiplier from the candidate: the mean over the cells of W'(u) + a at the cell means (heuristic) *)
nd2Lambda0[spec_, u_] := Module[{n = spec["n"], nc = Times @@ spec["Cells"], ub, a},
  ub = Flatten[Mean[vertexPart[u, n, #] & /@ Tuples[{0, 1}, n]]];
  a = If[spec["Av"] === None, 0., Flatten[Mean[Values[spec["Av"]]]]];
  Mean[nd2WMach[spec, ub, Range[nc]][[2]] + a]];

(* ------------------------------------------------------------------------------------------------ *)
(* the certificate                                                                                   *)
(* ------------------------------------------------------------------------------------------------ *)

certifyScalarFieldND2[params_Association, cand_, opts_List] := Module[
  {spec, t0 = AbsoluteTime[], mname, sub, verbose, refine, gfun, given, mode, u, ub, fields = {}, lbs, best, lo, hi, margin, qs = 16,
   trange, kfin, levels, scale, cal, rows = {}, known = <||>, lam0 = 0., vr},
  spec = nd2Spec[KeyDrop[params, "Field"]];
  vr = spec["VR"];
  {mname, sub} = ndMethod[OptionValue[CertifiedVariationalBound, opts, Method]];
  verbose = TrueQ[OptionValue[CertifiedVariationalBound, opts, "Verbose"]];
  refine = TrueQ[OptionValue[CertifiedVariationalBound, opts, "RefineCandidate"]];
  If[mname === "Convex",
    Throw[failure["NotImplemented", "Convex mode is not available with \"Source\", \"Mean\", \"PeriodicDirections\", \"ValueRange\" or \"DataTerm\"; use Method -> Automatic or \"Calibration\"."], "CertiField"]];
  If[spec["Linear"] && mname === "Calibration",
    Throw[failure["BadOption", "A potential that is at most linear in u is certified with an equilibrated flux: use Method -> Automatic."], "CertiField"]];
  If[! spec["Linear"] && vr === None && ! superlinearQ[spec["P"]],
    Throw[failure["Inconclusive", "The potential must be at most linear in u, have even degree with a positive leading coefficient, or have a value range."], "CertiField"]];
  given = Lookup[params, "Field", None];
  If[given =!= None,
    given = nd2PackField[spec, given];
    If[given === $Failed, Throw[failure["BadField", "\"Field\" must be a field on the grid of the model."], "CertiField"]]];
  mode = If[spec["Linear"], "Equilibrated", "Calibration"];
  gfun = With[{g = Lookup[params, "BoundaryData", 0]}, If[NumericQ[g] || ArrayQ[g, spec["n"]], None, g]];
  If[verbose, Print["  ScalarFieldND: ", ndGridString[spec["Cells"]], " cells, ", ToLowerCase[mode], " field"]];
  u = nd2Competitor[spec, cand, refine, gfun];
  ub = nd2UpperBound[spec, u];
  If[verbose, Print["  candidate: E <= ", ndNum[ub[[2]]], "   (", ndSec[AbsoluteTime[] - t0], " s)"]];
  If[AssociationQ[given], AppendTo[fields, given]];
  If[mode === "Equilibrated", AppendTo[fields, nd2MixedField[spec]]];
  If[mode === "Calibration" && (! AssociationQ[given] || mname === "Calibration"),
    qs = sub["Subintervals"];
    If[vr =!= None,
      If[sub["TRange"] =!= Automatic && exactify[sub["TRange"]] =!= vr,
        Throw[failure["BadOption", "With a value range the t-range is the value range."], "CertiField"]];
      trange = vr,
      trange = exactify[Replace[sub["TRange"], Automatic :> nd2AutoTRange[spec]]]];
    If[! (MatchQ[trange, {_Integer | _Rational, _Integer | _Rational}] && trange[[1]] < trange[[2]]),
      Throw[failure["BadOption", "\"TRange\" must be {tlo, thi} with rational tlo < thi."], "CertiField"]];
    kfin = Replace[sub["KnotIntervals"], Automatic -> 16];
    levels = Replace[sub["Levels"], Automatic :> If[IntegerQ[kfin] && kfin >= 1, ndAutoLevels[spec["Cells"], kfin], $Failed]];
    If[! ndLevelsQ[levels, spec["n"]] || levels[[-1, 1]] =!= spec["Cells"],
      Throw[failure["BadOption", "\"Levels\" must be a list of {cells, K} that doubles the cells and/or K from one level to the next and ends on the grid of the model; \"KnotIntervals\" must be a positive integer."], "CertiField"]];
    If[spec["TableMat"] =!= None && ! AllTrue[levels, Divisible[#[[2]] qs, spec["TM"]] &],
      Throw[failure["BadOption", "The table nodes of the data term must be subinterval ends: M must divide K times \"Subintervals\" on every level."], "CertiField"]];
    scale = Replace[sub["Scale"], Automatic :> nd2AutoScale[spec, trange]];
    If[! (NumericQ[scale] && scale > 0), Throw[failure["BadOption", "\"Scale\" must be a positive number."], "CertiField"]];
    lam0 = If[spec["Mean"] === None, 0., nd2Lambda0[spec, u]];
    If[verbose, Print["  t-range ", trange, ", scale ", scale, ", levels ", StringRiffle[(ndGridString[#[[1]]] <> " K" <> ToString[#[[2]]]) & /@ levels, " -> "]]];
    cal = nd2Calibrate[KeyDrop[params, "Field"], levels, trange, lam0, sub["Samples"], qs, scale,
      If[AssociationQ[sub["Schedule"]], Join[$ndSchedule, sub["Schedule"]], $ndSchedule],
      If[verbose, Function[s, Print["  level ", s["Level"], "/", s["Levels"], "  ", ndGridString[s["Cells"]], "  K ", s["K"], "  tau ", ToString[s["Tau"], InputForm],
        "  iterations ", s["Iterations"], "  certified lower bound ", ndNum[s["Best"]], "   (", ndSec[s["Seconds"]], " s)"]], None]];
    rows = cal[[3]];
    AppendTo[fields, cal[[2]]];
    known[Length[fields]] = cal[[4]]];
  If[fields === {}, Throw[failure["Inconclusive", "No lower bound is available."], "CertiField"]];
  lbs = Table[Lookup[known, j, nd2LowerBound[spec, fields[[j]], "Subintervals" -> qs]], {j, Length[fields]}];
  best = First[Ordering[N[lbs], -1]];
  lo = exactNumber[lbs[[best]]]; hi = exactNumber[ub[[2]]];
  margin = Max[Abs[lo], Abs[hi]]/10^9;
  <|"LowerBound" -> lo - margin, "UpperBound" -> hi + margin, "Seconds" -> AbsoluteTime[] - t0, "Field" -> fields[[best]],
    "Nodal" -> u, "BoundaryNodalValues" -> spec["G"], "SourceValues" -> spec["Lin"], "Subintervals" -> qs,
    "Mode" -> If[fields[[best]]["K"] == 0, "Equilibrated", "Calibration"], "Levels" -> rows, "Spec" -> spec|>];

(* re-verification from the stored nodal values, boundary data, source values and field *)
verifyScalarFieldND2[a_Association] := Module[{c = a["Certificate"], spec, g, s, n, lb, ub},
  spec = nd2Spec[KeyDrop[a["Parameters"], "Field"]];
  n = spec["n"];
  g = ndRealArray[Lookup[c, "BoundaryNodalValues", $Failed], spec["Cells"] + 1];
  If[g === $Failed || Max[Abs[(1 - spec["Fixed"]) g]] != 0. || Max[Abs[g - spec["G"]]] > 10.^-12 (1. + Max[Abs[spec["G"]]]) || ! nd2TiedQ[spec, g],
    Return[False]];
  spec = Append[spec, "G" -> g];
  If[spec["Lin"] =!= None,
    s = ndRealArray[Lookup[c, "SourceValues", $Failed], Dimensions[spec["Lin"]]];
    If[s === $Failed || Max[Abs[s - spec["Lin"]]] > 10.^-12 (1. + Max[Abs[spec["Lin"]]]), Return[False]];
    spec = Join[spec, <|"Lin" -> s, "Av" -> Association[# -> If[spec["LinNodal"], -vertexPart[s, n, #], -s] & /@ Tuples[{0, 1}, n]]|>]];
  lb = nd2LowerBound[spec, c["Field"], "Subintervals" -> Lookup[c, "Subintervals", 16]];
  ub = nd2UpperBound[spec, c["Nodal"]][[2]];
  NumericQ[lb] && NumericQ[ub] && lb >= exactNumber[a["LowerBound"]] && ub <= exactNumber[a["UpperBound"]]];

(* the potential as text *)
nd2PolyString[p_List] := Module[{terms},
  terms = DeleteCases[MapIndexed[If[#1 == 0, Nothing, Which[#2[[1]] == 1, exactString[#1], #2[[1]] == 2, exactString[#1] <> " t",
    True, exactString[#1] <> " t^" <> ToString[#2[[1]] - 1]]] &, p], Nothing];
  If[terms === {}, "0", StringReplace[StringRiffle[terms, " + "], "+ -" -> "- "]]];

nd2Description[params_, res_] := Module[{sp = res["Spec"], d = Max[displayDigits[res], 8], grid, nat, per, wtext, vr},
  vr = sp["VR"];
  grid = ndGridString[sp["Cells"]];
  nat = Flatten[Table[If[sp["Dirichlet"][[i, s]] || sp["Periodic"][[i]], Nothing, {i, If[s == 1, "Lower", "Upper"]}], {i, sp["n"]}, {s, 2}], 1];
  per = Flatten[Position[sp["Periodic"], True]];
  wtext = StringJoin[nd2PolyString[sp["P"]],
    If[vr =!= None && sp["C1"] != 0, " + " <> exactString[sp["C1"]] <> " (t - " <> exactString[vr[[1]]] <> ") Log(t - " <> exactString[vr[[1]]] <> ")", ""],
    If[vr =!= None && sp["C2"] != 0, " + " <> exactString[sp["C2"]] <> " (" <> exactString[vr[[2]]] <> " - t) Log(" <> exactString[vr[[2]]] <> " - t)", ""],
    If[sp["TableMat"] =!= None, " + D_c(t), D_c the piecewise linear interpolant in t, in each cell c, of the tabulated data term on the uniform grid of " <>
      ToString[sp["TM"]] <> " intervals of the value range", ""],
    If[sp["Lin"] =!= None, " - f(x) t, f " <> If[sp["LinNodal"], "the piecewise multilinear interpolant of the nodal source values", "constant on each cell"] <>
      " (machine numbers, stored in the certificate)", ""],
    If[vr =!= None, ", and W = +Infinity for t outside " <> exactString[vr], ""]];
  <|"Units" -> "energy units of the model",
    "Theorem" -> Row[{"Every u in H^1 of the box ", exactString[sp["Domain"]], If[per === {}, "", ", periodic in the directions " <> exactString[per]],
      ", whose trace on the Dirichlet faces is the piecewise multilinear interpolant of the nodal boundary values on the ", grid, " grid",
      If[sp["Mean"] === None, "", " and whose mean is " <> exactString[sp["Mean"]]],
      ", has E[u] \[GreaterEqual] ", lowerString[res["LowerBound"], d], ", where E[u] = \[Integral] (", exactString[sp["k"]], " |\[Del]u|^2 + W(x, u)) dx and W(x, t) = ",
      wtext, "; the competitor has E \[LessEqual] ", upperString[res["UpperBound"], d], "."}],
    "AdmissibleClass" -> Row[{"all u in H^1 of the box ", exactString[sp["Domain"]], If[per === {}, "", ", periodic in the directions " <> exactString[per]],
      ", whose trace on the Dirichlet faces is the piecewise multilinear interpolant of the nodal boundary values (machine numbers, stored in the certificate) on the ",
      grid, " grid", If[sp["Mean"] === None, "", "; mean " <> exactString[sp["Mean"]]], If[nat === {}, "", Row[{"; natural (free) faces: ", nat}]]}],
    "LowerBoundMethod" -> If[res["Mode"] === "Equilibrated",
      "equilibrated Raviart\[Dash]Thomas flux and the Poincar\[EAcute] inequality on boxes (docs/SCALAR_ND_2.md, Proposition 16); outward-rounded arithmetic",
      Row[{"calibration field (Raviart\[Dash]Thomas in x, quadratic B-splines in u, ", res["Field"]["K"], " knot intervals on ", exactString[res["Field"]["TRange"]],
        If[sp["Mean"] === None, "", ", a multiplier for the mean"], "); outward-rounded arithmetic"}]],
    "UpperBoundMethod" -> "energy of the piecewise multilinear function through the candidate's nodal values" <>
      If[sp["Mean"] === None, "", ", shifted to the exact mean"] <> " (exact gradient term, midpoint and vertex rules with remainders)",
    "TrustedBase" -> "IEEE-754 machine arithmetic (+, \[Minus], \[Times], \[Divide]), Wolfram Language exact rational arithmetic, and docs/SCALAR_ND.md and docs/SCALAR_ND_2.md (proofs checked by one author only); logarithms are enclosed from series with remainder bounds",
    "NotClaimed" -> "that the candidate is the minimiser or that the minimiser is unique; the gap bounds how far below the candidate any admissible field can be",
    "References" -> "docs/SCALAR_ND.md, docs/SCALAR_ND_2.md; calibrations and duality: Alberti\[Dash]Bouchitt\[EAcute]\[Dash]Dal Maso (2003), Bouchitt\[EAcute]\[Dash]Fragal\[AGrave] (2018); complementary energy: Prager\[Dash]Synge (1947)."|>];

End[];
