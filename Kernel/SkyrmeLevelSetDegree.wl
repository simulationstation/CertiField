(* ::Package:: *)

(* ::Title:: *)
(*CertiField: the level-set lower bound for the Skyrme energy in degree B*)

(* ::Text:: *)
(*Theorem (docs/SKYRME_DEGREE_TWO.md, Theorem 2).  For every finite-energy U: R^3 -> S^3 of degree B or -B, B >= 1,*)
(*    E(U) >= W_B(tau_1, v)  for some v >= 0,*)
(*where W_B is the value function of the relaxed level-set problem in which the flux through a level surface is*)
(*B c(tau) = 4 pi B sin^2 tau.  The degree of the map enters the analytic reduction in this one place.*)
(**)
(*levelSetLowerBoundDegree below is levelSetLowerBound of SkyrmeLevelSet.wl with that change and nothing else: the*)
(*three flux integrals of a step (C0, Ca, Cb) carry the factor B, and so do the two floating-point numbers ca, cb that*)
(*only steer heuristic choices.  The floors Xi and Xi_mu do not contain the flux and are unchanged; the Bogomolny*)
(*floor is accumulated from C0 and so becomes 6 B int c.  For B = 1 the two functions perform the same operations*)
(*and return the same numbers (Tests/SkyrmeDegree.wlt).  The body is kept textually parallel to the original so that*)
(*the two can be compared line by line.*)
(**)
(*Rigour: as in SkyrmeLevelSet.wl.  Grids are heuristic (Reference/SkyrmeDegree/make_grids.py) and do not affect*)
(*validity.*)

Begin["CertiField`Private`"];

(* ---- the verifier, with the flux of degree deg ---- *)
Options[levelSetLowerBoundDegree] = {"m1" -> 8., "m2" -> 3., "nmin" -> 6, "Monitor" -> None, "Trace" -> False};

levelSetLowerBoundDegree[deg_Integer?Positive, muq_, grid_Association, tab_Association, OptionsPattern[]] := Module[
  {cst = skyrmeConstants[], m1 = N@OptionValue["m1"], m2 = N@OptionValue["m2"], nmin = OptionValue["nmin"],
   monitor = OptionValue["Monitor"], mu2q = muq^2, muLo = fLo[muq], taus = grid["taus"], vref = grid["vref"],
   eug = grid["eug"], offs = grid["offs"], nLev, tq, S1, C1, S2, C2, CH, XiC, XmC, Vlo, Vhi, nodes, xiLo, ximLo,
   vb, x, xm, dvb, useM, YL, slp, Bog, T, PRE, nfar = 0, nfallback = 0, npairs = 0, t0 = AbsoluteTime[],
   ta, tb, dt, dq, halfDt, s2a, s2b, c2a, c2b, b0, bbq, baq, C0, Cb, Ca, VbLo, VbHi, VbIv, mterm, ca, cb, va, vL, vR,
   J, Jb, gL, gR, vmid, dA, wL1, wR1, end1, F1, suf, klo, shift, vtop, khi, cnt, starts, P, jj, kk, wlo, whi, valid,
   Yk, Sk, w0, wmid, PREk, lift, dB, dAj, cm, scale, dUp, KUp, KDn, Kt, alB, beB, khB, gam, al1, dk, alA, beA, khA,
   rhs, good, KIv, hb, alsum, BL, BR, alF, beF, khF, beFar, farStep, far1, suf2, base2, far6, f2L, f2R, use2, farL,
   farR, xa, xma, dvIv, dv, candY, candS, mid, best, idx, msk, YLn, Sn, endn, cmax, pre, ktop, fin, val, cert,
   endBound, vend, gend, trace = TrueQ[OptionValue["Trace"]], i0 = Lookup[grid["Meta"], "i0", None], traceRows = {}, kr},

  nLev = Length[taus] - 1;
  (* exact rational data at the grid points *)
  tq = exactOf /@ taus;
  If[! (tq[[-1]] < $piLo && tq[[1]] > 0 && OrderedQ[tq] && Length[vref] == nLev + 1 && Length[offs] == nLev + 1),
    Throw[failure["BadGrid", "The tau grid must be increasing and lie in (0, pi)."], "CertiField"]];
  S1 = sinBr /@ tq; C1 = cosBr /@ tq; S2 = sinBr /@ (2 tq); C2 = cosBr /@ (2 tq); CH = cosBr /@ (tq/2);
  XiC = fLo[4 $piLo ($piLo - #)^2] & /@ tq;                        (* <= 4 pi (pi - tau)^2 *)
  XmC = fLo[16 $piLo Max[#[[1]], 0]^2] & /@ CH;                     (* <= 16 pi cos^2(tau/2) *)
  Vlo = fLo[2 mu2q (1 + #[[1]])] & /@ C1;
  Vhi = fHi[2 mu2q (1 + #[[2]])] & /@ C1;

  nodes[i_] := Join[{0.}, vref[[i + 1]] eug[[offs[[i + 1]] + 1 ;;]]];         (* i: 0-based level *)
  xiLo[i_, v_] := Ramp[rDn[XiC[[i + 1]] radLoV[v, cst]]];
  ximLo[i_, v_] := With[{r = radLoV[v, cst]}, Ramp[rDn[XmC[[i + 1]] rDn[r + rDn[muLo rDn[r r]]]]]];

  (* terminal data: chords of Xi and Xi_mu on each cell *)
  vb = nodes[nLev];
  x = xiLo[nLev, vb]; xm = ximLo[nLev, vb];
  dvb = ivSub[ivPt[Rest[vb]], ivPt[Most[vb]]];
  useM = gtMask[Most[xm] + Rest[xm], Most[x] + Rest[x]];
  YL = vWhere[useM, Most[xm], Most[x]];
  slp = vWhere[useM, slopeLo[Most[xm], Rest[xm], dvb], slopeLo[Most[x], Rest[x], dvb]];
  Bog = 0.;
  T = Max[Last[x], Last[xm]];
  PRE = ConstantArray[-$big, Length[YL]];

  Do[
    ta = taus[[i + 1]]; tb = taus[[i + 2]]; dt = tb - ta;
    dq = tq[[i + 2]] - tq[[i + 1]];
    If[! (exactOf[dt] == dq && dq > 0), Throw[failure["BadGrid", "Step widths must be exact doubles."], "CertiField"]];
    halfDt = 0.5 dt;
    {s2a, s2b, c2a, c2b} = {S2[[i + 1]], S2[[i + 2]], C2[[i + 1]], C2[[i + 2]]};
    (* exact lower bounds of  int c,  int (1-t) c,  int t c  over the step, c = 4 pi deg sin^2 *)
    b0 = dq - (s2b[[2]] - s2a[[1]])/2;
    bbq = dq^2/2 - dq s2b[[2]]/2 - (c2b[[2]] - c2a[[1]])/4;
    baq = dq^2/2 + dq s2a[[1]]/2 + (c2b[[1]] - c2a[[2]])/4;
    If[! (b0 > 0), Throw[failure["Internal", "b0 <= 0"], "CertiField"]];
    C0 = fLo[2 deg $piLo b0]; Cb = fLo[2 deg $piLo Max[bbq, 0]/dq]; Ca = fLo[2 deg $piLo Max[baq, 0]/dq];
    VbLo = Vlo[[i + 2]]; VbHi = Vhi[[i + 2]]; VbIv = {VbLo, VbHi};
    mterm = fLo[2 mu2q Min[S1[[i + 1, 1]], S1[[i + 2, 1]]]];
    ca = deg $fourPi Sin[ta]^2; cb = deg $fourPi Sin[tb]^2;                  (* heuristics only *)
    va = nodes[i]; vL = Most[va]; vR = Rest[va];
    J = Length[vL]; Jb = Length[vb] - 1;
    gL = gLoV[vL, cst]; gR = gLoV[vR, cst];
    vmid = vWhere[gtMask[vL, 0.], Sqrt[vL vR], 0.5 vR];
    dA = tabDIndex[tab, ca/gFun[vmid]];
    wL1 = Most[vb]; wR1 = Rest[vb];
    end1 = First[ivAdd[ivPt[YL], ivMul[ivPt[slp], ivSub[ivPt[wR1], ivPt[wL1]]]]];
    F1 = vMax[vMin[YL, end1], PRE];
    suf = Append[suffixMin[F1], $big];
    klo = Clip[searchRight[vb, vL] - 1, {0, Jb - 1}];
    shift = Log[vref[[i + 2]]/vref[[i + 1]]];
    vtop = vMax[vR Exp[m1 shift], vR + m2 Max[ca, cb] dt];
    khi = Clip[searchRight[vb, vtop], {0, Jb - 1}];
    khi = khi + Ramp[Clip[klo + nmin, {0, Jb - 1}] - khi];
    cnt = khi - klo + 1;
    starts = Most[Prepend[Accumulate[cnt], 0]];
    P = Total[cnt];
    jj = repeatIndex[cnt];
    kk = klo[[jj + 1]] + (Range[0, P - 1] - starts[[jj + 1]]);
    wlo = vMax[wL1[[kk + 1]], vL[[jj + 1]]]; whi = wR1[[kk + 1]];
    valid = gtMask[whi, wlo];
    Yk = YL[[kk + 1]]; Sk = slp[[kk + 1]]; w0 = wL1[[kk + 1]];
    wmid = vWhere[gtMask[wlo, 0.], Sqrt[Abs[wlo whi]], 0.5 whi];
    PREk = PRE[[kk + 1]];
    lift = gtMask[PREk, Yk + Sk (0.5 (wlo + whi) - w0)];
    Yk = vWhere[lift, PREk, Yk]; Sk = vWhere[lift, 0., Sk];
    dB = tabDIndex[tab, cb/gFun[wmid]]; dAj = dA[[jj + 1]];
    (* heuristic choice of the multiplier (floating point) *)
    cm = 0.5 (ca + cb);
    scale = vWhere[gtMask[dAj, 0], tab["DG"][[dAj + 1]]/(cm dt), 1./vMax[gFun[vmid[[jj + 1]]] dt, 1.*^-300]];
    dUp = wlo - vR[[jj + 1]];
    KUp = vWhere[gtMask[dUp, 0.], tabKInv[tab, Ramp[dUp] scale, dAj], $big];
    KDn = tabKInv[tab, (whi - vL[[jj + 1]]) scale, dAj];
    Kt = vMin[vMax[vMax[Sk, 0.], KDn - VbHi], vMax[KUp - VbHi, 0.]];
    (* rigorous minorants at both ends of the step *)
    {alB, beB, khB} = tabMinorIv[tab, VbHi + Kt, dB];
    gam = gPrimeLoV[whi, cst];
    al1 = First[tabMinor[tab, khB[[2]] + dt (gam alB + mterm), dAj]];
    dk = 0.999 dt (gam vMin[al1, alB] + mterm);
    {alA, beA, khA} = tabMinorIv[tab, khB[[2]] + dk, dAj];
    rhs = rDn[dt rDn[rDn[gam vMin[alA, alB]] + mterm]];
    good = leMask[rUp[khA[[2]] - khB[[1]]], rhs];                            (* condition [*] *)
    alA = vWhere[good, alA, alB]; beA = vWhere[good, beA, beB];             (* fallback: m_a = m_b *)
    KIv = ivSub[khB, VbIv];
    hb = rDn[rDn[beA Ca] + rDn[beB Cb]];
    alsum = rDn[alA + alB];
    endBound[vEnd_, gEnd_] := Module[{ve = vEnd[[jj + 1]], hp, bst = ConstantArray[$big, P], lin, kd},
      hp = rDn[rDn[rDn[gEnd[[jj + 1]] halfDt] alsum] + hb];
      Do[
        lin = First[ivAdd[ivPt[Yk], ivMul[ivPt[Sk], ivSub[ivPt[w], ivPt[w0]]]]];
        kd = Last[ivMul[KIv, ivSub[ivPt[w], ivPt[ve]]]];
        bst = vMin[bst, rDn[rDn[hp + lin] - kd]],
        {w, {wlo, whi}}];
      segmentMin[vWhere[valid, bst, $big], cnt]];
    BL = endBound[vL, gL]; BR = endBound[vR, gR];
    npairs += P; nfallback += Total[valid (1 - good)];
    (* far targets: cells above the window and volumes above the last node *)
    {alF, beF, khF} = tabMinorIv[tab, {VbLo (1 - 1.*^-9)}, {tab["nD"] - 1}];
    beFar = If[khF[[2, 1]] <= VbLo, beF[[1]], 6.];                          (* (0, 6, 0) is always a minorant *)
    farStep = rDn[beFar C0];
    far1 = rDn[farStep + vMin[suf[[khi + 2]], T]];
    suf2 = Append[suffixMin[rDn[F1 + rDn[VbLo wL1]]], $big];
    base2 = vMin[suf2[[khi + 2]], rDn[T + rDn[VbLo Last[vb]]]];
    far6 = rDn[6. C0];
    f2L = rDn[rDn[far6 + base2] - rUp[VbLo vL]]; f2R = rDn[rDn[far6 + base2] - rUp[VbLo vR]];
    use2 = gtMask[f2L + f2R, 2 far1];
    farL = vWhere[use2, f2L, far1]; farR = vWhere[use2, f2R, far1];
    nfar += Total[ltMask[farL, BL]];
    BL = vMin[BL, farL]; BR = vMin[BR, farR];
    If[! (Max[Abs[BL]] < $big && Max[Abs[BR]] < $big), Throw[failure["Internal", "non-finite cell bound"], "CertiField"]];
    (* new level: per cell the best of DP chord, Xi chord, Bogomolny constant, Xi_mu chord *)
    Bog = rDn[Bog + rDn[6. C0]];
    xa = xiLo[i, va]; xma = ximLo[i, va];
    dvIv = ivSub[ivPt[vR], ivPt[vL]]; dv = vR - vL;
    If[Min[dvIv[[1]]] <= 0., Throw[failure["BadGrid", "degenerate volume cell"], "CertiField"]];
    candY = {BL, Most[xa], ConstantArray[Bog, J], Most[xma]};
    candS = {slopeLo[BL, BR, dvIv], slopeLo[Most[xa], Rest[xa], dvIv], ConstantArray[0., J], slopeLo[Most[xma], Rest[xma], dvIv]};
    mid = candY + 0.5 (# dv & /@ candS);
    best = mid[[1]]; idx = ConstantArray[0, J];
    Do[msk = gtMask[mid[[k]], best]; idx = vWhere[msk, k - 1, idx]; best = vWhere[msk, mid[[k]], best], {k, 2, 4}];
    YLn = Sum[(1 - Unitize[idx - (k - 1)]) candY[[k]], {k, 4}];
    Sn = Sum[(1 - Unitize[idx - (k - 1)]) candS[[k]], {k, 4}];
    endn = First[ivAdd[ivPt[YLn], ivMul[ivPt[Sn], dvIv]]];
    cmax = vMax[YLn, endn];
    pre = Prepend[Most[cumMax[cmax]], -$big];
    ktop = Clip[First[searchRight[vb, {Last[va]}]] - 1, {0, Jb}];
    T = Max[Last[xa], Last[xma], Bog, rDn[farStep + Min[suf[[ktop + 1]], T]], Max[cmax]];
    {YL, slp, vb, PRE} = {YLn, Sn, va, pre};
    (* optional record: the two one-sided values of the new function at the reference node (the node u = 0) *)
    If[trace && IntegerQ[i0] && (kr = i0 - offs[[i + 1]] + 2; 2 <= kr <= J),
      AppendTo[traceRows, {i, ta, vb[[kr]], endn[[kr - 1]], YL[[kr]]}]];
    If[monitor =!= None, monitor[<|"Level" -> i, "Levels" -> nLev, "Tau" -> ta, "Pairs" -> npairs,
      "Lambda0" -> YL[[1]]/(12 Pi^2), "Seconds" -> AbsoluteTime[] - t0|>]],
    {i, nLev - 1, 0, -1}];

  fin = cumMax[vMax[vMax[vMax[vMin[BL, BR], Most[xa]], Most[xma]], Bog]];
  val = Min[Min[fin], T];
  cert = rDn[val/cst["TWELVEPI2HI"]];
  <|"LowerBound" -> cert, "Degree" -> deg, "Levels" -> nLev, "VolumeNodes" -> Length[eug], "Pairs" -> npairs,
    "Fallbacks" -> nfallback, "FarBinding" -> nfar, "Seconds" -> AbsoluteTime[] - t0,
    If[trace, "Trace" -> traceRows, Nothing]|>];

(* ---- the model SkyrmeDegreeTwo ---- *)
(* Lower bound: the programme above with B = 2.  Upper bound: the trial map of SkyrmeDegreeOne is the vacuum outside a
   ball, so two translated copies with disjoint supports form one map of degree two whose energy is twice the energy
   of the trial; no inequality between the infima I_1 and I_2 is used. *)

$skyrmeDegreePresets = <|
  "Quick" -> <|"Table" -> {2000, 900, 1/10^10}|>,
  "Standard" -> <|"Table" -> {2000, 900, 1/10^10}|>|>;
$skyrmeDegreeMuNames = <|0 -> "0", 1 -> "1"|>;

CertiFieldModelInfo["SkyrmeDegreeTwo"] := <|
  "Name" -> "SkyrmeDegreeTwo",
  "Energy" -> "E(U) = \[Integral] |dU|^2 + |\[CapitalLambda]^2 dU|^2 + 2\[Mu]^2 (1 - cos f) over \[DoubleStruckCapitalR]^3, f = dist(U, 1); reported as E/(12\[Pi]^2)",
  "Parameters" -> <|"PionMass" -> "\[Mu] \[Element] {0, 1} (default 0)"|>,
  "AdmissibleClass" -> "finite-energy maps U: \[DoubleStruckCapitalR]^3 \[RightArrow] S^3 of degree 2 or \[Minus]2 with U - 1 \[Element] L^6",
  "LowerBound" -> "level-set lower bound with the flux of degree two (docs/SKYRME_DEGREE_TWO.md); certified dynamic programme",
  "UpperBound" -> "twice the certified energy of a hedgehog trial map that is the vacuum outside a ball (two disjoint copies form a map of degree two)",
  "Resolutions" -> Keys[$skyrmeDegreePresets],
  "NotClaimed" -> "that the bounds are close to the infimum: the interval is about 10% wide",
  "Status" -> "Supported"|>;

skyrmeDegreeGridDir[deg_Integer, muq_, preset_String] := FileNameJoin[{$CertiFieldDirectory, "Data", "SkyrmeDegree",
  "b" <> ToString[deg] <> "_mu" <> $skyrmeDegreeMuNames[muq] <> "_" <> ToLowerCase[preset]}];

skyrmeDegreeMu[params_] := Module[{mu = Lookup[params, "PionMass", 0]},
  mu = If[InexactNumberQ[mu], Rationalize[mu, 10^-12], mu];
  If[! KeyExistsQ[$skyrmeDegreeMuNames, mu],
    Throw[failure["BadModel", "Supported pion masses for degree two: 0 and 1 (Battye\[Dash]Sutcliffe units)."], "CertiField"]];
  mu];

(* a claim is the computed bound weakened by one part in 10^9: grids are placed with library functions (Log, Exp, cube
   roots), which may differ in the last digit between platforms and move the certified value by a few units of rounding *)
$skyrmeDegreeMargin = 1/10^9;

loadSkyrmeDegreeGrid[deg_Integer, muq_, preset_String] := Module[{dir = skyrmeDegreeGridDir[deg, muq, preset], grid},
  If[! DirectoryQ[dir], Throw[failure["MissingData", "Grid data not found: " <> dir], "CertiField"]];
  grid = loadSkyrmeGrid[dir];
  If[! (grid["Meta"]["degree"] === deg && grid["Meta"]["mu_exact"] === $skyrmeDegreeMuNames[muq]),
    Throw[failure["BadGrid", "The grid in " <> dir <> " was made for another degree or pion mass."], "CertiField"]];
  grid];

certifySkyrmeDegreeTwo[params_Association, cand_, opts_List] := Module[
  {muq = skyrmeDegreeMu[params], preset = OptionValue[CertifiedVariationalBound, opts, "Resolution"],
   verbose = TrueQ[OptionValue[CertifiedVariationalBound, opts, "Verbose"]], grid, tab, lbRes, trial, ubRes,
   t0 = AbsoluteTime[], monitor, lo, hi},
  If[! KeyExistsQ[$skyrmeDegreePresets, preset],
    Throw[failure["BadOption", "\"Resolution\" must be one of " <> ToString[Keys[$skyrmeDegreePresets]] <> " for this model."], "CertiField"]];
  grid = loadSkyrmeDegreeGrid[2, muq, preset];
  tab = minorantTable[$skyrmeDegreePresets[preset, "Table"]];
  monitor = If[verbose, Function[s, If[Mod[s["Level"], 250] == 0, Print["  level ", s["Level"], "/", s["Levels"], "   pairs ", s["Pairs"], "   ", Round[s["Seconds"], 0.1], " s"]]], None];
  lbRes = levelSetLowerBoundDegree[2, muq, grid, tab, "Monitor" -> monitor];
  trial = Which[
    cand === Automatic, Get[skyrmeTrialFile[muq]],
    AssociationQ[cand] && KeyExistsQ[cand, "r"], cand,
    True, Module[{f = cand, df}, df = Derivative[1][f]; hedgehogTrialData[f, df, muq]]];
  ubRes = hedgehogUpperBound[trial, muq];
  lo = exactOf[lbRes["LowerBound"]];
  hi = 2 exactNumber[ubRes["UpperBound"]];
  <|"LowerBound" -> lo - Abs[lo] $skyrmeDegreeMargin, "UpperBound" -> hi + Abs[hi] $skyrmeDegreeMargin,
    "Seconds" -> AbsoluteTime[] - t0, "Computed" -> {lbRes["LowerBound"], hi},
    "LevelSet" -> lbRes, "Hedgehog" -> KeyDrop[ubRes, "Enclosure"], "Trial" -> trial, "Preset" -> preset, "Mu" -> muq,
    "Degree" -> 2, "GridMeta" -> grid["Meta"], "TableInfo" -> tab["Info"]|>];

skyrmeDegreeTwoDescription[res_] := <|
  "Units" -> "12\[Pi]^2 (Faddeev\[Dash]Bogomolny units), Battye\[Dash]Sutcliffe normalisation",
  "Theorem" -> Row[{"Every finite-energy map U: \[DoubleStruckCapitalR]^3 \[RightArrow] S^3 of degree 2 or \[Minus]2 with pion mass \[Mu] = ", exactString[res["Mu"]],
    " has E(U)/(12\[Pi]^2) \[GreaterEqual] ", lowerString[res["LowerBound"], 8], ".  There is a map of degree two (two translated copies of a hedgehog trial map with disjoint supports) with E/(12\[Pi]^2) \[LessEqual] ",
    upperString[res["UpperBound"], 8], "."}],
  "AdmissibleClass" -> "all finite-energy maps \[DoubleStruckCapitalR]^3 \[RightArrow] S^3 of degree \[PlusMinus]2 with U \[Minus] 1 \[Element] L^6 (no symmetry assumed)",
  "LowerBoundMethod" -> Row[{"level-set bound with the flux of degree two (coarea, isoperimetry, affine minorants of the level energy), certified dynamic programme: ",
    res["LevelSet", "Levels"], " levels, ", res["LevelSet", "Pairs"], " transitions, outward-rounded arithmetic; the reported number is the computed one lowered by one part in 10^9"}],
  "UpperBoundMethod" -> Row[{"twice the certified energy of the hedgehog trial map (", res["Hedgehog", "Cells"], " cubic cells, rigorous quadrature); the trial map is the vacuum outside a ball, so two disjoint copies are one map of degree two"}],
  "TrustedBase" -> "IEEE-754 machine arithmetic (+, \[Minus], \[Times], \[Divide]), Wolfram Language Interval arithmetic, and the analytic reduction (docs/SKYRME_DEGREE_TWO.md, Theorems 1 and 2)",
  "NotClaimed" -> "that either number is close to the infimum I_2: published numerical values of I_2 lie about 6% above this lower bound (massless case); the upper bound is the trivial one from two separated Skyrmions and is not a competitor for the minimiser; nothing is claimed about the existence or the shape of a minimiser",
  "References" -> "docs/SKYRME_DEGREE_TWO.md; the level-set inequality for degree one and its certified programme are those of SkyrmeDegreeOne."|>;

verifySkyrmeDegreeTwo[a_Association] := Module[{c = a["Certificate"], lb, ub},
  lb = levelSetLowerBoundDegree[2, c["Mu"], loadSkyrmeDegreeGrid[2, c["Mu"], c["Preset"]],
    minorantTable[$skyrmeDegreePresets[c["Preset"], "Table"]]]["LowerBound"];
  ub = 2 exactNumber[hedgehogUpperBound[c["Trial"], c["Mu"]]["UpperBound"]];
  NumericQ[lb] && NumericQ[ub] && exactOf[lb] >= exactNumber[a["LowerBound"]] && ub <= exactNumber[a["UpperBound"]]];

(* ---- stored claims (Data/SkyrmeDegree/claims.json): recompute and compare ---- *)
skyrmeDegreeClaims[] := Import[FileNameJoin[{$CertiFieldDirectory, "Data", "SkyrmeDegree", "claims.json"}], "RawJSON"];

(* an exact rational from a decimal string such as "2.2084" or "0"; nothing else is accepted *)
skyrmeDegreeRational[s_String] := Module[{m = StringCases[s, RegularExpression["^(-?)([0-9]+)(?:\\.([0-9]+))?$"] -> {"$1", "$2", "$3"}]},
  If[m === {}, Throw[failure["BadData", "Not a decimal number: " <> s], "CertiField"]];
  With[{sg = If[m[[1, 1]] === "-", -1, 1], ip = FromDigits[m[[1, 2]]], fr = m[[1, 3]]},
    sg (ip + If[fr === "", 0, FromDigits[fr]/10^StringLength[fr]])]];
skyrmeDegreeRational[_] := Throw[failure["BadData", "A decimal string was expected."], "CertiField"];

(* True if the bound recomputed from the stored grid is at least the stored claim *)
verifySkyrmeDegreeClaim[name_String] := Catch[Module[{cl = skyrmeDegreeClaims[][name], muq, lb},
  muq = skyrmeDegreeRational[cl["mu"]];
  lb = levelSetLowerBoundDegree[cl["degree"], muq, loadSkyrmeDegreeGrid[cl["degree"], muq, cl["preset"]],
    minorantTable[$skyrmeDegreePresets[cl["preset"], "Table"]]]["LowerBound"];
  NumericQ[lb] && exactOf[lb] >= skyrmeDegreeRational[cl["lower"]]], "CertiField"];

End[];
