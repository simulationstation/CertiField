(* ::Package:: *)

(* ::Title:: *)
(*CertiField: the level-set lower bound for the Skyrme energy in degree one*)

(* ::Text:: *)
(*Wolfram Language re-implementation of the computer-assisted part of Theorem 3 of*)
(*"A level-set lower bound for the Skyrme energy" (note_levelset.tex, programme levelset_cert3.py).*)
(**)
(*Theorem (analytic part, proved in the note).  For every finite-energy U: R^3 -> S^3 of degree one,*)
(*    E(U) >= W(tau_1, v)  for some v >= 0,*)
(*where W is the value function of the relaxed level-set problem.  The programme below builds, backwards in tau,*)
(*functions Lambda_i <= W(tau_i, .) that are affine on the cells of a volume grid, using*)
(*  (M1)-(M3) certified affine minorants of h_kappa,*)
(*  the second-order step bound (Lemma 5, condition [*]),*)
(*  the closed-form floors Xi, Xi_mu, Bog (Lemma 4),*)
(*  monotonicity of W in the volume,*)
(*and returns  min_v Lambda_1(v) / (12 pi^2), rounded down.*)
(**)
(*Rigour: tau-dependent scalars are exact rationals (Taylor series with remainder, 35-digit pi); array arithmetic is*)
(*outward-rounded (RigorousArithmetic.wl); the minorant table is certified entry by entry (Minorants.wl); condition*)
(*[*] is checked with enclosures and the step falls back to the first-order bound when it fails.  Grids, windows,*)
(*tangent columns and multipliers are heuristic and do not affect validity.*)

Begin["CertiField`Private`"];

$fourPi = 4. Pi;
$g0 = Power[36. Pi, 1/3.];
gFun[v_] := $g0 Power[v, 1/3.]^2;           (* heuristic only *)

(* constants with exact enclosures (computed once) *)
skyrmeConstants[] := skyrmeConstants[] = Module[{g0, cr},
  g0 = cubeRootBounds[36 $piLo, 36 $piHi];
  cr = cubeRootBounds[3/(4 $piHi), 3/(4 $piLo)];
  <|"G0LO" -> g0[[1]], "CRLO" -> cr[[1]], "TWOTHIRDSLO" -> fLo[2/3], "TWELVEPI2HI" -> fHi[12 $piHi^2]|>];

gLoV[v_, c_] := With[{y = cbrtLo[v]}, Ramp[rDn[c["G0LO"] Ramp[rDn[y y]]]]];
gPrimeLoV[w_, c_] := Ramp[rDn[rDn[c["G0LO"] c["TWOTHIRDSLO"]]/cbrtHi[w]]];
radLoV[v_, c_] := Ramp[rDn[c["CRLO"] cbrtLo[v]]];

slopeLo[ya_, yb_, dvIv_] := First[ivDivPos[ivSub[ivPt[yb], ivPt[ya]], dvIv]];

$skyrmeMuExact = <|0 -> 0, 263/500 -> 263/500, 13/20 -> 13/20, 1 -> 1|>;

(* ---- grid data (heuristic, shipped as exact IEEE doubles) ---- *)
readF64[file_] := Developer`ToPackedArray[BinaryReadList[file, "Real64", ByteOrdering -> -1]];
readI64[file_] := Developer`ToPackedArray[BinaryReadList[file, "Integer64", ByteOrdering -> -1]];

loadSkyrmeGrid[dir_String] := Module[{meta},
  meta = Import[FileNameJoin[{dir, "meta.json"}], "RawJSON"];
  <|"Meta" -> meta, "taus" -> readF64[FileNameJoin[{dir, "taus.f64"}]], "vref" -> readF64[FileNameJoin[{dir, "vref.f64"}]],
    "eug" -> readF64[FileNameJoin[{dir, "eug.f64"}]], "offs" -> readI64[FileNameJoin[{dir, "offs.i64"}]]|>];

(* ---- the verifier ---- *)
Options[levelSetLowerBound] = {"m1" -> 8., "m2" -> 3., "nmin" -> 6, "Monitor" -> None};

levelSetLowerBound[muq_, grid_Association, tab_Association, OptionsPattern[]] := Module[
  {cst = skyrmeConstants[], m1 = N@OptionValue["m1"], m2 = N@OptionValue["m2"], nmin = OptionValue["nmin"],
   monitor = OptionValue["Monitor"], mu2q = muq^2, muLo = fLo[muq], taus = grid["taus"], vref = grid["vref"],
   eug = grid["eug"], offs = grid["offs"], nLev, tq, S1, C1, S2, C2, CH, XiC, XmC, Vlo, Vhi, nodes, xiLo, ximLo,
   vb, x, xm, dvb, useM, YL, slp, Bog, T, PRE, nfar = 0, nfallback = 0, npairs = 0, t0 = AbsoluteTime[],
   ta, tb, dt, dq, halfDt, s2a, s2b, c2a, c2b, b0, bbq, baq, C0, Cb, Ca, VbLo, VbHi, VbIv, mterm, ca, cb, va, vL, vR,
   J, Jb, gL, gR, vmid, dA, wL1, wR1, end1, F1, suf, klo, shift, vtop, khi, cnt, starts, P, jj, kk, wlo, whi, valid,
   Yk, Sk, w0, wmid, PREk, lift, dB, dAj, cm, scale, dUp, KUp, KDn, Kt, alB, beB, khB, gam, al1, dk, alA, beA, khA,
   rhs, good, KIv, hb, alsum, BL, BR, alF, beF, khF, beFar, farStep, far1, suf2, base2, far6, f2L, f2R, use2, farL,
   farR, xa, xma, dvIv, dv, candY, candS, mid, best, idx, msk, YLn, Sn, endn, cmax, pre, ktop, fin, val, cert,
   endBound, vend, gend},

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
    (* exact lower bounds of  int c,  int (1-t) c,  int t c  over the step, c = 4 pi sin^2 *)
    b0 = dq - (s2b[[2]] - s2a[[1]])/2;
    bbq = dq^2/2 - dq s2b[[2]]/2 - (c2b[[2]] - c2a[[1]])/4;
    baq = dq^2/2 + dq s2a[[1]]/2 + (c2b[[1]] - c2a[[2]])/4;
    If[! (b0 > 0), Throw[failure["Internal", "b0 <= 0"], "CertiField"]];
    C0 = fLo[2 $piLo b0]; Cb = fLo[2 $piLo Max[bbq, 0]/dq]; Ca = fLo[2 $piLo Max[baq, 0]/dq];
    VbLo = Vlo[[i + 2]]; VbHi = Vhi[[i + 2]]; VbIv = {VbLo, VbHi};
    mterm = fLo[2 mu2q Min[S1[[i + 1, 1]], S1[[i + 2, 1]]]];
    ca = $fourPi Sin[ta]^2; cb = $fourPi Sin[tb]^2;                          (* heuristics only *)
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
    If[monitor =!= None, monitor[<|"Level" -> i, "Levels" -> nLev, "Tau" -> ta, "Pairs" -> npairs,
      "Lambda0" -> YL[[1]]/(12 Pi^2), "Seconds" -> AbsoluteTime[] - t0|>]],
    {i, nLev - 1, 0, -1}];

  fin = cumMax[vMax[vMax[vMax[vMin[BL, BR], Most[xa]], Most[xma]], Bog]];
  val = Min[Min[fin], T];
  cert = rDn[val/cst["TWELVEPI2HI"]];
  <|"LowerBound" -> cert, "Levels" -> nLev, "VolumeNodes" -> Length[eug], "Pairs" -> npairs,
    "Fallbacks" -> nfallback, "FarBinding" -> nfar, "Seconds" -> AbsoluteTime[] - t0|>];

End[];
