(* ::Package:: *)

(* ::Title:: *)
(*CertiField: the Skyrme model in degree one for a general pion mass*)

(* ::Text:: *)
(*Statement and proofs: docs/SKYRME_ANY_MASS.md.*)
(**)
(*Everything numerical in this file is HEURISTIC.  It chooses the grids of the level-set programme and the data of the*)
(*hedgehog trial map for a pion mass for which no data are shipped.  The certified numbers are produced, unchanged,*)
(*by levelSetLowerBound (SkyrmeLevelSet.wl) and hedgehogUpperBound (SkyrmeHedgehog.wl), which accept any grid and*)
(*any rational trial data: a poor grid or profile gives a poor bound, never a wrong one.*)
(**)
(*The grids are those of the original verifier named in SkyrmeLevelSet.wl; its construction is reproduced here:*)
(*  1. the extremal of the relaxed level-set problem, by shooting on the central slope a.  With*)
(*     h_kappa(D) = 2 Sqrt[(1 + 2D)(D^2 + 2D + kappa)] and H_kappa its convex envelope,*)
(*         r' = dH/dkappa,   kappa' = -2 mu^2 Sin[tau] - (2/r) (H - D dH/dD),   D = Sin[tau]^2/r^2,*)
(*         r ~ tau/a,  kappa = a^4 - a^2  at the centre;*)
(*  2. the reference radius: monotone cubic interpolation of Log[r] in Log[tau/(Pi - tau)];*)
(*  3. levels tau_1 = 0.03 ... Pi - 0.02 with steps Min[0.004, delta/lambda], lambda = d Log[v_ref]/d tau;*)
(*  4. volume nodes v_ref(tau_i) Exp[u] on a fixed set of u, cut off below the radius 0.1 Sin[tau_i].*)
(*With h in place of H the system of step 1 is the hedgehog equation (tau = Pi - f), which gives the trial profile.*)

Begin["CertiField`Private`"];

(* ------------------------------------------------------------------------------------------------ *)
(* the extremal system                                                                               *)
(* ------------------------------------------------------------------------------------------------ *)

(* Right-hand side for the state {r, kappa, cost}.  h_kappa is convex iff (3 - Sqrt[5])/2 <= kappa <= (3 + Sqrt[5])/2.
   Otherwise the envelope is the tangent line through (0, 2 Sqrt[kappa]) below the tangent point D_t(kappa), the
   root of the increasing function  D^2 - 1 - kappa + Sqrt[kappa (4 D + 5)]  (the tangency condition
   h - D h' = 2 Sqrt[kappa] with the square root removed).  Where h - D h' < 0, that is D > D* = (1 + Sqrt[1 + 4 kappa])/2,
   the area constraint is inactive. *)
skRhs = With[{k1 = (3. - Sqrt[5.])/2., k2 = (3. + Sqrt[5.])/2., fourPi = 4. Pi},
  Compile[{{tau, _Real}, {y, _Real, 1}, {mu2, _Real}, {relaxed, True | False}},
    Module[{r = Max[y[[1]], 1.*^-12], s = Sin[tau], pot = 0., kap = 0., d = 0., dt = 0., two = False, lo = 0., hi = 0.,
        mid = 0., lam = 0., yk = 0., q = 0., hs = 0., hk = 0., ic = 0., ds = 0., rho = 0., rp = 0., kp = 0., cp = 0.},
      pot = 2. mu2 (1. + Cos[tau]);
      kap = Max[y[[2]], pot, 1.*^-300];
      d = s s/(r r);
      If[relaxed && ! (k1 <= kap && kap <= k2),
        If[kap <= 1.*^-16,
          dt = 1.; two = d < dt,
          If[d d - 1. - kap + Sqrt[kap (4. d + 5.)] < 0.,
            two = True; lo = d; hi = Sqrt[1. + kap];
            Do[mid = 0.5 (lo + hi); If[mid mid - 1. - kap + Sqrt[kap (4. mid + 5.)] < 0., lo = mid, hi = mid], {60}];
            dt = 0.5 (lo + hi)]]];
      If[two,
        lam = d/dt; yk = Sqrt[kap];
        q = (1. + 2. dt) (dt dt + 2. dt + kap);
        hs = (1. - lam) 2. yk + lam 2. Sqrt[q];
        hk = (1. - lam)/yk + lam Sqrt[(1. + 2. dt)/(dt dt + 2. dt + kap)];
        ic = 2. yk,
        q = (1. + 2. d) (d d + 2. d + kap);
        hs = 2. Sqrt[q];
        hk = Sqrt[(1. + 2. d)/(d d + 2. d + kap)];
        ic = hs - d (2. (d d + 2. d + kap) + (1. + 2. d) (2. d + 2.))/Sqrt[q]];
      If[ic < 0.,
        ds = 0.5 (1. + Sqrt[1. + 4. kap]);
        rho = 2. Sqrt[(1. + 2. ds) (ds ds + 2. ds + kap)]/ds;
        rp = s s Sqrt[(1. + 2. ds)/(ds ds + 2. ds + kap)]/ds/(r r);
        kp = -2. mu2 s;
        cp = fourPi s s rho - (kap - pot) fourPi r r rp,
        rp = hk;
        kp = -2. mu2 s - 2./r ic;
        cp = fourPi r r (hs - (kap - pot) rp)];
      {rp, kp, cp}]]];

skFiniteQ[v_] := VectorQ[v, Developer`MachineRealQ] && Max[Abs[v]] < 1.*^150;

(* one chunk [t0, t1] with the embedded Runge-Kutta pair of Dormand and Prince (orders 5 and 4), relative tolerance
   10^-10 and absolute tolerance 10^-13 as in the original (which uses LSODA).  Returns {y(t1), last step} or $Failed. *)
skIntegrate[mu2_, relaxed_, t0_, y0_, t1_, h0_] := Module[
  {t = t0, y = y0, h = Min[h0, t1 - t0], k1, k2, k3, k4, k5, k6, k7, yn, ev, en, steps = 0, hLast = h0, last = False, done = False},
  k1 = skRhs[t, y, mu2, relaxed];
  While[! done,
    If[++steps > 20000 || ! TrueQ[h > 1.*^-15 (1. + Abs[t])] || ! skFiniteQ[k1], Return[$Failed, Module]];
    last = t + h >= t1 - 1.*^-14 (t1 - t0);
    If[last, h = t1 - t];
    k2 = skRhs[t + 0.2 h, y + h (0.2 k1), mu2, relaxed];
    k3 = skRhs[t + 0.3 h, y + h (0.075 k1 + 0.225 k2), mu2, relaxed];
    k4 = skRhs[t + 0.8 h, y + h ((44./45.) k1 - (56./15.) k2 + (32./9.) k3), mu2, relaxed];
    k5 = skRhs[t + (8./9.) h, y + h ((19372./6561.) k1 - (25360./2187.) k2 + (64448./6561.) k3 - (212./729.) k4), mu2, relaxed];
    k6 = skRhs[t + h, y + h ((9017./3168.) k1 - (355./33.) k2 + (46732./5247.) k3 + (49./176.) k4 - (5103./18656.) k5), mu2, relaxed];
    yn = y + h ((35./384.) k1 + (500./1113.) k3 + (125./192.) k4 - (2187./6784.) k5 + (11./84.) k6);
    If[! skFiniteQ[yn], h = 0.2 h; Continue[]];
    k7 = skRhs[t + h, yn, mu2, relaxed];
    ev = h ((71./57600.) k1 - (71./16695.) k3 + (71./1920.) k4 - (17253./339200.) k5 + (22./525.) k6 - (1./40.) k7);
    If[! skFiniteQ[ev], h = 0.2 h; Continue[]];
    en = Sqrt[Total[(ev/(1.*^-13 + 1.*^-10 Max /@ Transpose[{Abs[y], Abs[yn]}]))^2]/3.];
    If[en <= 1.,
      t = If[last, t1, t + h]; y = yn; k1 = k7; hLast = h; If[last, done = True]];
    h = h Min[5., Max[0.2, 0.9 Max[en, 1.*^-10]^(-0.2)]]];
  {y, hLast}];

(* The path with central slope a, sampled at the ends of chunks of length Min[0.02, 0.2 (Pi - tau)], from tau = 10^-4
   to Pi - 2 10^-4.  "Early": the multiplier kappa - V reached 0, or the radius exceeded 10^7, before the end. *)
skShoot[a_?NumericQ, mu_?NumericQ, relaxed : (True | False)] := Quiet@Module[
  {mu2 = N[mu]^2, aa = N[a], tEnd = N[Pi] - 2.*^-4, t = 1.*^-4, y, h = 1.*^-5, early = False, dt, res, bag},
  y = {t/aa, aa^4 - aa^2, 0.};
  bag = Internal`Bag[{Prepend[y, t]}];
  While[t < tEnd - 1.*^-15,
    dt = Min[0.02, 0.2 (N[Pi] - t), tEnd - t];
    res = skIntegrate[mu2, relaxed, t, y, t + dt, h];
    If[res === $Failed, early = True; Break[]];
    {y, h} = res; t = t + dt;
    Internal`StuffBag[bag, Prepend[y, t]];
    If[! skFiniteQ[y] || y[[2]] - 2. mu2 (1. + Cos[t]) <= 0. || y[[1]] > 1.*^7, early = True; Break[]]];
  res = Developer`ToPackedArray[N[Internal`BagPart[bag, All]]];
  <|"a" -> aa, "t" -> res[[All, 1]], "r" -> res[[All, 2]], "kappa" -> res[[All, 3]], "cost" -> res[[All, 4]], "Early" -> early|>];

(* sign used for the bisection on the slope: -1 if the path ends early; for mu = 0 the sign of r Sqrt[kappa] - (Pi - tau)
   at the end (constant along the far tail of the massless relaxed problem, zero on the extremal); +1 otherwise *)
skDiscriminant[res_Association, mu_] := Which[
  res["Early"], -1.,
  mu == 0, Last[res["r"]] Sqrt[Max[Last[res["kappa"]], 0.]] - (N[Pi] - Last[res["t"]]),
  True, 1.];

skCriticalSlope[mu_, relaxed_] := Module[{grid, vals, i, lo, hi, vlo, mid, vm, f},
  f = skDiscriminant[skShoot[#, mu, relaxed], mu] &;
  grid = Subdivide[1., 7., 48];
  vals = {f[grid[[1]]]};
  i = 0;
  Do[AppendTo[vals, f[grid[[k]]]]; If[vals[[k - 1]] vals[[k]] < 0, i = k - 1; Break[]], {k, 2, Length[grid]}];
  If[i == 0, Return[$Failed]];
  {lo, hi, vlo} = {grid[[i]], grid[[i + 1]], vals[[i]]};
  Do[mid = 0.5 (lo + hi);
    If[mid == lo || mid == hi, Break[]];
    vm = f[mid];
    If[vm vlo > 0, lo = mid, hi = mid], {60}];
  0.5 (lo + hi)];

(* the extremal: the first path next to the critical slope that reaches the end *)
skExtremal[mu_, relaxed_] := skExtremal[mu, relaxed] = Module[{a = skCriticalSlope[mu, relaxed], res = $Failed, r},
  If[a =!= $Failed,
    Do[r = skShoot[a + da, mu, relaxed];
      If[! r["Early"] && OrderedQ[r["t"]] && OrderedQ[r["r"]] && Min[r["r"]] > 0., res = Append[r, "CriticalSlope" -> a]; Break[]],
      {da, {2.*^-9, -2.*^-9, 2.*^-8, -2.*^-8, 2.*^-7, -2.*^-7}}]];
  res];

(* ------------------------------------------------------------------------------------------------ *)
(* monotone cubic interpolation (the derivative rule of Fritsch and Butland, as in the original)        *)
(* ------------------------------------------------------------------------------------------------ *)

skPchipSlopes[x_, y_] := Module[{h = Differences[x], m, n = Length[x], d, w1, w2, edge},
  m = Differences[y]/h;
  d = ConstantArray[0., n];
  w1 = 2. Rest[h] + Most[h]; w2 = Rest[h] + 2. Most[h];
  Do[d[[k + 1]] = If[Sign[m[[k]]] != Sign[m[[k + 1]]] || m[[k]] == 0. || m[[k + 1]] == 0., 0.,
      (w1[[k]] + w2[[k]])/(w1[[k]]/m[[k]] + w2[[k]]/m[[k + 1]])], {k, n - 2}];
  edge[h0_, h1_, m0_, m1_] := Module[{dd = ((2. h0 + h1) m0 - h0 m1)/(h0 + h1)},
    Which[Sign[dd] != Sign[m0], 0., Sign[m0] != Sign[m1] && Abs[dd] > 3. Abs[m0], 3. m0, True, dd]];
  d[[1]] = edge[h[[1]], h[[2]], m[[1]], m[[2]]];
  d[[n]] = edge[h[[-1]], h[[-2]], m[[-1]], m[[-2]]];
  Developer`ToPackedArray[N[d]]];

(* {values, derivatives} of the interpolant at the points xq (a list inside [x_1, x_n]) *)
skPchipEval[x_, y_, d_, xq_List] := Module[{n = Length[x], idx, h, t, y0, y1, d0, d1},
  idx = Clip[searchRight[x, Developer`ToPackedArray[N[xq]]], {1, n - 1}];
  h = x[[idx + 1]] - x[[idx]]; t = (xq - x[[idx]])/h;
  y0 = y[[idx]]; y1 = y[[idx + 1]]; d0 = d[[idx]]; d1 = d[[idx + 1]];
  {(2. t^3 - 3. t^2 + 1.) y0 + (t^3 - 2. t^2 + t) h d0 + (-2. t^3 + 3. t^2) y1 + (t^3 - t^2) h d1,
   (6. t^2 - 6. t) (y0 - y1)/h + (3. t^2 - 4. t + 1.) d0 + (3. t^2 - 2. t) d1}];

(* ------------------------------------------------------------------------------------------------ *)
(* the grids                                                                                         *)
(* ------------------------------------------------------------------------------------------------ *)

(* Log[r_ref] as a function of s = Log[tau/(Pi - tau)], continued linearly outside the computed range *)
skReferencePath[mu_] := skReferencePath[mu] = Module[{ext = skExtremal[mu, True], s, y, d},
  If[ext === $Failed, Return[$Failed, Module]];
  s = Log[ext["t"]/(N[Pi] - ext["t"])]; y = Log[ext["r"]];
  d = skPchipSlopes[s, y];
  <|"s" -> s, "y" -> y, "d" -> d, "s0" -> First[s], "s1" -> Last[s], "Slope" -> ext["a"], "CriticalSlope" -> ext["CriticalSlope"],
    "Value" -> Last[ext["cost"]]/(12. Pi^2)|>];

skLnR[p_Association, tau_List] := Module[{s = Log[tau/(N[Pi] - tau)]},
  First[skPchipEval[p["s"], p["y"], p["d"], Clip[s, {p["s0"], p["s1"]}]]] + Last[p["d"]] Ramp[s - p["s1"]] - First[p["d"]] Ramp[p["s0"] - s]];

(* d Log[v_ref]/d tau *)
skLam[p_Association, tau_?NumericQ] := With[{s = Log[tau/(N[Pi] - tau)]},
  3. skPchipEval[p["s"], p["y"], p["d"], {Clip[s, {p["s0"], p["s1"]}]}][[2, 1]] (1./tau + 1./(N[Pi] - tau))];

skUGrid[du_, u1_, u2_, growth_, dumax_, ulo_, uhi_] := Module[{up = {0.}, dn = {0.}, e = du},
  While[Last[up] < uhi, AppendTo[up, Last[up] + e]; If[Last[up] >= u2, e = Min[e growth, dumax]]];
  e = du;
  While[Last[dn] > ulo, AppendTo[dn, Last[dn] - e]; If[Last[dn] <= -u1, e = Min[e growth, dumax]]];
  Developer`ToPackedArray[N[Join[Reverse[Rest[dn]], up]]]];

(* number of entries of the sorted list a that are < x *)
searchLeft[a_, x_] := Length[a] - searchRight[-Reverse[a], -x];

$skyrmeGridDefaults = <|"U1" -> 4., "U2" -> 2., "growth" -> 1.05, "dumax" -> 0.2, "ulo" -> -32., "uhi" -> 22., "t1" -> 0.03,
  "lN" -> 0.02, "eps1" -> 0.1, "dtmax" -> 0.004|>;
$skyrmeGridDelta = <|"Quick" -> 0.04, "Standard" -> 0.02, "Publication" -> 0.005|>;

(* the grid for the mass mu at the spacing delta, in the form returned by loadSkyrmeGrid; $Failed if the shooting fails *)
skyrmeGenerateGrid[mu_, delta_?NumericQ, preset_String : "Custom"] := Module[
  {p = skReferencePath[mu], g = $skyrmeGridDefaults, tN, taus, cur, ug, i0, lnr, vref, eug, ulow, offs},
  If[p === $Failed, Return[$Failed, Module]];
  tN = N[Pi] - g["lN"];
  taus = Internal`Bag[{g["t1"]}];
  cur = g["t1"];
  While[cur < tN - 1.*^-12,
    cur = Min[cur + Min[g["dtmax"], delta/Max[skLam[p, cur], 1.*^-6]], tN];
    Internal`StuffBag[taus, cur]];
  taus = Developer`ToPackedArray[N[Internal`BagPart[taus, All]]];
  ug = skUGrid[N[delta], g["U1"], g["U2"], g["growth"], g["dumax"], g["ulo"], g["uhi"]];
  i0 = First[Ordering[Abs[ug], 1]] - 1;
  lnr = skLnR[p, taus];
  vref = (4. Pi/3.) Exp[3. lnr]; eug = Exp[ug];
  ulow = 3. Log[g["eps1"] Sin[taus]/Exp[lnr]];
  offs = Clip[searchLeft[ug, ulow], {0, i0 - 1}];
  offs = Reverse[Rest[FoldList[Max, 0, Reverse[offs]]]];
  If[! (VectorQ[taus, Developer`MachineRealQ] && VectorQ[vref, Developer`MachineRealQ] && Min[vref] > 0. && Max[vref] < 1.*^200 &&
      OrderedQ[taus] && Min[Differences[taus]] > 0.), Return[$Failed, Module]];
  <|"Meta" -> <|"model" -> "SkyrmeDegreeOne", "mu_exact" -> ToString[mu, InputForm], "preset" -> ToLowerCase[preset],
      "params" -> Join[<|"du" -> N[delta], "cg" -> N[delta]|>, g], "N" -> Length[taus] - 1, "n_u" -> Length[ug], "i0" -> i0,
      "slope" -> p["Slope"], "note" -> "Generated by Kernel/SkyrmeAnyMass.wl. Grid data are heuristic; they affect the tightness of the bound, not its validity."|>,
    "taus" -> taus, "vref" -> Developer`ToPackedArray[vref], "eug" -> Developer`ToPackedArray[eug],
    "offs" -> Developer`ToPackedArray[offs]|>];

(* What levelSetLowerBound assumes about a grid without testing it.  The shipped grids satisfy this; a grid that comes
   from a certificate is tested here before it is used.  (The programme itself checks that the levels are increasing
   in (0, Pi) with exactly representable steps and that the volume cells of every level but the last are not
   degenerate.) *)
skyrmeGridQ[g_] := AssociationQ[g] && Module[{taus, vref, eug, offs, n, vb},
  {taus, vref, eug, offs} = Lookup[g, {"taus", "vref", "eug", "offs"}, None];
  VectorQ[taus, Developer`MachineRealQ] && VectorQ[vref, Developer`MachineRealQ] && VectorQ[eug, Developer`MachineRealQ] &&
    VectorQ[offs, IntegerQ] && (n = Length[taus]) >= 2 && Length[vref] == n && Length[offs] == n && Length[eug] >= 2 &&
    Min[vref] > 0. && Max[vref] < 1.*^200 && Min[eug] > 0. && Max[eug] < 1.*^200 && Min[Differences[eug]] > 0. &&
    Min[offs] >= 0 && Max[offs] <= Length[eug] - 1 &&
    (vb = Join[{0.}, Last[vref] eug[[Last[offs] + 1 ;;]]];
      Min[First[ivSub[ivPt[Rest[vb]], ivPt[Most[vb]]]]] > 0.)];

skyrmePackGrid[g_Association] := Join[g, <|"taus" -> Developer`ToPackedArray[g["taus"]], "vref" -> Developer`ToPackedArray[g["vref"]],
  "eug" -> Developer`ToPackedArray[g["eug"]], "offs" -> Developer`ToPackedArray[g["offs"]]|>];

(* a grid in the format of Data/Skyrme (raw little-endian doubles and 64-bit integers), e.g. for another verifier *)
exportSkyrmeGrid[dir_String, g_Association] := Module[{w},
  If[! DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];
  w[name_, data_, type_] := Module[{s = OpenWrite[FileNameJoin[{dir, name}], BinaryFormat -> True]},
    BinaryWrite[s, data, type, ByteOrdering -> -1]; Close[s]];
  w["taus.f64", g["taus"], "Real64"]; w["vref.f64", g["vref"], "Real64"]; w["eug.f64", g["eug"], "Real64"];
  w["offs.i64", g["offs"], "Integer64"];
  Export[FileNameJoin[{dir, "meta.json"}], g["Meta"], "RawJSON"];
  dir];

(* ------------------------------------------------------------------------------------------------ *)
(* the trial map                                                                                     *)
(* ------------------------------------------------------------------------------------------------ *)

(* The hedgehog profile f(r) = Pi - tau(r) and its derivative: the path of skExtremal[mu, False] inverted by monotone
   interpolation of Log[tau/(Pi - tau)] in Log[r], down to f = $skTailStart; beyond that radius the solution of the
   linearised equation, f proportional to (1 + mu r) Exp[-mu r]/r^2 (or 1/r^2 for mu = 0). *)
$skTailStart = 0.005;

skHedgehogProfile[mu_] := skHedgehogProfile[mu] = Module[
  {ext = skExtremal[mu, False], t, r, x, s, d, kc, rc, fc, r1, a, m = N[mu], tailT, tailD, tauOf, fFun, dfFun},
  If[ext === $Failed, Return[$Failed, Module]];
  t = ext["t"]; r = ext["r"]; a = ext["a"];
  x = Log[r]; s = Log[t/(N[Pi] - t)];
  d = skPchipSlopes[x, s];
  kc = Max[2, Total[UnitStep[N[Pi] - $skTailStart - t]]];
  rc = r[[kc]]; fc = N[Pi] - t[[kc]]; r1 = First[r];
  tailT[rr_] := If[m == 0., 1./rr^2, (1. + m rr) Exp[-m rr]/rr^2];
  tailD[rr_] := If[m == 0., -2./rr^3, -Exp[-m rr] (2. + 2. m rr + m^2 rr^2)/rr^3];
  tauOf[rr_] := Module[{ev = skPchipEval[x, s, d, {Log[rr]}], tau},
    tau = N[Pi]/(1. + Exp[-ev[[1, 1]]]);
    {tau, tau (1. - tau/N[Pi]) ev[[2, 1]]/rr}];
  fFun[rr_] := Which[rr <= r1, N[Pi] - a rr, rr <= rc, N[Pi] - First[tauOf[rr]], True, fc tailT[rr]/tailT[rc]];
  dfFun[rr_] := Which[rr <= r1, -a, rr <= rc, -Last[tauOf[rr]], True, fc tailD[rr]/tailT[rc]];
  <|"Mu" -> mu, "Slope" -> a, "TailRadius" -> rc, "f" -> Function[rr, fFun[N[rr]]], "df" -> Function[rr, dfFun[N[rr]]]|>];

(* rational trial data for hedgehogUpperBound; $Failed if the shooting fails *)
skyrmeGenerateTrial[mu_] := skyrmeGenerateTrial[mu] = Module[{p = skHedgehogProfile[mu], trial},
  If[p === $Failed, Return[$Failed, Module]];
  trial = Quiet@hedgehogTrialData[p["f"], p["df"], mu, "Radius" -> If[mu == 0, 60, Min[60, Ceiling[15/mu]]]];
  If[AssociationQ[trial] && AllTrue[Join[trial["r"], trial["p"], trial["dp"]], MatchQ[#, _Integer | _Rational] &], trial, $Failed]];

(* ------------------------------------------------------------------------------------------------ *)
(* which grid, which mass and which trial map for a requested mass (docs/SKYRME_ANY_MASS.md, Section 4) *)
(* ------------------------------------------------------------------------------------------------ *)

$skyrmeMuSup = 2;       (* masses mu with 0 <= mu < 2 are accepted *)
skyrmeMassQ[mu_] := MatchQ[mu, _Integer | _Rational] && 0 <= mu < $skyrmeMuSup;
skyrmeShippedMassQ[mu_] := KeyExistsQ[$skyrmeMuNames, mu];
skyrmeMonotoneMass[mu_] := Max[Select[Keys[$skyrmeMuNames], # <= mu &]];

skyrmeShippedGrid[mu_, preset_String] := With[{dir = skyrmeGridDir[mu, preset]},
  If[! DirectoryQ[dir], Throw[failure["MissingData", "Grid data not found: " <> dir], "CertiField"]];
  loadSkyrmeGrid[dir]];

skyrmeGeneratedGrid[mu_, preset_String] := skyrmeGeneratedGrid[mu, preset] =
  With[{g = Quiet@skyrmeGenerateGrid[mu, $skyrmeGridDelta[preset], preset]}, If[skyrmeGridQ[g], g, $Failed]];

(* "Mass": the mass at which the level-set programme is run (at most the requested one);
   "Extra": what is added to the certificate (nothing for a shipped grid at the requested mass) *)
skyrmeGridPlan[mu_, preset_String, method_, verbose_] := Module[{g, mu0},
  If[! MemberQ[{Automatic, "Monotonicity", "GeneratedGrid"}, method],
    Throw[failure["BadOption", "For SkyrmeDegreeOne, Method must be Automatic, \"Monotonicity\" or \"GeneratedGrid\"."], "CertiField"]];
  If[method === Automatic && skyrmeShippedMassQ[mu],
    Return[<|"Mass" -> mu, "Grid" -> skyrmeShippedGrid[mu, preset], "Extra" -> <||>|>, Module]];
  If[method =!= "Monotonicity",
    If[verbose, Print["  generating the grid for the pion mass ", mu]];
    g = skyrmeGeneratedGrid[mu, preset];
    If[AssociationQ[g],
      Return[<|"Mass" -> mu, "Grid" -> g, "Extra" -> <|"LowerBoundMass" -> mu, "GridSource" -> "Generated",
          "Grid" -> KeyTake[g, {"taus", "vref", "eug", "offs"}]|>|>, Module]];
    If[method === "GeneratedGrid",
      Throw[failure["Inconclusive", "The grid generator failed for this pion mass."], "CertiField"]];
    If[verbose, Print["  the grid generator failed; using the shipped grid of a smaller mass (monotonicity in the mass)"]]];
  mu0 = skyrmeMonotoneMass[mu];
  <|"Mass" -> mu0, "Grid" -> skyrmeShippedGrid[mu0, preset], "Extra" -> <|"LowerBoundMass" -> mu0, "GridSource" -> "Monotonicity"|>|>];

(* Any trial data give a valid upper bound at any mass.  If the profile solver fails, the shipped trial map of the
   nearest mass is used. *)
skyrmeTrial[mu_, verbose_] := Module[{trial},
  If[skyrmeShippedMassQ[mu], Return[Get[skyrmeTrialFile[mu]], Module]];
  If[verbose, Print["  computing a hedgehog profile for the pion mass ", mu]];
  trial = skyrmeGenerateTrial[mu];
  If[AssociationQ[trial], trial,
    If[verbose, Print["  the profile solver failed; using the shipped trial map of the nearest mass"]];
    Get[skyrmeTrialFile[First[SortBy[Keys[$skyrmeMuNames], Abs[# - mu] &]]]]]];

(* for VerifyCertificate: {mass of the lower bound, grid} from a certificate; the mass must not exceed the mass of
   the model, and a stored grid must pass skyrmeGridQ *)
skyrmeCertificateData[c_Association] := Module[{mu = c["Mu"], muL, g = Lookup[c, "Grid", None]},
  muL = Lookup[c, "LowerBoundMass", mu];
  If[! (MatchQ[mu, _Integer | _Rational] && MatchQ[muL, _Integer | _Rational] && 0 <= muL <= mu),
    Throw[failure["BadCertificate", "The lower bound was computed at a pion mass above that of the model."], "CertiField"]];
  If[g === None,
    If[! skyrmeShippedMassQ[muL], Throw[failure["BadCertificate", "No grid is stored and none is shipped for this mass."], "CertiField"]];
    {muL, skyrmeShippedGrid[muL, c["Preset"]]},
    If[! skyrmeGridQ[g], Throw[failure["BadCertificate", "The stored grid is not a valid grid."], "CertiField"]];
    {muL, skyrmePackGrid[g]}]];

End[];
