(* ::Package:: *)

(* ::Title:: *)
(*CertiField: certified upper bound for the Skyrme energy in degree one*)

(* ::Text:: *)
(*Units of note_levelset.tex: E(U) = Integrate[|dU|^2 + |Lambda^2 dU|^2 + 2 mu^2 (1 - cos f), R^3], reported as*)
(*E/(12 Pi^2).  For a hedgehog U = exp(i f(r) x.sigma/r),*)
(*    E = 4 Pi Integrate[ r^2 f'^2 + 2 Sin[f]^2 (1 + f'^2) + Sin[f]^4/r^2 + 2 mu^2 r^2 (1 - Cos[f]), {r, 0, oo} ].*)
(*Any admissible trial map gives an upper bound on the infimum I_1.  The trial used here is f = Pi - p with p a*)
(*piecewise cubic (Hermite) with rational data, p(0) = 0 exactly (so f(0) = Pi).  On the last cell p rises*)
(*linearly to Pi, and beyond it f == 0 (vacuum).  So the trial is continuous, has degree one and finite energy, and*)
(*E(trial) >= I_1 whatever the data.  The profile data are a numerical hedgehog; they affect only tightness.*)
(**)
(*Rigour:*)
(*  first cell [0, r1]: exact integral of the polynomial majorant (|sin p| <= |p|, 1 + cos p <= 2);*)
(*  middle cells: rigorous Simpson/first-order enclosures (RigorousQuadrature.wl);*)
(*  last cell: exact integral of the polynomial majorant (|sin f| <= |f|, 1 - cos f <= f^2/2, 1/r^2 <= 1/r_{n-1}^2),*)
(*  with |f| <= max |Pi - p_{n-1}| from a 35-digit rational enclosure of Pi.*)

Begin["CertiField`Private`"];

hedgehogDensityP[p_, dp_, r_, mu2_] := r^2 dp^2 + 2 Sin[p]^2 (1 + dp^2) + Sin[p]^4/r^2 + 2 mu2 r^2 (1 + Cos[p]);

hermiteCubic[{x0_, x1_}, {y0_, y1_}, {m0_, m1_}, x_] := Module[{h = x1 - x0, t},
  t = (x - x0)/h;
  Expand[(2 t^3 - 3 t^2 + 1) y0 + (t^3 - 2 t^2 + t) h m0 + (-2 t^3 + 3 t^2) y1 + (t^3 - t^2) h m1]];

(* ---- the hedgehog profile by shooting (heuristic) ---- *)
Options[hedgehogProfile] = {WorkingPrecision -> MachinePrecision, "RMax" -> Automatic, "Iterations" -> Automatic,
  "Slope" -> Automatic};

hedgehogProfile[mu_?NumericQ, OptionsPattern[]] := Module[
  {wp = OptionValue[WorkingPrecision], rmax, r0, iters, shoot, f, r, aLo = 1, aHi = 4, a, s, classify, sLo, sHi,
   good, rcut, tail, fcut, dfcut, prof, dprof},
  rmax = Replace[OptionValue["RMax"], Automatic :> If[mu == 0, 40, Min[36, 16/mu]]];
  iters = Replace[OptionValue["Iterations"], Automatic :> If[wp === MachinePrecision, 56, Ceiling[3.4 wp]]];
  r0 = SetPrecision[1/10^4, wp];
  shoot[aa_] := Quiet@NDSolveValue[{
      (r^2 + 2 Sin[f[r]]^2) f''[r] + 2 r f'[r] + Sin[2 f[r]] (f'[r]^2 - 1 - Sin[f[r]]^2/r^2) - mu^2 r^2 Sin[f[r]] == 0,
      f[r0] == Pi - aa r0, f'[r0] == -aa,
      WhenEvent[f[r] < 0, "StopIntegration"], WhenEvent[f'[r] > 0, "StopIntegration"]},
    f, {r, r0, rmax}, WorkingPrecision -> wp, MaxSteps -> Infinity, Method -> "StiffnessSwitching"];
  classify[sol_] := Module[{re = sol["Domain"][[1, 2]]}, Which[sol[re] < 0, 1, sol'[re] > 0, -1, True, 0]];
  Do[a = (aLo + aHi)/2; s = shoot[SetPrecision[a, wp]];
    Switch[classify[s], 1, aHi = a, -1, aLo = a, 0, aLo = a; aHi = a], {iters}];
  sLo = shoot[SetPrecision[aLo, wp]]; sHi = shoot[SetPrecision[aHi, wp]];
  (* reliable radius: where the bracketing solutions still agree to 1e-6 relative *)
  good = Select[Subdivide[r0, Min[sLo["Domain"][[1, 2]], sHi["Domain"][[1, 2]]], 4000],
    Abs[sLo[#] - sHi[#]] <= 10^-6 Abs[sLo[#]] &];
  rcut = 0.8 Max[good];
  fcut = sLo[rcut]; dfcut = sLo'[rcut];
  (* linearised tail beyond rcut: f ~ (1 + mu r) Exp[-mu r]/r^2 (mu > 0), ~ 1/r^2 (mu = 0) *)
  tail[rr_] := If[mu == 0, fcut (rcut/rr)^2, fcut (1 + mu rr) Exp[-mu rr] rcut^2/((1 + mu rcut) Exp[-mu rcut] rr^2)];
  prof[rr_?NumericQ] := Which[rr <= r0, Pi - aLo rr, rr <= rcut, sLo[rr], True, tail[rr]];
  dprof[rr_?NumericQ] := Which[rr <= r0, -aLo, rr <= rcut, sLo'[rr], True, (tail[rr + 10^-7] - tail[rr - 10^-7])/(2 10^-7)];
  <|"Mu" -> mu, "Slope" -> aLo, "ReliableRadius" -> rcut, "f" -> prof, "df" -> dprof|>];

(* ---- trial data from a profile: nodes, p = Pi - f, p' ---- *)
Options[hedgehogTrialData] = {"Radius" -> Automatic, "CoreStep" -> 1/25, "Growth" -> 104/100, "Digits" -> 18};
hedgehogTrialData[f_, df_, mu_, OptionsPattern[]] := Module[
  {R = Replace[OptionValue["Radius"], Automatic :> If[mu == 0, 60, Min[40, Ceiling[15/mu]]]], h0 = OptionValue["CoreStep"],
   g = OptionValue["Growth"], dig = OptionValue["Digits"], nodes = {0}, h, rat, pv, dpv},
  h = h0;
  While[Last[nodes] < R, If[Last[nodes] >= 3, h = Min[h g, 1/2]]; AppendTo[nodes, Min[Last[nodes] + h, R]]];
  rat[x_] := Rationalize[SetPrecision[x, dig + 4], 10^-dig];
  pv = Join[{0}, rat[Pi - f[N[#, dig + 4]]] & /@ nodes[[2 ;; -2]]];
  dpv = rat[-df[N[#, dig + 4]]] & /@ Most[nodes];
  <|"Mu" -> mu, "r" -> nodes, "p" -> pv, "dp" -> dpv|>];

(* ---- rigorous upper bound for the energy of the trial, in units of 12 Pi^2 ---- *)
(* "Tolerance": target width of the enclosure of Integrate[e, {r, 0, R}] (shared among the cells by length) *)
Options[hedgehogUpperBound] = {WorkingPrecision -> 30, "Tolerance" -> 10^-10};
hedgehogUpperBound[trial_Association, muq_, OptionsPattern[]] := Module[
  {r = trial["r"], p = trial["p"], dp = trial["dp"], n, x, mu2 = muq^2, poly, q, maj, first, middle, last, fmax, df,
   prec = OptionValue[WorkingPrecision], tol = OptionValue["Tolerance"], total, rn1, R, t0 = AbsoluteTime[], piIv},
  n = Length[r] - 1;
  If[! (r[[1]] === 0 && p[[1]] === 0 && OrderedQ[r] && Length[p] == n && Length[dp] == n &&
      AllTrue[Join[r, p, dp], MatchQ[#, _Integer | _Rational] &]),
    Throw[failure["BadCandidate", "The hedgehog trial must have exact rational data with p(0) = 0."], "CertiField"]];
  (* first cell: polynomial majorant, integrated exactly *)
  poly = hermiteCubic[r[[{1, 2}]], p[[{1, 2}]], dp[[{1, 2}]], x];
  q = Cancel[poly/x];
  If[! PolynomialQ[q, x], Throw[failure["Internal", "p(0) != 0"], "CertiField"]];
  maj = x^2 D[poly, x]^2 + 2 poly^2 (1 + D[poly, x]^2) + q^4 x^2 + 4 mu2 x^2;
  first = Integrate[maj, {x, 0, r[[2]]}];
  (* middle cells *)
  middle = Total[Table[
    poly = hermiteCubic[r[[{k, k + 1}]], p[[{k, k + 1}]], dp[[{k, k + 1}]], x];
    rigorousIntegral[hedgehogDensityP[poly, D[poly, x], x, mu2], {x, r[[k]], r[[k + 1]]},
      WorkingPrecision -> prec, "Cells" -> 1, "Tolerance" -> tol (r[[k + 1]] - r[[k]])/Last[r]],
    {k, 2, n - 1}]];
  (* last cell: f = (Pi - p_{n-1}) (R - x)/(R - r_{n-1}) *)
  rn1 = r[[n]]; R = r[[n + 1]];
  fmax = Max[Abs[$piHi - p[[n]]], Abs[$piLo - p[[n]]]];
  df = fmax/(R - rn1);
  maj = With[{fx = fmax (R - x)/(R - rn1)}, x^2 df^2 + 2 fx^2 (1 + df^2) + fx^4/rn1^2 + mu2 x^2 fx^2];
  last = Integrate[maj, {x, rn1, R}];
  total = toIv[first + last, prec] + middle;
  piIv = Interval[{$piLo, $piHi}];
  <|"UpperBound" -> Max[total/(3 piIv)], "Enclosure" -> total/(3 piIv), "Cells" -> n,
    "Seconds" -> AbsoluteTime[] - t0, "Radius" -> R|>];

End[];
