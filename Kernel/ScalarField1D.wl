(* ::Package:: *)

(* ::Title:: *)
(*CertiField: scalar field theories in one dimension (kinks, solitons, domain walls)*)

(* ::Text:: *)
(*Model.  E[u] = Integrate[k u'(x)^2 + V(u(x)), {x, a, b}], k > 0, a < b (either end may be infinite),*)
(*u(a) = uA, u(b) = uB (limits at infinite ends; there V must vanish).*)
(*Admissible class: absolutely continuous u with these boundary values and finite energy.*)
(**)
(*Lower bound (calibration).  If V + c >= 0 on the whole real line then pointwise*)
(*    k u'^2 + V(u) + c >= 2 Sqrt[k (V(u) + c)] |u'|,*)
(*and integrating gives, for every admissible u,*)
(*    E[u] >= LB(c) := 2 Sqrt[k] |Integrate[Sqrt[V(s) + c], {s, uA, uB}]| - c (b - a)*)
(*(c = 0 on infinite domains).  Equality holds for a monotone u with k u'^2 = V(u) + c, so maximising over c makes*)
(*the bound exact whenever the minimiser is monotone.  c is chosen numerically; the bound is then certified at that*)
(*exact rational c.*)
(**)
(*Upper bound.  The candidate is replaced by its cubic Hermite interpolant on a rational grid, with the boundary*)
(*values imposed exactly and continued by the vacua outside a finite window on infinite domains.  This competitor*)
(*is admissible.  Its energy is computed exactly (rational arithmetic) when V is a polynomial, and with rigorous*)
(*quadrature otherwise.*)

Begin["CertiField`Private`"];

(* ---- model normalisation ---- *)
(* Inexact numbers in the model are replaced by the exact rationals they denote (Rationalize[x, 0]), so the
   certificate is a statement about exactly that model. *)
scalarField1DSpec[params_Association] := Module[{V, k, dom, bv, s, vexpr},
  V = Lookup[params, "Potential", Missing[]];
  k = Lookup[params, "Stiffness", 1/2];
  dom = Lookup[params, "Domain", {-Infinity, Infinity}];
  bv = Lookup[params, "BoundaryValues", Missing[]];
  If[! MatchQ[V, _Function | _Symbol], Throw[failure["BadModel", "\"Potential\" must be a pure function of one variable."], "CertiField"]];
  If[! MatchQ[bv, {_?NumericQ, _?NumericQ}], Throw[failure["BadModel", "\"BoundaryValues\" must be a pair of numbers {uA, uB}."], "CertiField"]];
  If[! (MatchQ[dom, {_, _}] && (NumericQ[dom[[1]]] || dom[[1]] === -Infinity) && (NumericQ[dom[[2]]] || dom[[2]] === Infinity) && dom[[1]] < dom[[2]]),
    Throw[failure["BadModel", "\"Domain\" must be {a, b} with a < b; use -Infinity/Infinity for the real line."], "CertiField"]];
  If[! (NumericQ[k] && k > 0), Throw[failure["BadModel", "\"Stiffness\" must be a positive number."], "CertiField"]];
  vexpr = exactify[V[s]];
  If[! (NumericQ[vexpr /. s -> 1/3] && Im[N[vexpr /. s -> 1/3]] == 0),
    Throw[failure["BadModel", "\"Potential\" must be a real function of one real variable with no free parameters."], "CertiField"]];
  If[! analyticExpressionQ[vexpr],
    Throw[failure["BadModel", "\"Potential\" must be built from smooth functions (no Abs, Max, Piecewise, interpolation, ...)."], "CertiField"]];
  <|"V" -> Function @@ {s, vexpr}, "Vexpr" -> vexpr, "s" -> s, "k" -> exactify[k], "Domain" -> exactify /@ dom,
    "BoundaryValues" -> exactify /@ bv, "Polynomial" -> PolynomialQ[Expand[vexpr], s]|>];

exactify[expr_] := expr /. x_Real :> Rationalize[x, 0];

(* The smallest calibration constant c with V + c >= 0 on the whole real line, as an exact rational, together with a
   proof of V + c >= 0 by quantifier elimination (Resolve).  MinValue only proposes the value. *)
calibrationFloor[spec_] := Module[{m, c},
  m = Quiet@MinValue[spec["Vexpr"], spec["s"], Reals];
  If[! NumericQ[m] || ! TrueQ[Element[m, Reals]], m = If[provedNonnegative[spec, 0], 0, Missing["Unknown"]]];
  If[! NumericQ[m], Return[Missing["NotBoundedBelow"]]];
  c = Ceiling[-m 10^30]/10^30;                                              (* rational, c >= -m *)
  If[provedNonnegative[spec, c], c, Missing["NotProved"]]];

(* proof that V(s) + c >= 0 for all real s: quantifier elimination, directly or, for trigonometric polynomials,
   after Cos[s] -> cv, Sin[s] -> sv with cv^2 + sv^2 == 1 (a polynomial problem, decided by cylindrical decomposition).
   ForAll holds its arguments, so the expressions are injected with With. *)
provedNonnegative[spec_, c_] := Module[{x = spec["s"], e = spec["Vexpr"] + c, cv, sv, te},
  TrueQ[With[{xx = x, ee = e}, Quiet@Resolve[ForAll[xx, Element[xx, Reals], ee >= 0], Reals]]] ||
  (te = TrigExpand[e] /. {Cos[x] -> cv, Sin[x] -> sv};
   FreeQ[te, x] && PolynomialQ[te, {cv, sv}] &&
     TrueQ[With[{tt = te}, Quiet@Resolve[ForAll[{cv, sv}, cv^2 + sv^2 == 1, tt >= 0], Reals]]])];

(* ---- the lower bound LB(c) ---- *)
(* Simpson enclosures are used only on cells where V + c > 0 (there Sqrt[V + c] is analytic); elsewhere first order. *)
scalarLB[spec_, c_, prec_] := Module[{uA, uB, lo, hi, V = spec["V"], f, f4, s, enc, L, smooth},
  {uA, uB} = spec["BoundaryValues"];
  {lo, hi} = Sort[{uA, uB}];
  f = Function[u, Sqrt[nonNeg[V[u] + c]]];
  f4 = Function @@ {s, D[Sqrt[V[s] + c], {s, 4}]};
  smooth = Function[{l, r}, With[{g = asIv[Quiet[V[toIv[{l, r}, prec]]], prec] + toIv[c, prec]}, ivFinite[g] && Min[g] > 0]];
  enc = If[lo == hi, Interval[{0, 0}],
    rigorousIntegralFn[f, f4, {lo, hi}, WorkingPrecision -> prec, "Cells" -> 64, "Tolerance" -> 10^-13, "Smooth" -> smooth]];
  L = Subtract @@ Reverse[spec["Domain"]];
  If[c == 0, 2 Sqrt[toIv[spec["k"], prec]] enc, 2 Sqrt[toIv[spec["k"], prec]] enc - toIv[c L, prec]]];

nonNeg[x_Interval] := IntervalIntersection[x, Interval[{0, Infinity}]] /. Interval[] -> Interval[{0, 0}];
nonNeg[x_?NumericQ] := Max[x, 0];

(* numerical estimate of the optimal calibration constant (heuristic: any admissible c gives a valid bound) *)
(* LB(c) is concave in c (its second derivative is -Sqrt[k]/2 Integrate[(V + c)^(-3/2)] < 0) and decreasing for
   c > cmin + k (uB - uA)^2/L^2, so a golden-section search on that bracket finds its maximum. *)
optimalCalibration[spec_, cmin_] := Module[{uA, uB, L, k, lb, a, b, x1, x2, f1, f2, g = (Sqrt[5.] - 1)/2, sv},
  {uA, uB} = N[spec["BoundaryValues"]];
  L = N[Subtract @@ Reverse[spec["Domain"]]]; k = N[spec["k"]];
  lb[cc_?NumericQ] := 2 Sqrt[k] Abs[Quiet@NIntegrate[Sqrt[Max[spec["V"][sv] + cc, 0]], {sv, uA, uB}, AccuracyGoal -> 10]] - cc L;
  a = N[cmin]; b = a + 2 k (uB - uA)^2/L^2 + 1;
  x1 = b - g (b - a); x2 = a + g (b - a); f1 = lb[x1]; f2 = lb[x2];
  Do[If[f1 < f2,
      a = x1; x1 = x2; f1 = f2; x2 = a + g (b - a); f2 = lb[x2],
      b = x2; x2 = x1; f2 = f1; x1 = b - g (b - a); f1 = lb[x1]], {40}];
  If[NumericQ[a + b], Max[cmin, Rationalize[(a + b)/2, 10^-12]], cmin]];

(* ---- the competitor ---- *)
Options[scalarCompetitor] = {"Nodes" -> 2000, "Truncation" -> Automatic, "Digits" -> 20};

candidateFunctions[cand_] := Which[
  MatchQ[cand, {{_?NumericQ, _?NumericQ} ..}],
    With[{ip = Interpolation[N[cand], InterpolationOrder -> 3]}, {ip, ip', ip["Domain"][[1]]}],
  Head[cand] === InterpolatingFunction,
    {cand, Derivative[1][cand], cand["Domain"][[1]]},
  MatchQ[Head[cand], Function | Symbol | Composition],
    {cand, Derivative[1][cand], {-Infinity, Infinity}},
  True, Throw[failure["BadCandidate", "The candidate must be a function of x, an InterpolatingFunction, or a list of {x, u} pairs."], "CertiField"]];

(* The competitor: cubic Hermite interpolation of the candidate on a uniform rational grid, values and slopes rounded
   to multiples of 10^-digits (so all cells share one denominator), boundary values imposed exactly, constant
   (= the vacuum) outside the window on infinite domains.  Whatever the candidate, this is an admissible function;
   the candidate only affects how tight the upper bound is. *)
scalarCompetitor[spec_, cand_, OptionsPattern[]] := Module[
  {n = OptionValue["Nodes"], dig = OptionValue["Digits"], u, du, cdom, a, b, uA, uB, xs, us, ds, rat, L, tol = 10^-12, val, near},
  {u, du, cdom} = candidateFunctions[cand];
  {a, b} = spec["Domain"]; {uA, uB} = spec["BoundaryValues"];
  rat[x_] := If[NumericQ[x] && Im[x] == 0, Round[SetPrecision[x, dig + 10] 10^dig]/10^dig, $Failed];
  val[x_] := Quiet@Check[N[u[N[x, dig + 10]], dig + 10], $Failed];
  near[x_, target_] := With[{y = val[x]}, NumericQ[y] && Abs[y - target] < tol];
  (* finite window for infinite ends: where the candidate has reached the vacua (or the end of its own domain) *)
  If[a === -Infinity || b === Infinity,
    L = OptionValue["Truncation"];
    If[L === Automatic,
      L = 1;
      While[L < 2^14 && ((a === -Infinity && -L > cdom[[1]] && ! near[-L, uA]) || (b === Infinity && L < cdom[[2]] && ! near[L, uB])), L *= 2]];
    If[a === -Infinity, a = Max[-L, Ceiling[cdom[[1]]]]];
    If[b === Infinity, b = Min[L, Floor[cdom[[2]]]]]];
  If[! (a < b), Throw[failure["BadCandidate", "The candidate's domain does not overlap the model's domain."], "CertiField"]];
  xs = a + (b - a) Range[0, n]/n;
  us = rat[val[#]] & /@ xs;
  ds = rat[Quiet@Check[N[du[N[#, dig + 10]], dig + 10], $Failed]] & /@ xs;
  If[! AllTrue[Join[us, ds], MatchQ[#, _Integer | _Rational] &],
    Throw[failure["BadCandidate", "The candidate could not be evaluated to real numbers on the domain."], "CertiField"]];
  (* impose the boundary values exactly; flat ends on infinite domains *)
  us[[1]] = uA; us[[-1]] = uB;
  If[spec["Domain"][[1]] === -Infinity, ds[[1]] = 0];
  If[spec["Domain"][[2]] === Infinity, ds[[-1]] = 0];
  <|"x" -> xs, "u" -> us, "du" -> ds, "Window" -> {a, b}|>];

(* ---- energy of the competitor ---- *)

(* exact kinetic energy of a cubic Hermite cell of width h:  Integrate[p'(x)^2, {x, x0, x0 + h}] *)
hermiteKinetic[h_, u0_, u1_, d0_, d1_] := With[{dd = u1 - u0, m0 = h d0, m1 = h d1},
  (6/5 dd^2 - dd (m0 + m1)/5 + 2/15 (m0^2 + m1^2) - m0 m1/15)/h];

(* ascending t-coefficients of the Hermite cubic on [0, 1] with values u0, u1 and scaled slopes m0 = h d0, m1 = h d1 *)
hermiteCoefficients[u0_, u1_, m0_, m1_] := {u0, m0, 3 (u1 - u0) - 2 m0 - m1, 2 (u0 - u1) + m0 + m1};

(* polynomial V:  Integrate[V(p(t)), {t, 0, 1}]  as an exact polynomial in the cell data (computed once per model) *)
hermitePotentialForm[spec_] := Module[{t, a, b, c, d},
  Function @@ {{a, b, c, d}, Expand[Integrate[Expand[spec["V"][hermiteCoefficients[a, b, c, d] . t^Range[0, 3]]], {t, 0, 1}]]}];

(* other V:  Taylor expansion about p(1/2) to order K with Lagrange remainder,
     V(p(t)) = Sum[V^(k)(um) q(t)^k/k!, {k, 0, K - 1}] + V^(K)(xi) q(t)^K/K!,   q = p - um,
   the t-integrals of q^k are exact rationals and |V^(K)| is enclosed on an exact enclosure of the range of p.
   Cells where the remainder is too wide are bisected in t. *)
$taylorOrder = 8;
polyMul[a_, b_] := ListConvolve[a, b, {1, -1}, 0];
polyInt01[c_] := Total[c/Range[Length[c]]];
polyRange01[c_] := Fold[#1 Interval[{0, 1}] + #2 &, 0, Reverse[c]];          (* exact interval containing p([0, 1]) *)
shiftHalf[c_] := Table[Sum[c[[j + 1]] Binomial[j, i]/2^j, {j, i, Length[c] - 1}], {i, 0, Length[c] - 1}];   (* p(1/2 + s/2) *)

cubicPotentialIv[vd_, cf_, tol_, prec_, depth_] := Module[{K = $taylorOrder, um, q, ints, rng, umI, mk, enc},
  um = cf . (1/2)^Range[0, 3];
  q = cf - {um, 0, 0, 0};
  ints = polyInt01 /@ NestList[polyMul[#, q] &, {1}, K];
  rng = um + polyRange01[q];
  umI = toIv[um, prec];
  mk = Quiet[vd[[K + 1]][N[rng, prec]]];
  enc = If[! ivFinite[mk], Interval[{-Infinity, Infinity}],
    Sum[vd[[k + 1]][umI] toIv[ints[[k + 1]]/k!, prec], {k, 0, K - 1}] + toIv[ints[[K + 1]]/K!, prec] Interval[{-1, 1}] Max[Abs[mk]]];
  If[depth < 12 && (! ivFinite[enc] || ivWidth[enc] > tol),
    (cubicPotentialIv[vd, cf/2^Range[0, 3], tol, prec, depth + 1] + cubicPotentialIv[vd, shiftHalf[cf], tol, prec, depth + 1])/2,
    enc]];

Options[competitorEnergy] = {"Tolerance" -> 10^-12};
(* exact (polynomial V) or rigorous (other V) energy of the competitor; returns a rational or an Interval *)
competitorEnergy[spec_, comp_, prec_, OptionsPattern[]] := Module[
  {xs = comp["x"], us = comp["u"], ds = comp["du"], hs, kin, pot, form, vd, cfs, tolCell},
  hs = Differences[xs];
  kin = spec["k"] Total[hermiteKinetic[hs, Most[us], Rest[us], Most[ds], Rest[ds]]];
  If[spec["Polynomial"],
    form = hermitePotentialForm[spec];
    kin + Total[hs form[Most[us], Rest[us], hs Most[ds], hs Rest[ds]]],
    vd = Table[Function @@ {spec["s"], D[spec["Vexpr"], {spec["s"], k}]}, {k, 0, $taylorOrder}];
    cfs = Transpose[hermiteCoefficients[Most[us], Rest[us], hs Most[ds], hs Rest[ds]]];
    tolCell = OptionValue["Tolerance"]/(Last[xs] - First[xs]);
    pot = Total[MapThread[toIv[#1, prec] cubicPotentialIv[vd, #2, tolCell, prec, 0] &, {hs, cfs}]];
    toIv[kin, prec] + pot]];

(* ---- the certificate ---- *)
certifyScalarField1D[params_Association, cand_, opts_List] := Module[
  {spec = scalarField1DSpec[params], prec = OptionValue[CertifiedVariationalBound, opts, WorkingPrecision], cmin, c, lbIv,
   comp, eu, lb, ub, nodes, t0 = AbsoluteTime[], infL, infR},
  nodes = Replace[OptionValue[CertifiedVariationalBound, opts, "Nodes"], Automatic -> 2000];
  infL = spec["Domain"][[1]] === -Infinity; infR = spec["Domain"][[2]] === Infinity;
  If[(infL && ! PossibleZeroQ[spec["V"][spec["BoundaryValues"][[1]]]]) || (infR && ! PossibleZeroQ[spec["V"][spec["BoundaryValues"][[2]]]]),
    Throw[failure["BadModel", "At an infinite end the boundary value must be a zero of the potential (a vacuum)."], "CertiField"]];
  cmin = calibrationFloor[spec];
  If[! NumericQ[cmin],
    Throw[failure["Inconclusive", "Could not prove that the potential is bounded below (V + c >= 0 on the real line); no lower bound is available."], "CertiField"]];
  If[(infL || infR) && cmin > 0, Throw[failure["Inconclusive", "On an infinite domain the potential must be nonnegative."], "CertiField"]];
  c = If[infL || infR, 0, Max[optimalCalibration[spec, cmin], cmin]];
  lbIv = scalarLB[spec, c, prec];
  comp = scalarCompetitor[spec, cand, "Nodes" -> nodes];
  eu = competitorEnergy[spec, comp, prec];
  lb = Min[lbIv]; ub = If[Head[eu] === Interval, Max[eu], eu];
  <|"LowerBound" -> lb, "UpperBound" -> ub, "Seconds" -> AbsoluteTime[] - t0,
    "Calibration" -> c, "Competitor" -> comp, "Exact" -> spec["Polynomial"], "Spec" -> spec|>];

End[];
