(* ::Package:: *)

(* ::Title:: *)
(*CertiField: rigorous one-dimensional quadrature*)

(* ::Text:: *)
(*Enclosure of Integrate[f[x], {x, a, b}] using Interval arithmetic at a fixed working precision.  On a cell [l, r]*)
(*with h = r - l and m its midpoint, both of the following contain the integral:*)
(*  first order :  h f([l, r])*)
(*  Simpson     :  h/6 (f(l) + 4 f(m) + f(r)) - h^5/2880 f''''([l, r])      (f in C^4 on the closed cell)*)
(*The cell enclosure is their intersection.  Cells that are too wide are bisected.  Wherever an interval*)
(*evaluation is unbounded (a singularity, or a fourth derivative that cannot be enclosed), only the other*)
(*enclosure is used.  If neither bounds the cell, the result is unbounded and the caller reports Inconclusive.*)
(*Interval evaluation uses the Wolfram Language's outward-rounded Interval arithmetic.  Endpoints are converted*)
(*with N[Interval[{q, q}], prec], which rounds outward.*)

Begin["CertiField`Private`"];

ivFinite[iv_Interval] := MatchQ[iv, Interval[{_?NumericQ, _?NumericQ}]] && FreeQ[iv, DirectedInfinity | Indeterminate];
ivFinite[_] := False;
ivWidth[Interval[{a_, b_}]] := b - a;
ivWidth[_] := Infinity;

toIv[q_, prec_] := N[Interval[{q, q}], prec];
toIv[{l_, r_}, prec_] := N[Interval[{l, r}], prec];

(* An evaluation result as an interval: exact real numbers become (outward-rounded) point intervals; anything that
   is not a single finite real interval becomes the whole line.  (Symbolic derivatives can simplify to plain numbers.) *)
asIv[x_Interval, _] := If[ivFinite[x], x, Interval[{-Infinity, Infinity}]];
asIv[x_?NumericQ, prec_] := If[Im[x] == 0, toIv[x, prec], Interval[{-Infinity, Infinity}]];
asIv[_, _] := Interval[{-Infinity, Infinity}];

(* f, f4: pure functions of one argument (the integrand and its fourth derivative).
   smoothQ: a predicate on the cell [l, r]; the Simpson enclosure (which needs f in C^4 on the closed cell) is used
   only where it returns True.  The first-order enclosure needs only continuity. *)
cellEnclosure[f_, f4_, l_, r_, prec_, smoothQ_ : (True &)] := Module[{h = r - l, m = (l + r)/2, hI, fo, si, rem, fl, fm, fr, enc},
  hI = toIv[h, prec];
  fo = hI asIv[Quiet[f[toIv[{l, r}, prec]]], prec];
  si = If[TrueQ[smoothQ[l, r]],
    {fl, fm, fr} = asIv[Quiet[f[toIv[#, prec]]], prec] & /@ {l, m, r};
    rem = asIv[Quiet[f4[toIv[{l, r}, prec]]], prec];
    If[AllTrue[{fl, fm, fr, rem}, ivFinite], hI/6 (fl + 4 fm + fr) - toIv[h^5/2880, prec] rem, Interval[{-Infinity, Infinity}]],
    Interval[{-Infinity, Infinity}]];
  enc = Which[
    ivFinite[fo] && ivFinite[si], IntervalIntersection[fo, si],
    ivFinite[si], si,
    ivFinite[fo], fo,
    True, Interval[{-Infinity, Infinity}]];
  If[enc === Interval[], Throw[failure["Internal", "Empty intersection in quadrature."], "CertiField"]];
  enc];

Options[rigorousIntegralFn] = {WorkingPrecision -> 30, "Cells" -> 64, "Tolerance" -> 10^-12, "MaxDepth" -> 48,
  "MaxCells" -> 20000, "Smooth" -> (True &)};

(* f, f4: pure functions (integrand and its fourth derivative); {a, b}: exact (rational) endpoints, a < b.
   Globally adaptive: while the total width exceeds the tolerance, every cell wider than tolerance/(2 #cells) is
   bisected.  Cells near a singularity (where only the first-order enclosure is available) are refined
   geometrically, so they do not dominate the cost.  The result is a rigorous enclosure whatever the stopping point. *)
rigorousIntegralFn[f_, f4_, {a_, b_}, OptionsPattern[]] := Module[
  {prec = OptionValue[WorkingPrecision], n = OptionValue["Cells"], tol = OptionValue["Tolerance"],
   maxDepth = OptionValue["MaxDepth"], maxCells = OptionValue["MaxCells"], smooth = OptionValue["Smooth"],
   cells, widths, thr, split, enc, it = 0},
  enc[l_, r_] := cellEnclosure[f, f4, l, r, prec, smooth];
  cells = {#[[1]], #[[2]], 0, enc[#[[1]], #[[2]]]} & /@ Partition[a + (b - a) Range[0, n]/n, 2, 1];
  While[True,
    widths = ivWidth /@ cells[[All, 4]];
    If[Total[widths] <= tol || Length[cells] >= maxCells || ++it > 4 maxDepth, Break[]];
    thr = tol/(2 Length[cells]);
    split = MapThread[#2 > thr && #1[[3]] < maxDepth &, {cells, widths}];
    If[! MemberQ[split, True], Break[]];
    cells = Join @@ MapThread[
      If[#2, With[{l = #1[[1]], r = #1[[2]], d = #1[[3]] + 1, mid = (#1[[1]] + #1[[2]])/2},
          {{l, mid, d, enc[l, mid]}, {mid, r, d, enc[mid, r]}}], {#1}] &, {cells, split}]];
  Total[cells[[All, 4]]]];

Options[rigorousIntegral] = Options[rigorousIntegralFn];
(* expr: an expression in the symbol x *)
rigorousIntegral[expr_, {x_Symbol, a_, b_}, opts : OptionsPattern[]] :=
  rigorousIntegralFn[Function @@ {x, expr}, Function @@ {x, D[expr, {x, 4}]}, {a, b}, opts];

(* ---- smoothness: expressions built only from functions that are analytic wherever they are finite and real ---- *)
$nonSmoothHeads = Abs | RealAbs | Sign | Max | Min | Clip | Ramp | Floor | Ceiling | Round | IntegerPart | FractionalPart |
  Mod | Quotient | UnitStep | HeavisideTheta | Piecewise | If | Which | Switch | Boole | KroneckerDelta | DiracDelta |
  Re | Im | Arg | Conjugate | Interpolation | InterpolatingFunction | Unitize;
analyticExpressionQ[expr_] := FreeQ[expr, $nonSmoothHeads] && FreeQ[expr, _Function | _Derivative];

(* helpers for safe display / export of bounds *)
(* exact directed rounding to a number of significant digits (never rounds a bound in the unsafe direction) *)
exactNumber[q_Real] := SetPrecision[q, Infinity];
exactNumber[q_] := q;
lowerDecimal[q_, digits_] := Module[{e = exactNumber[q], s},
  If[e == 0, 0, s = 10^(digits - 1 - Floor[Log10[Abs[N[e]]]]); Floor[e s]/s]];
upperDecimal[q_, digits_] := Module[{e = exactNumber[q], s},
  If[e == 0, 0, s = 10^(digits - 1 - Floor[Log10[Abs[N[e]]]]); Ceiling[e s]/s]];

End[];
