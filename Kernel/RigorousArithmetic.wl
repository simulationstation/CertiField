(* ::Package:: *)

(* ::Title:: *)
(*CertiField: rigorous arithmetic*)

(* ::Text:: *)
(*Two layers, both used by the certificate verifiers.*)
(**)
(*1. Exact enclosures of scalars.  Rationals bracket sin and cos at rational points (alternating Taylor series*)
(*   with the remainder), pi (35 digits), and cube roots (checked by cubing exactly).  Doubles are rounded in*)
(*   the safe direction by exact comparison.*)
(**)
(*2. Outward-rounded machine arithmetic on packed arrays.  Every +, -, *, / result is widened by a relative*)
(*   2^-49 (eight units in the last place) plus an absolute 10^-300.  This encloses the exact real result as long*)
(*   as each operation is IEEE-754 with at most two roundings (x/y may be evaluated as x*y^-1).  Comparisons and*)
(*   selections are done with 0/1 masks, which are exact.  Packed arrays cannot hold infinities, so +/-10^300*)
(*   stands in for them; sentinels are only ever compared or selected, never multiplied.*)
(**)
(*Elementary functions from the math library are never trusted inside a certificate: square and cube roots are*)
(*verified by squaring/cubing with directed rounding.*)

Begin["CertiField`Private`"];

(* ------------------------------------------------------------------------------------------------ *)
(* exact scalars                                                                                     *)
(* ------------------------------------------------------------------------------------------------ *)

$piLo = 314159265358979323846264338327950288/10^35;   (* pi = 3.14159265358979323846264338327950288419... *)
$piHi = $piLo + 1/10^35;
$eps40 = 1/10^40;

exactOf[x_Real] := SetPrecision[x, Infinity];
exactOf[x_Integer | x_Rational] := x;

binadeExponent[x_Real] := Last[MantissaExponent[Abs[x], 2]] - 1;   (* 2^e <= |x| < 2^(e+1) *)

nextDown[x_Real] := Module[{xe = exactOf[x], e},
  Which[
    xe > 0, e = binadeExponent[x]; N[xe - If[xe == 2^e, 2^(e - 53), 2^(e - 52)]],
    xe < 0, e = binadeExponent[x]; N[xe - 2^(e - 52)],
    True, -$MinMachineNumber]];
nextUp[x_Real] := -nextDown[-x];

(* largest (smallest) machine number <= (>=) the exact number q *)
fLo[q_] := Module[{x = N[q] + 0.}, While[exactOf[x] > q, x = nextDown[x]]; x];
fHi[q_] := Module[{x = N[q] + 0.}, While[exactOf[x] < q, x = nextUp[x]]; x];

(* {lo, hi} rationals with lo <= Sin[x] <= hi, x rational in [0, 7]: alternating series with decreasing terms *)
sinBr[x_] := Module[{x2 = x x, term = x, s = x, k = 0, nxt, res = None},
  While[res === None,
    nxt = -term x2/((2 k + 2) (2 k + 3));
    If[(2 k + 2) (2 k + 3) > x2 && Abs[nxt] < $eps40,
      res = {Min[s, s + nxt], Max[s, s + nxt]},
      s += nxt; term = nxt; k++]];
  res];

cosBr[x_] := Module[{x2 = x x, term = 1, s = 1, k = 0, nxt, res = None},
  While[res === None,
    nxt = -term x2/((2 k + 1) (2 k + 2));
    If[(2 k + 1) (2 k + 2) > x2 && Abs[nxt] < $eps40,
      res = {Min[s, s + nxt], Max[s, s + nxt]},
      s += nxt; term = nxt; k++]];
  res];

(* machine numbers {a, b} with a^3 <= qlo and b^3 >= qhi, checked exactly *)
cubeRootBounds[qlo_, qhi_] := Module[{a = N[qlo^(1/3)], b = N[qhi^(1/3)]},
  While[exactOf[a]^3 > qlo, a = nextDown[a]];
  While[exactOf[b]^3 < qhi, b = nextUp[b]];
  {a, b}];

(* ------------------------------------------------------------------------------------------------ *)
(* outward-rounded machine arithmetic                                                                *)
(* ------------------------------------------------------------------------------------------------ *)

$rEps = 2.^-49;
$rTiny = 1.*^-300;
$big = 1.*^300;          (* finite stand-in for +Infinity *)

rUp[x_] := x + (Abs[x] $rEps + $rTiny);
rDn[x_] := x - (Abs[x] $rEps + $rTiny);

(* exact elementwise selections (0/1 masks; no rounding is involved) *)
vWhere[m_, x_, y_] := m x + (1 - m) y;
vMax[a_, b_] := vWhere[UnitStep[a - b], a, b];
vMin[a_, b_] := vWhere[UnitStep[b - a], a, b];
gtMask[x_, y_] := 1 - UnitStep[y - x];     (* x > y  *)
geMask[x_, y_] := UnitStep[x - y];         (* x >= y *)
leMask[x_, y_] := UnitStep[y - x];         (* x <= y *)
ltMask[x_, y_] := 1 - UnitStep[x - y];     (* x < y  *)

(* intervals are pairs {lo, hi} of machine numbers or packed arrays *)
ivPt[x_] := {x, x};
ivAdd[{a_, b_}, {c_, d_}] := {rDn[a + c], rUp[b + d]};
ivSub[{a_, b_}, {c_, d_}] := {rDn[a - d], rUp[b - c]};
ivMul[{a_, b_}, {c_, d_}] := Module[{p1 = a c, p2 = a d, p3 = b c, p4 = b d},
  {rDn[vMin[vMin[p1, p2], vMin[p3, p4]]], rUp[vMax[vMax[p1, p2], vMax[p3, p4]]]}];
(* division by an interval with strictly positive lower end *)
ivDivPos[{a_, b_}, {c_, d_}] := Module[{q1 = a/c, q2 = a/d, q3 = b/c, q4 = b/d},
  {rDn[vMin[vMin[q1, q2], vMin[q3, q4]]], rUp[vMax[vMax[q1, q2], vMax[q3, q4]]]}];

(* verified roots of nonnegative arrays: lower/upper bounds checked by squaring/cubing with directed rounding *)
fixDown[y0_, bad_] := Module[{y = y0, b = bad[y0], it = 0},
  While[Total[b] > 0,
    If[++it > 80, Throw[$Failed, "CertiFieldArithmetic"]];
    y = vWhere[b, y (1 - 2.^-48), y]; b = bad[y]];
  y];
fixUp[y0_, bad_] := Module[{y = y0, b = bad[y0], it = 0},
  While[Total[b] > 0,
    If[++it > 80, Throw[$Failed, "CertiFieldArithmetic"]];
    y = vWhere[b, y (1 + 2.^-48) + $rTiny, y]; b = bad[y]];
  y];

(* Below the thresholds the absolute margin 10^-300 dominates, so those entries get direct valid bounds instead:
   0 <= Sqrt[v] <= 10^-100 for 0 < v <= 10^-200, and 0 <= v^(1/3) <= 10^-83 for 0 < v <= 10^-250. *)
sqrtLo[v_] := vWhere[gtMask[v, 1.*^-200],
  fixDown[Sqrt[Ramp[v]] (1 - 2.^-50), Function[y, gtMask[v, 1.*^-200] gtMask[rUp[y y], v]]], 0.];
sqrtHi[v_] := vWhere[gtMask[v, 1.*^-200],
  fixUp[Sqrt[Ramp[v]] (1 + 2.^-50), Function[y, gtMask[v, 1.*^-200] ltMask[rDn[y y], v]]],
  vWhere[gtMask[v, 0.], 1.*^-100, 0.]];
cbrtLo[v_] := vWhere[gtMask[v, 1.*^-250],
  fixDown[Power[Ramp[v], 1/3.] (1 - 2.^-50), Function[y, gtMask[v, 1.*^-250] gtMask[rUp[rUp[y y] y], v]]], 0.];
cbrtHi[v_] := vWhere[gtMask[v, 1.*^-250],
  fixUp[Power[Ramp[v], 1/3.] (1 + 2.^-50), Function[y, gtMask[v, 1.*^-250] ltMask[rDn[rDn[y y] y], v]]],
  vWhere[gtMask[v, 0.], 1.*^-83, 0.]];

(* ------------------------------------------------------------------------------------------------ *)
(* array utilities (compiled to the Wolfram Virtual Machine; no C compiler needed)                    *)
(* ------------------------------------------------------------------------------------------------ *)

(* numpy.searchsorted(a, x, side='right'): number of entries of the sorted list a that are <= x *)
searchRight = Compile[{{a, _Real, 1}, {x, _Real, 1}},
  Module[{n = Length[a], lo = 0, hi = 0, mid = 0, xi = 0.},
    Table[
      xi = x[[i]]; lo = 0; hi = n;
      While[lo < hi, mid = Quotient[lo + hi, 2]; If[a[[mid + 1]] <= xi, lo = mid + 1, hi = mid]];
      lo, {i, Length[x]}]],
  RuntimeOptions -> "Speed"];

(* numpy.minimum.reduceat for consecutive segments of the given (positive) lengths *)
segmentMin = Compile[{{v, _Real, 1}, {cnt, _Integer, 1}},
  Module[{out = Table[0., {Length[cnt]}], p = 1, m = 0.},
    Do[
      m = v[[p]];
      Do[If[v[[q]] < m, m = v[[q]]], {q, p + 1, p + cnt[[k]] - 1}];
      out[[k]] = m; p += cnt[[k]], {k, Length[cnt]}];
    out],
  RuntimeOptions -> "Speed"];

(* numpy.repeat(range(n), cnt) *)
repeatIndex = Compile[{{cnt, _Integer, 1}},
  Module[{out = Table[0, {Total[cnt]}], p = 1},
    Do[Do[out[[p]] = k - 1; p++, {cnt[[k]]}], {k, Length[cnt]}];
    out],
  RuntimeOptions -> "Speed"];

cumMax = Compile[{{v, _Real, 1}},
  Module[{out = v, m = v[[1]]}, Do[If[v[[i]] > m, m = v[[i]]]; out[[i]] = m, {i, Length[v]}]; out],
  RuntimeOptions -> "Speed"];
suffixMin = Compile[{{v, _Real, 1}},
  Module[{out = v, m = v[[-1]]}, Do[If[v[[i]] < m, m = v[[i]]]; out[[i]] = m, {i, Length[v], 1, -1}]; out],
  RuntimeOptions -> "Speed"];

End[];
