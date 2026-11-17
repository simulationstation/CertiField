(* ::Package:: *)

(* ::Title:: *)
(*CertiField: the CertifiedBound object*)

Begin["CertiField`Private`"];

(* ---- derived quantities (rigorous: computed exactly from the stored bounds, then rounded up) ---- *)
gapOf[a_Association] := If[NumericQ[a["LowerBound"]] && NumericQ[a["UpperBound"]],
  exactNumber[a["UpperBound"]] - exactNumber[a["LowerBound"]], Missing["NotAvailable"]];
relGapOf[a_Association] := With[{g = gapOf[a]}, If[NumericQ[g] && a["UpperBound"] > 0, g/exactNumber[a["UpperBound"]], Missing["NotAvailable"]]];

$derivedProperties = {"Gap", "RelativeGap", "Report", "Association", "Properties", "Summary"};

CertifiedBound[a_Association]["Properties"] := Sort@DeleteDuplicates@Join[Keys[a], $derivedProperties];
CertifiedBound[a_Association]["Association"] := a;
CertifiedBound[a_Association]["Gap"] := N[upperDecimal[gapOf[a], 6]] /; NumericQ[gapOf[a]];
CertifiedBound[a_Association]["RelativeGap"] := N[upperDecimal[relGapOf[a], 6]] /; NumericQ[relGapOf[a]];
CertifiedBound[a_Association]["Report"] := certReport[a];
CertifiedBound[a_Association]["Summary"] := certSummaryString[a];
CertifiedBound[a_Association][key_String] := Lookup[a, key, Missing["NotAvailable", key]];
CertifiedBound[a_Association][keys_List] := CertifiedBound[a] /@ keys;

(* ---- display ---- *)
$statusColor = <|"Certified" -> RGBColor[0.13, 0.55, 0.13], "Inconclusive" -> RGBColor[0.85, 0.45, 0.0]|>;
statusStyle[s_String] := Style[s, Bold, Lookup[$statusColor, s, Gray]];

(* ---- safe decimal strings: lower bounds rounded down, upper bounds rounded up, exactly ---- *)
(* number of significant digits that makes the gap visible (two more than needed to separate the bounds) *)
displayDigits[a_Association] := With[{g = relGapOf[a]},
  If[NumericQ[g] && g > 0, Clip[Ceiling[-Log10[N[g]]] + 2, {6, 16}], 10]];

(* q: exact number that is a multiple of 10^-j (as returned by lowerDecimal/upperDecimal) *)
decimalString[q_, j_Integer] := Module[{n = Round[q 10^Max[j, 0]], ds, jj = Max[j, 0]},
  ds = IntegerString[Abs[n], 10, jj + 1];
  If[n < 0, "-", ""] <> StringDrop[ds, -jj] <> If[jj > 0, "." <> StringTake[ds, -jj], ""]];
decimalPlaces[q_, digits_] := If[exactNumber[q] == 0, 0, digits - 1 - Floor[Log10[Abs[N[exactNumber[q]]]]]];
lowerString[q_, digits_] := decimalString[lowerDecimal[q, digits], decimalPlaces[q, digits]];
upperString[q_, digits_] := decimalString[upperDecimal[q, digits], decimalPlaces[q, digits]];

(* a relative gap as a percentage, rounded up *)
percentString[g_] := Module[{p = 100 upperDecimal[g, 3], e},
  If[p == 0, "0%",
    e = Floor[Log10[N[p]]];
    If[e >= -3,
      decimalString[upperDecimal[p, 3], Max[0, 2 - e]] <> "%",
      decimalString[upperDecimal[p, 3]/10^e, 2] <> "\[Times]10^" <> ToString[e] <> " %"]]];

boundsString[a_Association] := With[{d = displayDigits[a]},
  lowerString[a["LowerBound"], d] <> " \[LessEqual] inf E \[LessEqual] " <> upperString[a["UpperBound"], d]];

boundsDisplay[a_Association] := If[NumericQ[a["LowerBound"]] && NumericQ[a["UpperBound"]],
  boundsString[a], Row[{a["LowerBound"], " \[LessEqual] inf E \[LessEqual] ", a["UpperBound"]}]];

gapDisplay[a_Association] := With[{g = relGapOf[a]}, If[NumericQ[g], percentString[g], "\[LongDash]"]];

certIcon[status_] := Graphics[{
    {GrayLevel[0.7], AbsoluteThickness[1], Line[{{0, 0}, {1, 0}}]},
    {Lookup[$statusColor, status, Gray], AbsoluteThickness[4], CapForm["Butt"], Line[{{0.55, 0}, {0.82, 0}}]},
    {Lookup[$statusColor, status, Gray], AbsoluteThickness[1.5],
     Line[{{0.55, -0.18}, {0.55, 0.18}}], Line[{{0.82, -0.18}, {0.82, 0.18}}]},
    {GrayLevel[0.35], Disk[{0.78, 0}, 0.045]}},
  ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}],
  PlotRange -> {{-0.05, 1.05}, {-0.35, 0.35}}, Background -> GrayLevel[0.97], AspectRatio -> 1];

CertifiedBound /: MakeBoxes[obj : CertifiedBound[a_Association], fmt_] /; KeyExistsQ[a, "Status"] :=
  BoxForm`ArrangeSummaryBox[CertifiedBound, obj, certIcon[a["Status"]],
    {
      {BoxForm`SummaryItem[{"Model: ", a["Model"]}], BoxForm`SummaryItem[{"Status: ", statusStyle[a["Status"]]}]},
      {BoxForm`SummaryItem[{"Bounds: ", boundsDisplay[a]}], BoxForm`SummaryItem[{"Relative gap: ", gapDisplay[a]}]}
    },
    {
      {BoxForm`SummaryItem[{"Theorem: ", Lookup[a, "Theorem", "\[LongDash]"]}]},
      {BoxForm`SummaryItem[{"Admissible class: ", Lookup[a, "AdmissibleClass", "\[LongDash]"]}]},
      {BoxForm`SummaryItem[{"Lower bound by: ", Lookup[a, "LowerBoundMethod", "\[LongDash]"]}]},
      {BoxForm`SummaryItem[{"Upper bound by: ", Lookup[a, "UpperBoundMethod", "\[LongDash]"]}]},
      {BoxForm`SummaryItem[{"Verified in: ", Row[{NumberForm[Lookup[a, "Seconds", 0.], {6, 1}], " s, CertiField ", Lookup[a, "Version", ""]}]}]}
    },
    fmt, "Interpretable" -> Automatic];

certReport[a_Association] := Grid[
  Select[{
    {"Model", a["Model"]},
    {"Parameters", Lookup[a, "Parameters", <||>]},
    {"Status", statusStyle[a["Status"]]},
    {"Lower bound", lowerString[a["LowerBound"], Max[displayDigits[a], 10]]},
    {"Upper bound", upperString[a["UpperBound"], Max[displayDigits[a], 10]]},
    {"Relative gap", gapDisplay[a]},
    {"Units", Lookup[a, "Units", None]},
    {"Theorem", Lookup[a, "Theorem", None]},
    {"Admissible class", Lookup[a, "AdmissibleClass", None]},
    {"Lower bound method", Lookup[a, "LowerBoundMethod", None]},
    {"Upper bound method", Lookup[a, "UpperBoundMethod", None]},
    {"Trusted base", Lookup[a, "TrustedBase", None]},
    {"Not claimed", Lookup[a, "NotClaimed", None]},
    {"References", Lookup[a, "References", None]},
    {"Verification time", Row[{NumberForm[Lookup[a, "Seconds", 0.], {6, 1}], " s"}]},
    {"CertiField version", Lookup[a, "Version", None]}},
   #[[2]] =!= None &],
  Alignment -> {{Right, Left}, Top}, Spacings -> {2, 0.8}, Dividers -> {False, All},
  FrameStyle -> GrayLevel[0.85], ItemStyle -> {{Directive[Bold, GrayLevel[0.3]], Automatic}},
  ItemSize -> {{Automatic, 60}}];

certSummaryString[a_Association] := With[{d = displayDigits[a]}, StringJoin[
  a["Model"], ": ", a["Status"], "  ",
  lowerString[a["LowerBound"], d], " <= inf E <= ", upperString[a["UpperBound"], d],
  If[NumericQ[relGapOf[a]], "  (relative gap <= " <> StringReplace[percentString[relGapOf[a]], {"\[Times]" -> "x"}] <> ")", ""]]];

End[];
