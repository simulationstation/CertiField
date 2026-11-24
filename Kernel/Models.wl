(* ::Package:: *)

(* ::Title:: *)
(*CertiField: model registry and the public entry points*)

Begin["CertiField`Private`"];

$CertiFieldModels = {"ScalarField1D", "Phi4Kink", "SineGordonKink", "MagneticDomainWall", "ScalarFieldND", "LinearElasticity", "SkyrmeDegreeOne"};

(* ---- model specifications ---- *)
CertiFieldModel[name_String, rules___Rule] := CertiFieldModel[name, Association[rules]];
CertiFieldModel[name_String, rules : {___Rule}] := CertiFieldModel[name, Association[rules]];

CertiFieldModel /: MakeBoxes[m : CertiFieldModel[name_String, p_Association], fmt_] :=
  BoxForm`ArrangeSummaryBox[CertiFieldModel, m, None,
    {{BoxForm`SummaryItem[{"Model: ", name}]}, {BoxForm`SummaryItem[{"Parameters: ", Short[p, 2]}]}}, {}, fmt,
    "Interpretable" -> Automatic];

(* named presets reduce to ScalarField1D *)
presetSpec["Phi4Kink", p_] := With[{lam = Lookup[p, "Lambda", 1], v = Lookup[p, "Vacuum", 1]},
  <|"Potential" -> Function[u, lam/4 (u^2 - v^2)^2], "Stiffness" -> 1/2, "Domain" -> Lookup[p, "Domain", {-Infinity, Infinity}],
    "BoundaryValues" -> Lookup[p, "BoundaryValues", {-v, v}]|>];
presetSpec["SineGordonKink", p_] := With[{m = Lookup[p, "Mass", 1]},
  <|"Potential" -> Function[u, m^2 (1 - Cos[u])], "Stiffness" -> 1/2, "Domain" -> Lookup[p, "Domain", {-Infinity, Infinity}],
    "BoundaryValues" -> Lookup[p, "BoundaryValues", {0, 2 Pi}]|>];
presetSpec["MagneticDomainWall", p_] := With[{A = Lookup[p, "A", 1], K = Lookup[p, "K", 1]},
  <|"Potential" -> Function[t, K Sin[t]^2], "Stiffness" -> A, "Domain" -> Lookup[p, "Domain", {-Infinity, Infinity}],
    "BoundaryValues" -> Lookup[p, "BoundaryValues", {0, Pi}]|>];
presetSpec["ScalarField1D", p_] := p;

presetCandidate["Phi4Kink", p_] := With[{lam = Lookup[p, "Lambda", 1], v = Lookup[p, "Vacuum", 1]}, Function[x, v Tanh[v Sqrt[lam/2] x]]];
presetCandidate["SineGordonKink", p_] := With[{m = Lookup[p, "Mass", 1]}, Function[x, 4 ArcTan[Exp[m x]]]];
presetCandidate["MagneticDomainWall", p_] := With[{A = Lookup[p, "A", 1], K = Lookup[p, "K", 1]}, Function[x, 2 ArcTan[Exp[x Sqrt[K/A]]]]];

(* ---- model documentation ---- *)
CertiFieldModelInfo["ScalarField1D"] := <|
  "Name" -> "ScalarField1D",
  "Energy" -> "E[u] = \[Integral] (k u'(x)^2 + V(u(x))) dx over a < x < b",
  "Parameters" -> <|"Potential" -> "pure function V", "Stiffness" -> "k > 0 (default 1/2)", "Domain" -> "{a, b}, ends may be infinite",
    "BoundaryValues" -> "{uA, uB}; at an infinite end it must be a zero of V"|>,
  "AdmissibleClass" -> "absolutely continuous u with the boundary values (limits at infinite ends) and finite energy",
  "LowerBound" -> "calibration: E \[GreaterEqual] 2\[Sqrt]k |\[Integral] \[Sqrt](V+c) ds| - c (b - a), for every c with V + c \[GreaterEqual] 0 on \[DoubleStruckCapitalR]; exact when the minimiser is monotone",
  "UpperBound" -> "energy of the cubic Hermite interpolant of the candidate (exact rational arithmetic for polynomial V)",
  "Status" -> "Supported"|>;
CertiFieldModelInfo["Phi4Kink"] := <|"Name" -> "Phi4Kink", "Energy" -> "E[u] = \[Integral] (u'^2/2 + \[Lambda]/4 (u^2 - v^2)^2) dx", "Parameters" -> <|"Lambda" -> 1, "Vacuum" -> 1|>,
  "Reduces to" -> "ScalarField1D", "Minimum (exact)" -> "(2\[Sqrt]2/3) \[Sqrt]\[Lambda] v^3 on the line", "Status" -> "Supported"|>;
CertiFieldModelInfo["SineGordonKink"] := <|"Name" -> "SineGordonKink", "Energy" -> "E[u] = \[Integral] (u'^2/2 + m^2 (1 - cos u)) dx", "Parameters" -> <|"Mass" -> 1|>,
  "Reduces to" -> "ScalarField1D", "Minimum (exact)" -> "8 m on the line", "Status" -> "Supported"|>;
CertiFieldModelInfo["MagneticDomainWall"] := <|"Name" -> "MagneticDomainWall",
  "Energy" -> "\[Sigma][\[Theta]] = \[Integral] (A \[Theta]'^2 + K sin^2 \[Theta]) dx  (Bloch wall, uniaxial anisotropy; energy per unit wall area)",
  "Parameters" -> <|"A" -> "exchange stiffness", "K" -> "anisotropy constant"|>, "Reduces to" -> "ScalarField1D",
  "Minimum (exact)" -> "4 \[Sqrt](A K) on the line", "Status" -> "Supported"|>;
CertiFieldModelInfo["ScalarFieldND"] := <|
  "Name" -> "ScalarFieldND",
  "Energy" -> "E[u] = \[Integral] (k |\[Del]u|^2 + W(u)) dx over a box in two or three dimensions",
  "Parameters" -> <|"Potential" -> "pure function W, a polynomial with rational coefficients (need not be convex)",
    "Stiffness" -> "k > 0 (default 1/2)", "Domain" -> "{{a1, b1}, {a2, b2}} or {{a1, b1}, {a2, b2}, {a3, b3}}, rational",
    "Cells" -> "grid, an integer or one integer per direction (default 32 in 2D, 16 in 3D)",
    "BoundaryData" -> "a number, a function of the coordinates, or an array of nodal values; the certified datum is the piecewise multilinear interpolant of the nodal values on the boundary grid",
    "NaturalFaces" -> "list of {direction, \"Lower\" | \"Upper\"}: faces where u is free (default: none)",
    "Field" -> "optional: a calibration field to certify, or a directory written by the Python reference; it replaces the package's own choice of field"|>,
  "Options" -> <|"Method" -> "Automatic (convex mode if W is convex, calibration mode otherwise), \"Convex\", \"Calibration\", or {\"Calibration\", suboptions}",
    "Suboptions of calibration mode" -> <|"TRange" -> "{tlo, thi}, rational: the range of the B-splines in u (default: the hull of the Dirichlet data and of the critical points of W)",
      "KnotIntervals" -> "K on the final grid (default 16)", "Levels" -> "{{cells, K}, ...}, coarse to fine, ending on the grid of the model (default: from Cells and K)",
      "Scale" -> "multiplies the smoothing parameters (default: max(1, range of W on the t-range))",
      "Samples" -> "samples per knot interval in the optimiser (default 4)", "Subintervals" -> "subintervals per knot interval in the certified evaluation (default 16)",
      "Schedule" -> "<|\"First\" | \"Later\" | \"Last\" -> {{tau, iterations}, ...}|> for the first, the intermediate and the last level"|>,
    "RefineCandidate" -> "True: the given candidate is only the starting point of the package's own minimisation (default False: it is certified as given)",
    "Verbose" -> "True prints the progress of calibration mode"|>,
  "AdmissibleClass" -> "all u in H^1 of the box with the given trace on the Dirichlet faces",
  "LowerBound" -> "E \[GreaterEqual] \[Integral] S(x,g)\[CenterDot]\[Nu] - \[Integral] sup_t (div S + |\[PartialD]_t S|^2/(4k) - W) for any Raviart\[Dash]Thomas/B-spline field S, evaluated with outward-rounded arithmetic. Convex mode: the flux of a finite-volume solve (seconds, tight for convex W). Calibration mode: a field found by a multilevel optimiser (minutes, for W that is not convex)",
  "UpperBound" -> "energy of the piecewise multilinear function through the candidate's nodal values (rigorous enclosure)",
  "Reference" -> "docs/SCALAR_ND.md",
  "Status" -> "Supported"|>;
CertiFieldModelInfo["SkyrmeDegreeOne"] := <|
  "Name" -> "SkyrmeDegreeOne",
  "Energy" -> "E(U) = \[Integral] |dU|^2 + |\[CapitalLambda]^2 dU|^2 + 2\[Mu]^2 (1 - cos f) over \[DoubleStruckCapitalR]^3, f = dist(U, 1); reported as E/(12\[Pi]^2)",
  "Parameters" -> <|"PionMass" -> "a rational \[Mu] with 0 \[LessEqual] \[Mu] < 2 (default 1). Grids and trial maps are shipped for 0, 263/500, 13/20 and 1; for other masses they are computed on request"|>,
  "Options" -> <|"Resolution" -> "\"Quick\", \"Standard\" or \"Publication\": the spacing of the grids of the level-set programme",
    "Method" -> "Automatic; \"Monotonicity\": the lower bound of the largest shipped mass not above \[Mu] (valid because the energy is non-decreasing in the mass; no grid is generated); \"GeneratedGrid\": a generated grid even where one is shipped"|>,
  "AdmissibleClass" -> "finite-energy maps U: \[DoubleStruckCapitalR]^3 \[RightArrow] S^3 of degree one with U - 1 \[Element] L^6 (the class in which minimisers exist)",
  "LowerBound" -> "level-set lower bound (coarea + isoperimetry + convexification; certified dynamic programme)",
  "UpperBound" -> "energy of an explicit hedgehog trial map (rigorous quadrature)",
  "Resolutions" -> Keys[$skyrmePresets],
  "Reference" -> "docs/SKYRME_ANY_MASS.md (masses other than the four shipped ones)",
  "Status" -> "Supported"|>;
CertiFieldModelInfo[other_] := Missing["UnknownModel", other];

(* ---- Skyrme: presets of the level-set certificate (grids shipped as data) ---- *)
$skyrmePresets = <|
  "Quick" -> <|"Table" -> {2000, 900, 1/10^10}|>,
  "Standard" -> <|"Table" -> {2000, 900, 1/10^10}|>,
  "Publication" -> <|"Table" -> {2000, 900, 1/10^10}|>|>;
$skyrmeMuNames = <|0 -> "0", 263/500 -> "263o500", 13/20 -> "13o20", 1 -> "1"|>;

skyrmeGridDir[muq_, preset_String] := FileNameJoin[{$CertiFieldDirectory, "Data", "Skyrme", "mu" <> $skyrmeMuNames[muq] <> "_" <> ToLowerCase[preset]}];
skyrmeTrialFile[muq_] := FileNameJoin[{$CertiFieldDirectory, "Data", "Skyrme", "hedgehog_mu" <> $skyrmeMuNames[muq] <> ".wl"}];

(* the certified minorant table is built once per kernel session (and cached on disk, re-certified on load) *)
(* The cache stores only the raw table; every entry is re-certified whenever it is loaded. *)
minorantTable[{nK_, nD_, kmin_}] := minorantTable[{nK, nD, kmin}] = Module[{cache, data, tab = $Failed},
  cache = FileNameJoin[{$UserBaseDirectory, "ApplicationData", "CertiField",
    "minorants_" <> ToString[nK] <> "_" <> ToString[nD] <> "_v" <> $CertiFieldVersion <> ".mx"}];
  If[FileExistsQ[cache],
    data = Quiet@Import[cache, "MX"];
    If[AssociationQ[data],
      tab = Catch[makeMinorantTable[data["KN"], data["DG"], data["AL"], data["BE"], data["Info"]], "CertiFieldTable"]]];
  If[! AssociationQ[tab],
    tab = Catch[buildMinorantTable[nK, nD, N[kmin]], "CertiFieldTable"];
    If[! AssociationQ[tab], Throw[failure["Internal", "Could not certify the minorant table."], "CertiField"]];
    Quiet@CreateDirectory[DirectoryName[cache]];
    Quiet@Export[cache, <|"KN" -> tab["KN"], "DG" -> tab["DG"], "AL" -> Partition[tab["ALf"], tab["nD"]],
      "BE" -> Partition[tab["BEf"], tab["nD"]], "Info" -> tab["Info"]|>, "MX"]];
  tab];

skyrmeMu[params_] := Module[{mu = Lookup[params, "PionMass", 1]},
  mu = If[InexactNumberQ[mu], Rationalize[mu, 10^-12], mu];
  If[! skyrmeMassQ[mu],
    Throw[failure["BadModel", "The pion mass must be a rational number with 0 \[LessEqual] \[Mu] < 2 (Battye\[Dash]Sutcliffe units). Grids are shipped for 0, 0.526, 0.65 and 1 and generated for other masses."], "CertiField"]];
  mu];

(* For a mass without shipped data the grid and the trial map come from Kernel/SkyrmeAnyMass.wl (skyrmeGridPlan,
   skyrmeTrial); both are heuristic.  The lower bound may be computed at a smaller mass (docs/SKYRME_ANY_MASS.md). *)
certifySkyrme[params_Association, cand_, opts_List] := Module[
  {muq = skyrmeMu[params], preset = OptionValue[CertifiedVariationalBound, opts, "Resolution"], verbose = TrueQ[OptionValue[CertifiedVariationalBound, opts, "Verbose"]],
   plan, grid, tab, lbRes, trial, ubRes, t0 = AbsoluteTime[], monitor},
  If[! KeyExistsQ[$skyrmePresets, preset], Throw[failure["BadOption", "\"Resolution\" must be one of " <> ToString[Keys[$skyrmePresets]] <> "."], "CertiField"]];
  plan = skyrmeGridPlan[muq, preset, OptionValue[CertifiedVariationalBound, opts, Method], verbose];
  grid = plan["Grid"];
  tab = minorantTable[$skyrmePresets[preset, "Table"]];
  monitor = If[verbose, Function[s, If[Mod[s["Level"], 250] == 0, Print["  level ", s["Level"], "/", s["Levels"], "   pairs ", s["Pairs"], "   ", Round[s["Seconds"], 0.1], " s"]]], None];
  lbRes = levelSetLowerBound[plan["Mass"], grid, tab, "Monitor" -> monitor];
  trial = Which[
    cand === Automatic, skyrmeTrial[muq, verbose],
    AssociationQ[cand] && KeyExistsQ[cand, "r"], cand,
    True, Module[{f = cand, df}, df = Derivative[1][f]; hedgehogTrialData[f, df, muq]]];
  ubRes = hedgehogUpperBound[trial, muq];
  Join[<|"LowerBound" -> lbRes["LowerBound"], "UpperBound" -> ubRes["UpperBound"], "Seconds" -> AbsoluteTime[] - t0,
    "LevelSet" -> lbRes, "Hedgehog" -> KeyDrop[ubRes, "Enclosure"], "Trial" -> trial, "Preset" -> preset, "Mu" -> muq,
    "GridMeta" -> grid["Meta"], "TableInfo" -> tab["Info"]|>, plan["Extra"]]];

(* ---- assembling the CertifiedBound ---- *)
(* Method and "RefineCandidate" are used by ScalarFieldND only (Kernel/ScalarFieldND.wl) *)
Options[CertifiedVariationalBound] = {"Resolution" -> "Quick", "Nodes" -> Automatic, "GapTolerance" -> None,
  "Verbose" -> False, WorkingPrecision -> 30, Method -> Automatic, "RefineCandidate" -> False};

(* convenience: CertifiedVariationalBound["Name", candidate?, modelParam -> v, ..., option -> v, ...] *)
$cvbOptionNames := ToString /@ Keys[Options[CertifiedVariationalBound]];
splitRules[rules_List] := {Select[rules, ! MemberQ[$cvbOptionNames, ToString[First[#]]] &], Select[rules, MemberQ[$cvbOptionNames, ToString[First[#]]] &]};
CertifiedVariationalBound[name_String, cand_ /; ! OptionQ[cand], rules___?OptionQ] :=
  With[{sp = splitRules[Flatten[{rules}]]}, CertifiedVariationalBound[CertiFieldModel[name, Association[sp[[1]]]], cand, Sequence @@ sp[[2]]]];
CertifiedVariationalBound[name_String, rules___?OptionQ] :=
  With[{sp = splitRules[Flatten[{rules}]]}, CertifiedVariationalBound[CertiFieldModel[name, Association[sp[[1]]]], Automatic, Sequence @@ sp[[2]]]];
CertifiedVariationalBound[m : CertiFieldModel[_String, _Association], opts : OptionsPattern[]] := CertifiedVariationalBound[m, Automatic, opts];
CertifiedVariationalBound[CertiFieldModel[name_String, params_Association], cand_ /; ! OptionQ[cand], opts : OptionsPattern[]] :=
  Catch[assembleCertificate[name, params, cand, {opts}], "CertiField"];

(* The reported bounds are machine numbers rounded outward from the computed enclosures (lower bound down, upper bound
   up), so they are safe to use in further floating-point work and are reproduced exactly by VerifyCertificate.
   The raw enclosure ends are kept in "Certificate". *)
assembleCertificate[name_, params_, cand_, opts_] := Module[{res, raw, base, status, tolerance, rel},
  If[! MemberQ[$CertiFieldModels, name], Message[CertifiedVariationalBound::nomodel, name, $CertiFieldModels];
    Throw[failure["UnknownModel", "Unknown model: " <> name], "CertiField"]];
  tolerance = OptionValue[CertifiedVariationalBound, opts, "GapTolerance"];
  raw = Which[
    name === "SkyrmeDegreeOne", certifySkyrme[params, cand, opts],
    name === "ScalarFieldND", certifyScalarFieldND[params, cand, opts],
    name === "LinearElasticity", certifyLinearElasticity[params, cand, opts],
    True, certifyScalarField1D[presetSpec[name, params], If[cand === Automatic && name =!= "ScalarField1D", presetCandidate[name, params], cand], opts]];
  If[! (NumericQ[raw["LowerBound"]] && NumericQ[raw["UpperBound"]]),
    Throw[failure["Inconclusive", "A bound could not be enclosed."], "CertiField"]];
  res = Join[raw, <|"LowerBound" -> fLo[exactNumber[raw["LowerBound"]]], "UpperBound" -> fHi[exactNumber[raw["UpperBound"]]],
    "RawBounds" -> {raw["LowerBound"], raw["UpperBound"]}|>];
  If[res["LowerBound"] > res["UpperBound"], Throw[failure["Internal", "Lower bound exceeds upper bound."], "CertiField"]];
  base = Which[
    name === "SkyrmeDegreeOne", skyrmeDescription[res],
    name === "ScalarFieldND", scalarNDDescription[params, res],
    name === "LinearElasticity", elasticityDescription[params, res],
    True, scalarDescription[name, params, res]];
  rel = N[relGapOf[res]];
  status = If[tolerance === None || rel <= tolerance, "Certified", "Inconclusive"];
  CertifiedBound[Join[base, <|"Status" -> status, "LowerBound" -> res["LowerBound"], "UpperBound" -> res["UpperBound"],
    "Seconds" -> res["Seconds"], "Version" -> $CertiFieldVersion, "GapTolerance" -> tolerance,
    "Certificate" -> KeyDrop[res, {"LowerBound", "UpperBound", "Seconds", "Spec"}], "Model" -> name, "Parameters" -> params|>]]];

(* exact numbers in theorem statements, as unambiguous linear text: 1/2, 2*Pi, -Infinity *)
exactString[x_] := ToString[x, InputForm];

skyrmeDescription[res_] := <|
  "Units" -> "12\[Pi]^2 (Faddeev\[Dash]Bogomolny units), Battye\[Dash]Sutcliffe normalisation",
  "Theorem" -> Row[{"Every finite-energy degree-one map U: \[DoubleStruckCapitalR]^3 \[RightArrow] S^3 with pion mass \[Mu] = ", exactString[res["Mu"]],
    " has E(U)/(12\[Pi]^2) \[GreaterEqual] ", lowerString[res["LowerBound"], 8], ".  The hedgehog trial map has E/(12\[Pi]^2) \[LessEqual] ",
    upperString[res["UpperBound"], 8], ".",
    If[Lookup[res, "LowerBoundMass", res["Mu"]] =!= res["Mu"],
      "  (The lower bound is certified at the pion mass " <> exactString[res["LowerBoundMass"]] <> " and holds here because the energy is non-decreasing in the pion mass.)", ""]}],
  "AdmissibleClass" -> "all finite-energy maps \[DoubleStruckCapitalR]^3 \[RightArrow] S^3 of degree one (no symmetry assumed)",
  "LowerBoundMethod" -> Row[{"level-set bound (coarea, isoperimetry, convexified level energy), certified dynamic programme: ",
    res["LevelSet", "Levels"], " levels, ", res["LevelSet", "Pairs"], " transitions, outward-rounded arithmetic"}],
  "UpperBoundMethod" -> Row[{"hedgehog trial map with ", res["Hedgehog", "Cells"], " cubic cells, rigorous quadrature"}],
  "TrustedBase" -> "IEEE-754 machine arithmetic (+, \[Minus], \[Times], \[Divide]), Wolfram Language Interval arithmetic, and the analytic reduction (note_levelset.tex, Theorems 1, 3, 8)",
  "NotClaimed" -> "that the hedgehog is the minimiser; the gap is a certified bound on how far below the hedgehog any degree-one map can be",
  "References" -> "Level-set lower bound: note_levelset.tex (2026), Theorem 3; computer-assisted part re-implemented in Wolfram Language."|>;

scalarDescription[name_, params_, res_] := Module[{sp = res["Spec"], dom, bv, d = Max[displayDigits[res], 10]},
  dom = sp["Domain"]; bv = sp["BoundaryValues"];
  <|"Units" -> "energy units of the model",
    "Theorem" -> Row[{"Every admissible u with u(", exactString[dom[[1]]], ") = ", exactString[bv[[1]]], ", u(", exactString[dom[[2]]], ") = ", exactString[bv[[2]]],
      " has E[u] \[GreaterEqual] ", lowerString[res["LowerBound"], d], "; the competitor has E \[LessEqual] ", upperString[res["UpperBound"], d], "."}],
    "AdmissibleClass" -> Row[{"absolutely continuous u on (", exactString[dom[[1]]], ", ", exactString[dom[[2]]], ") with u(", exactString[dom[[1]]], ") = ",
      exactString[bv[[1]]], ", u(", exactString[dom[[2]]], ") = ", exactString[bv[[2]]], " and finite energy"}],
    "LowerBoundMethod" -> Row[{"calibration k u'^2 + V(u) + c \[GreaterEqual] 2\[Sqrt](k(V+c)) |u'| with c = ", res["Calibration"], "; rigorous quadrature"}],
    "UpperBoundMethod" -> Row[{"cubic Hermite interpolant of the candidate on ", Length[res["Competitor", "x"]] - 1, " cells",
      If[TrueQ[res["Exact"]], ", energy in exact rational arithmetic", ", rigorous quadrature"]}],
    "TrustedBase" -> "Wolfram Language exact rational arithmetic, Interval arithmetic, and quantifier elimination (Resolve) for V + c \[GreaterEqual] 0",
    "References" -> "Bogomolny (1976); Modica\[Dash]Mortola (1977); the calibration family with c is standard (energy integral of the Euler\[Dash]Lagrange equation)."|>];

(* ---- independent re-verification from the stored certificate ----
   Recomputes every rigorous step from the recorded inputs (grid preset, minorant-table parameters, trial map,
   calibration constant, competitor data); nothing numerical is copied from the stored result. *)
VerifyCertificate[CertifiedBound[a_Association]] := TrueQ[Catch[verifyCertificate[a], "CertiField"]];

verifyCertificate[a_Association] := Module[{c = a["Certificate"], name = a["Model"], lb, ub, spec, prec = 30},
  If[name === "ScalarFieldND", Return[verifyScalarFieldND[a]]];
  If[name === "LinearElasticity", Return[verifyLinearElasticity[a]]];
  Which[
    name === "SkyrmeDegreeOne",
      (* the mass at which the lower bound was computed (at most the mass of the model) and its grid, shipped or stored *)
      With[{d = skyrmeCertificateData[c]},
        lb = levelSetLowerBound[d[[1]], d[[2]], minorantTable[$skyrmePresets[c["Preset"], "Table"]]]["LowerBound"]];
      ub = hedgehogUpperBound[c["Trial"], c["Mu"]]["UpperBound"],
    True,
      spec = scalarField1DSpec[presetSpec[name, a["Parameters"]]];
      If[! provedNonnegative[spec, c["Calibration"]], Return[False]];
      lb = Min[scalarLB[spec, c["Calibration"], prec]];
      ub = Max[competitorEnergy[spec, c["Competitor"], prec]]];
  NumericQ[lb] && NumericQ[ub] && lb >= a["LowerBound"] && ub <= a["UpperBound"]];

End[];
