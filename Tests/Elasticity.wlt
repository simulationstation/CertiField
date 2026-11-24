(* ::Package:: *)

(* LinearElasticity suite (under a minute).  Run with:  wolframscript -file tools/run_tests.wls --elasticity
   The expected numbers of the test vectors come from the Python reference (Reference/Elasticity).  Three kinds of
   inequality must never be relaxed: "lower <= exact minimum <= upper", "certified lower bound <= exact value of the
   bound" and "lower bound of any admissible stress field <= energy of any admissible displacement".  Tolerances on
   agreement with the reference are stated with their reason. *)

Needs["CertiField`"];

elRef = FileNameJoin[{CertiField`Private`$CertiFieldDirectory, "Reference", "Elasticity"}];
elTV = FileNameJoin[{elRef, "testvectors"}];
elIndex = Import[FileNameJoin[{elTV, "index.json"}], "RawJSON"];
elNames = {"dirichlet_p1", "bodyforce_p2", "block_p1", "threefree", "random", "corner"};
elClose[a_, b_, tol_] := Abs[a - b] <= tol (1 + Abs[b]);
elQ[x_] := SetPrecision[x, Infinity];

(* a test vector: the model built in the Wolfram Language, with the reference's boundary data, candidate and stress field *)
elCase[name_String] := elCase[name] = Module[{dir = FileNameJoin[{elTV, name}], spec, g},
  spec = CertiField`Private`linearElasticitySpec[ToExpression[elIndex[name, "wolfram_model"]]];
  g = CertiField`Private`elImportArray[spec, dir, "gval.f64"];
  <|"SameData" -> Max[Abs[g - spec["G"]]] <= 10.^-14,
    "Spec" -> Append[spec, "G" -> g],
    "Field" -> CertiField`Private`elPackField[spec, CertiField`Private`elImportField[spec, dir]],
    "U" -> CertiField`Private`elImportArray[spec, dir, "U.f64"],
    "Expected" -> elIndex[name, "expected"]|>];

(* ---------- an independent exact evaluation ----------
   Both fields are rebuilt cell by cell as polynomials with rational coefficients and everything is done with
   symbolic polynomial algebra (Expand, D, Integrate).  No code is shared with the evaluators of Kernel/Elasticity.wl
   except elTraces, whose result (the nodal values on the traction-free faces) is checked here: the traction of the
   rebuilt field must vanish identically on those faces. *)
elHermPoly[xi_] := {1 - 3 xi^2 + 2 xi^3, xi - 2 xi^2 + xi^3, 3 xi^2 - 2 xi^3, -xi^2 + xi^3};
elLagPoly[1, xi_] := {1 - xi, xi};
elLagPoly[2, xi_] := {(1 - xi) (1 - 2 xi), 4 xi (1 - xi), xi (2 xi - 1)};

elExactFields[spec_, u_, f_] := Module[{n1, n2, h = spec["h"], p = spec["p"], vn = spec["VNodes"], b = spec["b"], tr, xx, par, uu},
  {n1, n2} = spec["Cells"];
  tr = CertiField`Private`elTraces[spec];
  xx = Map[elQ, Normal[f["X"]], {3}];
  If[tr["Nodes"] === {},
    Do[xx[[t, 1, 1]] = 0, {t, 3}],
    par = elQ /@ Normal[f["Parameters"]];
    Do[xx[[t, tr["Nodes"][[q, 1]], tr["Nodes"][[q, 2]]]] = tr["Tb"][[q, t]] + Sum[par[[k]] tr["T"][[k, q, t]], {k, Length[par]}], {q, Length[tr["Nodes"]]}, {t, 4}]];
  uu = Map[elQ, Normal[u], {3}];
  Table[
    Module[{x0 = vn[[1, c1]], y0 = vn[[2, c2]], hx, hy, lx, ly, phi},
      hx = elHermPoly[(\[FormalX] - x0)/h[[1]]]; hy = elHermPoly[(\[FormalY] - y0)/h[[2]]];
      lx = elLagPoly[p, (\[FormalX] - x0)/h[[1]]]; ly = elLagPoly[p, (\[FormalY] - y0)/h[[2]]];
      phi = Sum[xx[[1 + d1 + 2 d2, c1 + e1, c2 + e2]] h[[1]]^d1 h[[2]]^d2 hx[[2 e1 + d1 + 1]] hy[[2 e2 + d2 + 1]], {e1, 0, 1}, {d1, 0, 1}, {e2, 0, 1}, {d2, 0, 1}];
      <|"s" -> Expand[{D[phi, {\[FormalY], 2}] - b[[1]] \[FormalX], D[phi, {\[FormalX], 2}] - b[[2]] \[FormalY], -D[phi, \[FormalX], \[FormalY]]}],
        "u" -> Expand[Table[Sum[uu[[a, p (c1 - 1) + i + 1, p (c2 - 1) + j + 1]] lx[[i + 1]] ly[[j + 1]], {i, 0, p}, {j, 0, p}], {a, 2}]],
        "box" -> {x0, x0 + h[[1]], y0, y0 + h[[2]]}|>],
    {c1, n1}, {c2, n2}]];

(* <|"Admissible" -> the hypotheses of Theorem 1 hold exactly and the candidate is continuous, "LB" -> the exact value
   of the bound, "E" -> the exact energy of the candidate, "Identity" -> E - LB = 1/2 Integrate[(e - S s).D.(e - S s)]|> *)
elExactCheck[spec_, u_, f_] := Module[
  {cells = elExactFields[spec, u, f], n1, n2, b = spec["b"], dm = spec["D"], sm = spec["S"], dir = spec["Dirichlet"], x = \[FormalX], y = \[FormalY],
   trac, strain, int2, zero, equil, cont, free, en, comp, defect, bnd},
  {n1, n2} = spec["Cells"];
  trac[s_, i_] := If[i == 1, {s[[1]], s[[3]]}, {s[[3]], s[[2]]}];
  strain[v_] := {D[v[[1]], x], D[v[[2]], y], D[v[[1]], y] + D[v[[2]], x]};
  int2[e_, {x0_, x1_, y0_, y1_}] := Integrate[Expand[e], {x, x0, x1}, {y, y0, y1}];
  zero[e_] := AllTrue[Flatten[{Expand[e]}], # === 0 &];
  equil = And @@ Flatten[Table[With[{s = cells[[c1, c2]]["s"]},
      zero[D[s[[1]], x] + D[s[[3]], y] + b[[1]]] && zero[D[s[[3]], x] + D[s[[2]], y] + b[[2]]]], {c1, n1}, {c2, n2}]];
  cont = And @@ Flatten[{
      Table[With[{a = cells[[c1, c2]], c = cells[[c1 + 1, c2]], xf = cells[[c1, c2]]["box"][[2]]},
        zero[(trac[a["s"], 1] - trac[c["s"], 1]) /. x -> xf] && zero[(a["u"] - c["u"]) /. x -> xf]], {c1, n1 - 1}, {c2, n2}],
      Table[With[{a = cells[[c1, c2]], c = cells[[c1, c2 + 1]], yf = cells[[c1, c2]]["box"][[4]]},
        zero[(trac[a["s"], 2] - trac[c["s"], 2]) /. y -> yf] && zero[(a["u"] - c["u"]) /. y -> yf]], {c1, n1}, {c2, n2 - 1}]}];
  free = And @@ Flatten[{
      If[dir[[1, 1]], {}, Table[zero[trac[cells[[1, c2]]["s"], 1] /. x -> cells[[1, c2]]["box"][[1]]], {c2, n2}]],
      If[dir[[1, 2]], {}, Table[zero[trac[cells[[n1, c2]]["s"], 1] /. x -> cells[[n1, c2]]["box"][[2]]], {c2, n2}]],
      If[dir[[2, 1]], {}, Table[zero[trac[cells[[c1, 1]]["s"], 2] /. y -> cells[[c1, 1]]["box"][[3]]], {c1, n1}]],
      If[dir[[2, 2]], {}, Table[zero[trac[cells[[c1, n2]]["s"], 2] /. y -> cells[[c1, n2]]["box"][[4]]], {c1, n1}]]}];
  en = Sum[With[{c = cells[[c1, c2]]}, With[{e = strain[c["u"]]}, int2[e . dm . e/2 - b . c["u"], c["box"]]]], {c1, n1}, {c2, n2}];
  comp = Sum[With[{c = cells[[c1, c2]]}, int2[c["s"] . sm . c["s"]/2, c["box"]]], {c1, n1}, {c2, n2}];
  defect = Sum[With[{c = cells[[c1, c2]]}, With[{r = strain[c["u"]] - sm . c["s"]}, int2[r . dm . r/2, c["box"]]]], {c1, n1}, {c2, n2}];
  (* on the Dirichlet faces the trace of the candidate is the datum *)
  bnd = Total[Flatten[{
      If[dir[[1, 1]], Table[With[{c = cells[[1, c2]]}, -Integrate[Expand[(trac[c["s"], 1] . c["u"]) /. x -> c["box"][[1]]], {y, c["box"][[3]], c["box"][[4]]}]], {c2, n2}], {}],
      If[dir[[1, 2]], Table[With[{c = cells[[n1, c2]]}, Integrate[Expand[(trac[c["s"], 1] . c["u"]) /. x -> c["box"][[2]]], {y, c["box"][[3]], c["box"][[4]]}]], {c2, n2}], {}],
      If[dir[[2, 1]], Table[With[{c = cells[[c1, 1]]}, -Integrate[Expand[(trac[c["s"], 2] . c["u"]) /. y -> c["box"][[3]]], {x, c["box"][[1]], c["box"][[2]]}]], {c1, n1}], {}],
      If[dir[[2, 2]], Table[With[{c = cells[[c1, n2]]}, Integrate[Expand[(trac[c["s"], 2] . c["u"]) /. y -> c["box"][[4]]], {x, c["box"][[1]], c["box"][[2]]}]], {c1, n1}], {}]}]];
  <|"Admissible" -> equil && cont && free, "LB" -> bnd - comp, "E" -> en, "Identity" -> (en - (bnd - comp) == defect && defect >= 0)|>];

(* the exact evaluation against the two certified evaluators, for given fields *)
elExactAgainstCertified[spec_, u_, f_] := Module[{ex = elExactCheck[spec, u, f], lb, ub},
  lb = CertiField`Private`elLowerBound[spec, f]["Lower"];
  ub = CertiField`Private`elUpperBound[spec, u];
  {ex["Admissible"], ex["Identity"], lb <= ex["LB"], ex["LB"] - lb <= 10^-9 (1 + Abs[ex["LB"]]), ub[[1]] <= ex["E"] <= ub[[2]], ub[[2]] - ub[[1]] <= 10^-9 (1 + Abs[ex["E"]])}];

(* random admissible fields: any nodal values and parameters give a statically admissible stress field, any nodal
   values with the Dirichlet values imposed give an admissible displacement *)
elRandomField[spec_, f_, amp_] := <|"X" -> f["X"] + amp RandomReal[{-1, 1}, Dimensions[f["X"]]],
  "Parameters" -> If[f["Parameters"] === {}, {}, f["Parameters"] + amp RandomReal[{-1, 1}, Length[f["Parameters"]]]]|>;
elRandomDisplacement[spec_, u_, amp_] := CertiField`Private`elImposeDirichlet[spec, u + amp RandomReal[{-1, 1}, Dimensions[u]]];

(* ---------- primitives ---------- *)

(* the exact one-dimensional integrals against symbolic integration *)
VerificationTest[
  Module[{xi, herm, l1, l2, hc = CertiField`Private`$elHermite},
    herm = elHermPoly[xi]; l1 = elLagPoly[1, xi]; l2 = elLagPoly[2, xi];
    {Table[CertiField`Private`elIMat[hc, a, hc, c] === Table[Integrate[D[herm[[i]], {xi, a}] D[herm[[j]], {xi, c}], {xi, 0, 1}], {i, 4}, {j, 4}], {a, 0, 2}, {c, 0, 2}],
     Table[CertiField`Private`elIMat[hc, a, CertiField`Private`elLagrange[2], 0, k] === Table[Integrate[xi^k D[herm[[i]], {xi, a}] l2[[j]], {xi, 0, 1}], {i, 4}, {j, 3}], {a, 0, 2}, {k, 0, 1}],
     Table[CertiField`Private`elIMat[CertiField`Private`elLagrange[1], a, CertiField`Private`elLagrange[1], c] === Table[Integrate[D[l1[[i]], {xi, a}] D[l1[[j]], {xi, c}], {xi, 0, 1}], {i, 2}, {j, 2}], {a, 0, 1}, {c, 0, 1}],
     Table[CertiField`Private`elIMat[CertiField`Private`elLagrange[2], a, CertiField`Private`elLagrange[2], c] === Table[Integrate[D[l2[[i]], {xi, a}] D[l2[[j]], {xi, c}], {xi, 0, 1}], {i, 3}, {j, 3}], {a, 0, 1}, {c, 0, 1}]}],
  {ConstantArray[True, {3, 3}], ConstantArray[True, {3, 2}], ConstantArray[True, {2, 2}], ConstantArray[True, {2, 2}]},
  TestID -> "el-exact-one-dimensional-integrals"]

(* the exact element matrices are the cell integrals they stand for: random rational nodal values on one cell *)
VerificationTest[
  Module[{x, y, h = {2/3, 5/7}, x0 = {1/3, -2}, dm = {{4, 1, 1/2}, {1, 3, -1/4}, {1/2, -1/4, 2}}, sm, b = {1/5, -2/3}, loc, vals, phi, s, sb, res = {}},
    sm = Inverse[dm];
    SeedRandom[11];
    (* the stress space *)
    loc = CertiField`Private`elAiryLocal[h, sm, b];
    vals = RandomInteger[{-9, 9}, 16]/RandomInteger[{1, 9}, 16];
    phi = Sum[With[{k = loc["Dofs"][[j]]}, vals[[j]] h[[1]]^Mod[k[[1]], 2] h[[2]]^Mod[k[[2]], 2] elHermPoly[(x - x0[[1]])/h[[1]]][[k[[1]] + 1]] elHermPoly[(y - x0[[2]])/h[[2]]][[k[[2]] + 1]]], {j, 16}];
    s = {D[phi, {y, 2}], D[phi, {x, 2}], -D[phi, x, y]}; sb = {-b[[1]] x, -b[[2]] y, 0};
    AppendTo[res, Integrate[Expand[s . sm . s], {x, x0[[1]], x0[[1]] + h[[1]]}, {y, x0[[2]], x0[[2]] + h[[2]]}] == vals . loc["A"] . vals];
    AppendTo[res, Integrate[Expand[s . sm . sb], {x, x0[[1]], x0[[1]] + h[[1]]}, {y, x0[[2]], x0[[2]] + h[[2]]}] == vals . (loc["w"] . {1, x0[[1]], x0[[2]]})];
    (* the displacement spaces *)
    Do[
      Module[{st = CertiField`Private`elStiffnessLocal[p, h, dm, b], uv, u, e},
        uv = RandomInteger[{-9, 9}, Length[st["Dofs"]]]/RandomInteger[{1, 9}, Length[st["Dofs"]]];
        u = Table[Sum[With[{d = st["Dofs"][[j]]}, If[d[[1]] == a, uv[[j]] elLagPoly[p, (x - x0[[1]])/h[[1]]][[d[[2]] + 1]] elLagPoly[p, (y - x0[[2]])/h[[2]]][[d[[3]] + 1]], 0]], {j, Length[uv]}], {a, 2}];
        e = {D[u[[1]], x], D[u[[2]], y], D[u[[1]], y] + D[u[[2]], x]};
        AppendTo[res, Integrate[Expand[e . dm . e], {x, x0[[1]], x0[[1]] + h[[1]]}, {y, x0[[2]], x0[[2]] + h[[2]]}] == uv . st["K"] . uv];
        AppendTo[res, Integrate[Expand[b . u], {x, x0[[1]], x0[[1]] + h[[1]]}, {y, x0[[2]], x0[[2]] + h[[2]]}] == uv . st["f"]]],
      {p, 2}];
    res],
  {True, True, True, True, True, True}, TestID -> "el-exact-element-matrices"]

(* the families of polynomials on the traction-free faces: number of nodes and of parameters *)
VerificationTest[
  Module[{sp},
    sp[free_, b_] := CertiField`Private`linearElasticitySpec[<|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {1/3, 2}}, "Cells" -> {3, 4}, "TractionFreeFaces" -> free, "BodyForce" -> b|>];
    Table[With[{tr = CertiField`Private`elTraces[sp[fr, {1/2, -1/3}]]}, {Length[tr["Nodes"]], Length[tr["T"]]}],
      {fr, {{}, {{1, "Lower"}}, {{2, "Upper"}}, {{1, "Lower"}, {2, "Lower"}}, {{1, "Lower"}, {1, "Upper"}}, {{2, "Lower"}, {2, "Upper"}}, {{1, "Upper"}, {2, "Lower"}, {2, "Upper"}}}}]],
  {{0, 0}, {5, 0}, {4, 0}, {8, 0}, {10, 3}, {8, 3}, {11, 0}}, TestID -> "el-traction-free-families"]

(* ---------- the rigorous evaluators against the reference, on identical inputs ----------
   Tolerance 10^-9: the reference widens every operation by one unit in the last place, the paclet by 2^-49 relative,
   and the port sums the cells exactly; the observed differences are below 10^-11. *)

Do[
  With[{name = nm},
    VerificationTest[elCase[name]["SameData"], True, TestID -> "el-" <> name <> "-boundary-data-agree"];
    VerificationTest[
      With[{tr = CertiField`Private`elTraces[elCase[name]["Spec"]], ex = elCase[name]["Expected"]}, {Length[tr["T"]], Length[tr["Nodes"]]} === {ex["parameters"], ex["traction_free_nodes"]}],
      True, TestID -> "el-" <> name <> "-parameters"];
    VerificationTest[
      Module[{c = elCase[name], ub},
        ub = CertiField`Private`elUpperBound[c["Spec"], c["U"]];
        {ub[[1]] <= ub[[2]], elClose[N[ub[[1]]], c["Expected", "upper_enclosure"][[1]], 10^-9], elClose[N[ub[[2]]], c["Expected", "upper_enclosure"][[2]], 10^-9]}],
      {True, True, True}, TestID -> "el-" <> name <> "-competitor-energy"];
    VerificationTest[
      Module[{c = elCase[name], lb},
        lb = CertiField`Private`elLowerBound[c["Spec"], c["Field"]];
        {elClose[N[lb["Lower"]], c["Expected", "lower"], 10^-9],
         elClose[N[lb["Boundary"][[1]]], c["Expected", "boundary"][[1]], 10^-9], elClose[N[lb["Boundary"][[2]]], c["Expected", "boundary"][[2]], 10^-9],
         elClose[N[lb["IntQ"][[1]]], c["Expected", "intQ"][[1]], 10^-9], elClose[N[lb["IntQ"][[2]]], c["Expected", "intQ"][[2]], 10^-9]}],
      {True, True, True, True, True}, TestID -> "el-" <> name <> "-lower-bound"];
    (* the same fields in exact arithmetic: static admissibility, both exact values inside the enclosures *)
    VerificationTest[
      With[{c = elCase[name]}, elExactAgainstCertified[c["Spec"], c["U"], c["Field"]]],
      {True, True, True, True, True, True}, TestID -> "el-" <> name <> "-exact-evaluation"];
    (* and for fields that no solver produced *)
    VerificationTest[
      Module[{c = elCase[name]},
        SeedRandom[Hash[name]];
        elExactAgainstCertified[c["Spec"], elRandomDisplacement[c["Spec"], c["U"], 0.3], elRandomField[c["Spec"], c["Field"], 0.3]]],
      {True, True, True, True, True, True}, TestID -> "el-" <> name <> "-exact-evaluation-random-fields"]],
  {nm, elNames}]

(* the stored certificates of the reference are recomputed (both bounds at least as strong as the stored claims) *)
VerificationTest[
  Table[CertiField`Private`elVerifyReferenceCertificate[FileNameJoin[{elRef, "results", nm}]]["Verified"], {nm, {"bending", "block_p1_16", "block_p1_32", "block_p2_16", "block_p2_32", "gravity"}}],
  {True, True, True, True, True, True}, TestID -> "el-reference-certificates-verify"]

(* the stress field found by the Wolfram Language solve is as good as the reference's, and so is the candidate *)
VerificationTest[
  Table[
    Module[{c = elCase[nm], spec, lb, ub},
      spec = c["Spec"];
      lb = CertiField`Private`elLowerBound[spec, CertiField`Private`elSolveStress[spec]]["Lower"];
      ub = CertiField`Private`elUpperBound[spec, CertiField`Private`elSolveCandidate[spec]][[2]];
      If[nm === "random", lb > c["Expected", "lower"] && ub < c["Expected", "upper_enclosure"][[2]],
        elClose[N[lb], c["Expected", "lower"], 10^-8] && elClose[N[ub], c["Expected", "upper_enclosure"][[2]], 10^-8]]],
    {nm, elNames}],
  {True, True, True, True, True, True}, TestID -> "el-solves-match-reference"]

(* ---------- E-T1, E-T2: minima known exactly ----------
   Data, grids and compliance are dyadic, so the stored nodal values are exactly the data.  The exact fields lie in
   the discrete spaces, so both bounds must equal the minimum to rounding. *)

elKnown = <|
  (* homogeneous strain: u = A x + c on all faces; minimum |Omega| e0.D.e0/2 = 101/256 *)
  "affine-p1" -> {<|"Lame" -> {1, 1}, "Domain" -> {{0, 2}, {0, 1}}, "Cells" -> {4, 2}, "Degree" -> 1,
      "BoundaryData" -> Function[{x, y}, {0.125 x - 0.1875 y + 0.5, 0.5 x + 0.25 y - 0.25}]|>, 101/256},
  "affine-p2" -> {<|"Lame" -> {1, 1}, "Domain" -> {{0, 2}, {0, 1}}, "Cells" -> {4, 2}, "Degree" -> 2,
      "BoundaryData" -> Function[{x, y}, {0.125 x - 0.1875 y + 0.5, 0.5 x + 0.25 y - 0.25}]|>, 101/256},
  (* uniaxial stress (1/2, 0, 0) with free sides: u = (3 x/16, -y/16); minimum |Omega| s^2 S11/2 = 3/32 *)
  "tension-p1" -> {<|"Lame" -> {1, 1}, "Domain" -> {{0, 2}, {0, 1}}, "Cells" -> {4, 2}, "Degree" -> 1, "TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}},
      "BoundaryData" -> Function[{x, y}, {0.1875 x, -0.0625 y}]|>, 3/32},
  "tension-p2" -> {<|"Lame" -> {1, 1}, "Domain" -> {{0, 2}, {0, 1}}, "Cells" -> {4, 2}, "Degree" -> 2, "TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}},
      "BoundaryData" -> Function[{x, y}, {0.1875 x, -0.0625 y}]|>, 3/32},
  (* a quadratic displacement with the constant body force that equilibrates it; the minimum is computed below *)
  "quadratic" -> {<|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> {2, 2}, "Degree" -> 2, "BodyForce" -> {-1, -1},
      "BoundaryData" -> Function[{x, y}, {x^2/4 + x y/8, y^2/8 - x y/4 + x/2}]|>, Automatic},
  (* pure bending, stress (y, 0, 0), free sides: u = (3 x y/8, -y^2/16 - 3 x^2/16); minimum S11/12 = 1/32 *)
  "bending" -> {<|"Lame" -> {1, 1}, "Domain" -> {{0, 2}, {-1/2, 1/2}}, "Cells" -> {4, 2}, "Degree" -> 2, "TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}},
      "BoundaryData" -> Function[{x, y}, {0.375 x y, -0.0625 y y - 0.1875 x x}]|>, 1/32}|>;

(* the quadratic displacement: its body force and its energy, by symbolic integration *)
elQuadraticCase = Module[{x, y, dm = {{3, 1, 0}, {1, 3, 0}, {0, 0, 1}}, u, e, s, b},
  u = {x^2/4 + x y/8, y^2/8 - x y/4 + x/2};
  e = {D[u[[1]], x], D[u[[2]], y], D[u[[1]], y] + D[u[[2]], x]}; s = dm . e;
  b = -{D[s[[1]], x] + D[s[[3]], y], D[s[[3]], x] + D[s[[2]], y]};
  <|"BodyForce" -> b, "Minimum" -> Integrate[e . s/2 - b . u, {x, 0, 1}, {y, 0, 1}]|>];

VerificationTest[elQuadraticCase, <|"BodyForce" -> {-1, -1}, "Minimum" -> 121/192|>, TestID -> "el-quadratic-case-data"]

elKnownMinimum[name_] := Replace[elKnown[name][[2]], Automatic :> elQuadraticCase["Minimum"]];

Do[
  With[{name = nm},
    (* the two evaluators: the exact minimum is enclosed, and reached to rounding *)
    VerificationTest[
      Module[{spec = CertiField`Private`linearElasticitySpec[elKnown[name][[1]]], ex = elKnownMinimum[name], lb, ub},
        lb = CertiField`Private`elLowerBound[spec, CertiField`Private`elSolveStress[spec]]["Lower"];
        ub = CertiField`Private`elUpperBound[spec, CertiField`Private`elSolveCandidate[spec]][[2]];
        {lb <= ex <= ub, ex - lb <= 10^-10 (1 + Abs[ex]), ub - ex <= 10^-10 (1 + Abs[ex])}],
      {True, True, True}, TestID -> "el-known-minimum-" <> name <> "-reached-by-both-bounds"];
    (* the public interface *)
    VerificationTest[
      With[{c = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", elKnown[name][[1]]]], ex = elKnownMinimum[name]},
        {c["Status"], c["LowerBound"] <= ex <= c["UpperBound"], c["RelativeGap"] < 10^-8, VerifyCertificate[c]}],
      {"Certified", True, True, True}, TestID -> "el-known-minimum-" <> name <> "-certificate"]],
  {nm, Keys[elKnown]}]

(* ---------- E-T3: falsification ----------
   Every statically admissible stress field and every admissible displacement of the same problem must satisfy
   LB(stress) <= E[displacement].  Random fields, perturbed optimal fields and optimal fields of other grids are
   paired; one violation would expose an error in the theorem or in the evaluation. *)

elProblems = <|
  "block" -> Function[{n, p}, <|"Lame" -> {3/2, 1}, "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> n, "Degree" -> p, "TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}},
      "BoundaryData" -> Function[{x, y}, {0.1 x, 0}]|>],
  "dirichlet-body-force" -> Function[{n, p}, <|"Lame" -> {3/2, 1}, "Domain" -> {{0, 2}, {0, 1}}, "Cells" -> n, "Degree" -> p, "BodyForce" -> {1/2, -1/3},
      "BoundaryData" -> Function[{x, y}, {0.1 x - 0.15 y + 0.3 x y, 0.25 x + 0.2 y y - 0.1}]|>],
  "three-free-anisotropic" -> Function[{n, p}, <|"StiffnessMatrix" -> {{4, 1, 1/2}, {1, 3, -1/4}, {1/2, -1/4, 2}}, "Domain" -> {{1/3, 1}, {0, 1}}, "Cells" -> n, "Degree" -> p,
      "BodyForce" -> {1/5, -2/3}, "TractionFreeFaces" -> {{1, "Upper"}, {2, "Lower"}, {2, "Upper"}}, "BoundaryData" -> Function[{x, y}, {0.1 x - 0.15 y + 0.3 x y, 0.25 x + 0.2 y y - 0.1}]|>],
  "corner" -> Function[{n, p}, <|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> n, "Degree" -> p, "BodyForce" -> {1, 1},
      "TractionFreeFaces" -> {{1, "Lower"}, {2, "Lower"}}, "BoundaryData" -> Function[{x, y}, {0.1 x - 0.15 y + 0.3 x y, 0.25 x + 0.2 y y - 0.1}]|>]|>;

(* within one problem instance: the smallest E - LB over all pairs of an admissible displacement and an admissible
   stress field (optimal ones and perturbations of five sizes).  The energy of every displacement is its lower
   enclosure end here and every bound is certified, so a negative value would be a certified contradiction. *)
elSmallestMargin[params_] := Module[{spec = CertiField`Private`linearElasticitySpec[params], u, f, lbs, ens},
  u = CertiField`Private`elSolveCandidate[spec]; f = CertiField`Private`elSolveStress[spec];
  lbs = Join[{CertiField`Private`elLowerBound[spec, f]["Lower"]},
    Flatten[Table[CertiField`Private`elLowerBound[spec, elRandomField[spec, f, amp]]["Lower"], {amp, {10.^-4, 10.^-3, 10.^-2, 0.1, 1.}}, {3}]]];
  ens = Join[{CertiField`Private`elUpperBound[spec, u][[1]]},
    Flatten[Table[CertiField`Private`elUpperBound[spec, elRandomDisplacement[spec, u, amp]][[1]], {amp, {10.^-4, 10.^-3, 10.^-2, 0.1, 1.}}, {3}]]];
  Min[ens] - Max[lbs]];

Do[
  With[{name = nm},
    VerificationTest[
      (SeedRandom[Hash[name]];
       Min[Table[elSmallestMargin[elProblems[name][n, p]], {n, {{2, 2}, {5, 3}, {8, 8}}}, {p, 2}]] >= 0),
      True, TestID -> "el-falsification-" <> name]],
  {nm, Keys[elProblems]}]

(* the data of the block are affine on the Dirichlet faces, so the certified problem is the same on every grid and
   for both degrees: all lower bounds lie below all upper bounds, and both improve under refinement *)
elBlock = Association[Table[{n, p} -> CertifiedVariationalBound[CertiFieldModel["LinearElasticity", elProblems["block"][n, p]]], {n, {4, 8, 16, 32}}, {p, 2}]];

VerificationTest[AllTrue[Values[elBlock], #["Status"] === "Certified" &], True, TestID -> "el-block-status"]
VerificationTest[Max[#["LowerBound"] & /@ Values[elBlock]] <= Min[#["UpperBound"] & /@ Values[elBlock]], True, TestID -> "el-block-bounds-of-all-grids-never-cross"]
VerificationTest[
  {Less @@ Table[elBlock[{n, 1}]["LowerBound"], {n, {4, 8, 16, 32}}], Greater @@ Table[elBlock[{n, 1}]["UpperBound"], {n, {4, 8, 16, 32}}],
   Greater @@ Table[elBlock[{n, 2}]["UpperBound"], {n, {4, 8, 16, 32}}], And @@ Table[elBlock[{n, 2}]["UpperBound"] < elBlock[{n, 1}]["UpperBound"], {n, {4, 8, 16, 32}}]},
  {True, True, True, True}, TestID -> "el-block-bounds-improve-under-refinement"]

(* ---------- E-T4: the block clamped on one face, pulled on the opposite one, free elsewhere ---------- *)

VerificationTest[
  Module[{t, c},
    {t, c} = AbsoluteTiming[CertifiedVariationalBound[CertiFieldModel["LinearElasticity", elProblems["block"][16, 1]]]];
    {c["RelativeGap"] < 0.01, t < 60, c["LowerBound"] > 0.01501, c["UpperBound"] < 0.01507}],
  {True, True, True, True}, TestID -> "el-block-gap-below-one-percent-in-under-a-minute"]
VerificationTest[
  {elBlock[{32, 2}]["RelativeGap"] < 0.0005, elBlock[{32, 1}]["RelativeGap"] < 0.0013, elBlock[{16, 2}]["RelativeGap"] < 0.0012},
  {True, True, True}, TestID -> "el-block-gaps-as-in-the-reference"]

(* ---------- certificates ---------- *)

elCert = elBlock[{16, 2}];
elTamper[key_, fun_] := With[{a = elCert["Association"]}, CertifiedBound[Append[a, "Certificate" -> Append[a["Certificate"], key -> fun[a["Certificate"][key]]]]]];

VerificationTest[VerifyCertificate[elCert], True, TestID -> "el-verify-genuine-certificate"]
VerificationTest[VerifyCertificate[CertifiedBound[Append[elCert["Association"], "LowerBound" -> 0.01502]]], False, TestID -> "el-verify-rejects-tampered-lower-bound"]
VerificationTest[VerifyCertificate[CertifiedBound[Append[elCert["Association"], "UpperBound" -> 0.01503]]], False, TestID -> "el-verify-rejects-tampered-upper-bound"]
VerificationTest[VerifyCertificate[elTamper["StressField", Append[#, "X" -> 1.001 #["X"]] &]], False, TestID -> "el-verify-rejects-tampered-stress-field"]
VerificationTest[VerifyCertificate[elTamper["StressField", Append[#, "Parameters" -> 1.001 #["Parameters"]] &]], False, TestID -> "el-verify-rejects-tampered-parameters"]
VerificationTest[VerifyCertificate[elTamper["Nodal", ReplacePart[#, {1, 5, 5} -> #[[1, 5, 5]] + 0.001] &]], False, TestID -> "el-verify-rejects-tampered-candidate"]
VerificationTest[VerifyCertificate[elTamper["Nodal", ReplacePart[#, {1, 1, 5} -> 0.001] &]], False, TestID -> "el-verify-rejects-candidate-off-the-datum"]
VerificationTest[VerifyCertificate[elTamper["BoundaryNodalValues", 1.001 # &]], False, TestID -> "el-verify-rejects-tampered-boundary-values"]
VerificationTest[
  VerifyCertificate[CertifiedBound[Append[elCert["Association"], "Parameters" -> Append[elCert["Parameters"], "Lame" -> {3/2, 101/100}]]]],
  False, TestID -> "el-verify-rejects-another-material"]
(* a tampered claim that is weaker than the truth is still a true statement *)
VerificationTest[VerifyCertificate[CertifiedBound[Append[elCert["Association"], "LowerBound" -> 0.015]]], True, TestID -> "el-verify-accepts-a-weaker-claim"]

(* ---------- candidates ---------- *)

(* a rough candidate gives an honest wide gap; "RefineCandidate" replaces it by the finite-element solution *)
VerificationTest[
  Module[{m = CertiFieldModel["LinearElasticity", elProblems["block"][16, 1]], f = Function[{x, y}, {0.1 x + 0.02 Sin[Pi x] Cos[3 y], 0.03 x (1 - x)}], a, b},
    a = CertifiedVariationalBound[m, f]; b = CertifiedVariationalBound[m, f, "RefineCandidate" -> True];
    {a["Status"], a["RelativeGap"] > 0.05, elClose[a["LowerBound"], elBlock[{16, 1}]["LowerBound"], 10^-8], VerifyCertificate[a], b["UpperBound"] == elBlock[{16, 1}]["UpperBound"]}],
  {"Certified", True, True, True, True}, TestID -> "el-rough-candidate-gives-honest-wide-gap"]

(* nodal arrays as the candidate, and a gap tolerance that the rough candidate does not meet *)
VerificationTest[
  Module[{params = elProblems["block"][8, 2], spec, u, c},
    spec = CertiField`Private`linearElasticitySpec[params];
    u = Table[Outer[If[a == 1, 0.1 #1 &, 0. #1 &], spec["UAxes"][[1]], spec["UAxes"][[2]]], {a, 2}];
    c = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", params], u, "GapTolerance" -> 10^-4];
    {c["Status"], c["UpperBound"] >= elBlock[{8, 2}]["UpperBound"], N[c["UpperBound"]] - 0.0175 < 10^-6, VerifyCertificate[c]}],
  {"Inconclusive", True, True, True}, TestID -> "el-nodal-array-candidate-and-gap-tolerance"]

(* a finite-element solution from NDSolve as the candidate *)
VerificationTest[
  Module[{u1, u2, x, y, lam = 3/2, mu = 1, sol, m, c},
    sol = NDSolveValue[{
        Inactive[Div][{{-(lam + 2 mu), 0}, {0, -mu}} . Inactive[Grad][u1[x, y], {x, y}], {x, y}] + Inactive[Div][{{0, -lam}, {-mu, 0}} . Inactive[Grad][u2[x, y], {x, y}], {x, y}] == 0,
        Inactive[Div][{{0, -mu}, {-lam, 0}} . Inactive[Grad][u1[x, y], {x, y}], {x, y}] + Inactive[Div][{{-mu, 0}, {0, -(lam + 2 mu)}} . Inactive[Grad][u2[x, y], {x, y}], {x, y}] == 0,
        DirichletCondition[{u1[x, y] == 0, u2[x, y] == 0}, x == 0], DirichletCondition[{u1[x, y] == 1/10, u2[x, y] == 0}, x == 1]},
      {u1, u2}, {x, 0, 1}, {y, 0, 1}];
    m = CertiFieldModel["LinearElasticity", elProblems["block"][32, 2]];
    c = CertifiedVariationalBound[m, Function[{s, t}, {sol[[1]][s, t], sol[[2]][s, t]}]];
    {c["Status"], c["RelativeGap"] < 0.002, elClose[c["LowerBound"], elBlock[{32, 2}]["LowerBound"], 10^-8], c["UpperBound"] >= elBlock[{32, 2}]["UpperBound"]}],
  {"Certified", True, True, True}, TestID -> "el-ndsolve-candidate"]

(* a supplied stress field: the reference's, read from its directory; and a poor one (all zeros), which gives a poor
   but valid bound *)
VerificationTest[
  Module[{params = ToExpression[elIndex["block_p1", "wolfram_model"]], a, b},
    a = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", Append[params, "StressField" -> FileNameJoin[{elTV, "block_p1"}]]]];
    b = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", Append[params, "StressField" -> <|"X" -> ConstantArray[0, {4, 5, 5}], "Parameters" -> {0, 0, 0}|>]]];
    {a["Status"], elClose[a["LowerBound"], elIndex["block_p1", "expected", "lower"], 10^-8], VerifyCertificate[a], b["LowerBound"] <= 0 < b["UpperBound"], VerifyCertificate[b]}],
  {"Certified", True, True, True, True}, TestID -> "el-supplied-stress-field"]

(* only the boundary values of the boundary function are used: it may be singular inside the rectangle *)
VerificationTest[
  With[{c = CertifiedVariationalBound[CertiFieldModel["LinearElasticity",
      <|"Lame" -> {1, 1}, "Domain" -> {{-1, 1}, {-1, 1}}, "Cells" -> 8, "BoundaryData" -> Function[{x, y}, {Log[x^2 + y^2], x/(x^2 + y^2)}]|>]]},
    ! FailureQ[c] && c["Status"] === "Certified" && 0 < c["LowerBound"] <= c["UpperBound"] && VerifyCertificate[c]],
  True, TestID -> "el-boundary-function-singular-inside-the-rectangle"]

(* grids with a single cell; a constant datum (rigid translation), whose minimum is 0 *)
VerificationTest[
  Module[{a, b, c, d},
    a = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", <|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {0, 2}}, "Cells" -> 1, "BoundaryData" -> Function[{x, y}, {x y, x - y}]|>]];
    b = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", <|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {0, 2}}, "Cells" -> {1, 2}, "Degree" -> 2, "BodyForce" -> {1, -1},
        "TractionFreeFaces" -> {{1, "Upper"}, {2, "Lower"}, {2, "Upper"}}, "BoundaryData" -> Function[{x, y}, {y, 0}]|>]];
    c = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", <|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {0, 2}}, "Cells" -> {1, 1},
        "TractionFreeFaces" -> {{1, "Upper"}, {2, "Lower"}, {2, "Upper"}}, "BoundaryData" -> Function[{x, y}, {y, 0}]|>]];
    d = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", <|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {0, 2}}, "Cells" -> 3, "BoundaryData" -> {0.25, -1.5}|>]];
    {a["LowerBound"] <= a["UpperBound"] && VerifyCertificate[a], b["LowerBound"] <= b["UpperBound"] && VerifyCertificate[b],
     c["LowerBound"] <= c["UpperBound"] && VerifyCertificate[c], d["LowerBound"] <= 0 <= d["UpperBound"] && d["UpperBound"] - d["LowerBound"] < 10^-12 && VerifyCertificate[d]}],
  {True, True, True, True}, TestID -> "el-single-cell-grids-and-rigid-translation"]

(* ---------- presentation ---------- *)

VerificationTest[
  {StringQ[elCert["Summary"]], StringStartsQ[elCert["Summary"], "LinearElasticity: Certified"], StringContainsQ[ToString[elCert["Theorem"]], "piecewise polynomial of degree 2 through the nodal boundary values on the 16x16 grid"],
   StringContainsQ[ToString[elCert["Theorem"]], "{{7/2, 3/2, 0}, {3/2, 7/2, 0}, {0, 0, 1}}"], StringContainsQ[ToString[elCert["Theorem"]], "no traction on the faces x2 = 0, x2 = 1"],
   StringContainsQ[ToString[elCert["NotClaimed"]], "minimiser"], StringContainsQ[ToString[elCert["TrustedBase"]], "docs/ELASTICITY.md"]},
  {True, True, True, True, True, True, True}, TestID -> "el-theorem-states-the-certified-problem"]
VerificationTest[Head[ToBoxes[elCert]], InterpretationBox | TagBox | TemplateBox, SameTest -> MatchQ, TestID -> "el-summary-box-renders"]
VerificationTest[
  {MemberQ[$CertiFieldModels, "LinearElasticity"], CertiFieldModelInfo["LinearElasticity"]["Status"], CertiFieldModelInfo["LinearElasticity"]["Reference"],
   KeyExistsQ[CertiFieldModelInfo["LinearElasticity"]["Parameters"], #] & /@ {"Lame", "StiffnessMatrix", "BodyForce", "Domain", "Cells", "Degree", "BoundaryData", "TractionFreeFaces", "StressField"}},
  {True, "Supported", "docs/ELASTICITY.md", ConstantArray[True, 9]}, TestID -> "el-model-is-registered"]

(* ---------- failure modes ---------- *)

elFails[params_, rest___] := FailureQ[CertifiedVariationalBound[CertiFieldModel["LinearElasticity", params], rest]];
elBase = <|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 4|>;

VerificationTest[
  elFails[Append[elBase, "TractionFreeFaces" -> {{1, "Lower"}, {1, "Upper"}, {2, "Lower"}, {2, "Upper"}}]],
  True, TestID -> "el-no-dirichlet-face-refused"]
VerificationTest[
  {elFails[Append[elBase, "Lame" -> {1, 0}]], elFails[Append[elBase, "Lame" -> {-2, 1}]], elFails[Append[elBase, "Lame" -> {1, -1}]],
   elFails[Append[KeyDrop[elBase, "Lame"], "StiffnessMatrix" -> {{1, 2, 0}, {2, 1, 0}, {0, 0, 1}}]],
   elFails[Append[KeyDrop[elBase, "Lame"], "StiffnessMatrix" -> {{2, 1, 0}, {0, 2, 0}, {0, 0, 1}}]],
   elFails[Append[elBase, "StiffnessMatrix" -> {{2, 0, 0}, {0, 2, 0}, {0, 0, 1}}]], elFails[KeyDrop[elBase, "Lame"]]},
  {True, True, True, True, True, True, True}, TestID -> "el-stiffness-not-positive-definite-or-missing-refused"]
VerificationTest[
  Module[{lam},
    {elFails[Append[elBase, "Lame" -> {lam, 1}]], elFails[Append[elBase, "Lame" -> {Infinity, 1}]], elFails[Append[elBase, "Lame" -> {1 + I, 1}]],
     elFails[Append[elBase, "BodyForce" -> {lam, 0}]], elFails[Append[elBase, "BodyForce" -> {0, 0, 1}]], elFails[Append[elBase, "Degree" -> 3]],
     elFails[Append[elBase, "Cells" -> {4, 0}]], elFails[Append[elBase, "Domain" -> {{0, 1}, {1, 1}}]], elFails[Append[elBase, "TractionFreeFaces" -> {{3, "Lower"}}]]}],
  ConstantArray[True, 9], TestID -> "el-symbolic-or-malformed-model-refused"]
VerificationTest[
  {elFails[Append[elBase, "BoundaryData" -> Function[{x, y}, {1/x, 0}]]], elFails[Append[elBase, "BoundaryData" -> Function[{x, y}, {I x, 0}]]],
   elFails[Append[elBase, "BoundaryData" -> Function[{x, y}, x]]], elFails[Append[elBase, "BoundaryData" -> ConstantArray[0., {2, 4, 4}]]],
   elFails[Append[elBase, "BoundaryData" -> {Infinity, 0}]]},
  {True, True, True, True, True}, TestID -> "el-non-finite-or-malformed-boundary-data-refused"]
VerificationTest[
  {elFails[elBase, Function[{x, y}, {Sqrt[x - 2], 0}]], elFails[elBase, Function[{x, y}, x]], elFails[elBase, ConstantArray[0., {2, 4, 4}]],
   elFails[Append[elBase, "StressField" -> <|"X" -> ConstantArray[0., {4, 4, 4}], "Parameters" -> {}|>]],
   elFails[Append[elBase, "StressField" -> <|"X" -> ConstantArray[0., {4, 5, 5}], "Parameters" -> {1.}|>]],
   elFails[Append[elBase, "StressField" -> "no such directory"]]},
  {True, True, True, True, True, True}, TestID -> "el-bad-candidate-or-stress-field-refused"]
VerificationTest[
  With[{f = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", <|"Lame" -> {1, 1}, "Domain" -> {{0, 1}, {0, 1}, {0, 1}}, "Cells" -> 4|>]]},
    {FailureQ[f], f["Tag"], StringContainsQ[f["MessageTemplate"], "three-dimensional"]}],
  {True, "NotImplemented", True}, TestID -> "el-three-dimensional-domain-refused-with-a-clear-message"]
(* machine overflow in the certified evaluation is detected, not rounded away *)
VerificationTest[
  With[{f = CertifiedVariationalBound[CertiFieldModel["LinearElasticity", Append[elBase, "BoundaryData" -> Function[{x, y}, {10.^200 x, 0}]]]]},
    {FailureQ[f], f["Tag"]}],
  {True, "Inconclusive"}, TestID -> "el-overflow-is-inconclusive"]
VerificationTest[
  Module[{spec = CertiField`Private`linearElasticitySpec[elBase], x},
    x = ConstantArray[0., {4, 5, 5}]; x[[1, 3, 3]] = 10.^160;
    Catch[CertiField`Private`elLowerBound[spec, <|"X" -> x, "Parameters" -> {}|>], "CertiField"]["Tag"]],
  "Inconclusive", TestID -> "el-overflow-in-the-lower-bound-is-inconclusive"]
