(* ::Package:: *)

(* ::Title:: *)
(*CertiField: linear elasticity on boxes in three dimensions*)

(* ::Text:: *)
(*The model LinearElasticity takes this path when "Domain" has three intervals.  Wolfram Language port of*)
(*Reference/Elasticity/elasticity3d.py; statements and proofs in docs/ELASTICITY_3D.md.  Conventions as in*)
(*Kernel/Elasticity.wl.*)
(**)
(*Model.  E[u] = Integrate[e(u).D.e(u)/2 - b.u] over a box, e(u) = (d1 u1, d2 u2, d3 u3, d2 u3 + d3 u2, d1 u3 + d3 u1,*)
(*d1 u2 + d2 u1), D a constant symmetric positive definite 6x6 matrix with rational entries, b a constant rational body*)
(*force, u = g on the Dirichlet faces, no traction on the other faces.  The Dirichlet datum of the certified problem is*)
(*the trace of the continuous function that is a polynomial of degree p (1 or 2) in each variable on every cell and*)
(*takes the nodal boundary values.*)
(**)
(*Lower bound (Prager-Synge).  The stress field is s_ij = eps_ikl eps_jmn d_k d_m Phi_ln (Beltrami) of six functions*)
(*Phi11, Phi22, Phi33, Phi23, Phi13, Phi12 in the space of C^1 piecewise tricubic functions (eight nodal values per*)
(*node: the derivatives d1^d1 d2^d2 d3^d3 Phi, d in {0, 1}^3, stored with the type t = d1 + 2 d2 + 4 d3), plus*)
(*-diag(b1 x1, b2 x2, b3 x3).  Any machine numbers give a stress field that satisfies div s + b = 0 and has a*)
(*single-valued traction.  At the nodes of the traction-free faces the nodal values that enter the traction are those of*)
(*polynomials psi_b + Sum[P_k psi_k] found in exact rational arithmetic (el3States), with machine numbers P_k: the*)
(*traction vanishes there exactly.  el3LowerBound evaluates the bound, a polynomial of degree two in the nodal values*)
(*with exact rational coefficients.*)
(**)
(*Upper bound.  The energy of the Q_p function with the candidate's nodal values (el3UpperBound).*)
(**)
(*Rigorous: el3UpperBound, el3LowerBound, el3NodalEnclosure, el3States and the exact local matrices.  The two solves*)
(*(el3SolveCandidate, el3SolveStress) are ordinary machine arithmetic and are never trusted.*)

Begin["CertiField`Private`"];

(* ------------------------------------------------------------------------------------------------ *)
(* tables                                                                                            *)
(* ------------------------------------------------------------------------------------------------ *)

(* position of sigma_ij among the six components (11, 22, 33, 23, 13, 12) *)
$el3Voigt = {{1, 6, 5}, {6, 2, 4}, {5, 4, 3}};
(* {m, a, i}: the strain e_m contains d_i u_a *)
$el3Strain = {{1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 3, 2}, {4, 2, 3}, {5, 3, 1}, {5, 1, 3}, {6, 2, 1}, {6, 1, 2}};
(* the Beltrami field: s_m = Sum[c d^alpha Phi_f] over the triples {f, alpha, c}; f = 1, ..., 6 for Phi11, Phi22,
   Phi33, Phi23, Phi13, Phi12 (docs/ELASTICITY_3D.md, Section 3.2) *)
$el3Beltrami = {
  {{2, {0, 0, 2}, 1}, {3, {0, 2, 0}, 1}, {4, {0, 1, 1}, -2}},
  {{1, {0, 0, 2}, 1}, {3, {2, 0, 0}, 1}, {5, {1, 0, 1}, -2}},
  {{1, {0, 2, 0}, 1}, {2, {2, 0, 0}, 1}, {6, {1, 1, 0}, -2}},
  {{1, {0, 1, 1}, -1}, {4, {2, 0, 0}, -1}, {5, {1, 1, 0}, 1}, {6, {1, 0, 1}, 1}},
  {{2, {1, 0, 1}, -1}, {4, {1, 1, 0}, 1}, {5, {0, 2, 0}, -1}, {6, {0, 1, 1}, 1}},
  {{3, {1, 1, 0}, -1}, {4, {1, 0, 1}, 1}, {5, {0, 1, 1}, 1}, {6, {0, 0, 2}, -1}}};

(* scaled one-dimensional tables of a direction with cell length h:
   el3HTab[h][[a + 1, c + 1, k + 1, l + 1]] = h^(1 + Mod[k, 2] + Mod[l, 2] - a - c) Integrate[H_k^(a) H_l^(c)],
   el3JTab[h][[q + 1, a + 1, k + 1]] = h^(1 + Mod[k, 2] - a) Integrate[xi^q H_k^(a)] *)
el3HTab[h_] := Table[Table[h^(1 + Mod[k, 2] + Mod[l, 2] - a - c), {k, 0, 3}, {l, 0, 3}] elIMat[$elHermite, a, $elHermite, c], {a, 0, 2}, {c, 0, 2}];
el3JTab[h_] := Table[h^(1 + Mod[k, 2] - a) elPInt[elPMul[UnitVector[q + 1, q + 1], elPDer[$elHermite[[k + 1]], a]]], {q, 0, 1}, {a, 0, 2}, {k, 0, 3}];

(* the local unknowns of the stress functions: {f, k1, k2, k3}, k_i = 2 (end of the cell in direction i) + (1 for a
   derivative in direction i); local index 64 (f - 1) + 16 k1 + 4 k2 + k3 + 1 *)
$el3Dofs = Flatten[Table[{f, k1, k2, k3}, {f, 6}, {k1, 0, 3}, {k2, 0, 3}, {k3, 0, 3}], 3];
(* after the affine part of each function at the lower corner is subtracted, the value and the three first derivatives
   at that corner are zero and are left out *)
$el3Keep = Flatten[Position[$el3Dofs, {_, k1_, k2_, k3_} /; ! MemberQ[{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {k1, k2, k3}], {1}, Heads -> False]];

(* exact element matrix A of Integrate[s.S.s] over a cell for the stress of the 384 local basis functions.  The stress of
   an affine function is zero, so A annihilates the nodal values of 1, x1, x2, x3 in each function; the certified
   evaluation relies on this (it works with differences), so it is checked here in exact arithmetic, with the
   symmetry of A. *)
el3StressA[h_List, sm_List] := el3StressA[h, sm] = Module[{th = el3HTab /@ h, blocks, am},
  blocks = Table[ConstantArray[0, {64, 64}], {6}, {6}];
  Do[
    If[sm[[m, n]] != 0,
      Do[
        blocks[[t1[[1]], t2[[1]]]] += sm[[m, n]] t1[[3]] t2[[3]] KroneckerProduct[
          th[[1, t1[[2, 1]] + 1, t2[[2, 1]] + 1]], th[[2, t1[[2, 2]] + 1, t2[[2, 2]] + 1]], th[[3, t1[[2, 3]] + 1, t2[[2, 3]] + 1]]],
        {t1, $el3Beltrami[[m]]}, {t2, $el3Beltrami[[n]]}]],
    {m, 6}, {n, 6}];
  am = ArrayFlatten[blocks];
  If[am =!= Transpose[am] || ! AllTrue[Flatten[Normal[el3Affine[h] . am]], # == 0 &],
    Throw[failure["Internal", "The element matrix of the stress functions is not symmetric or does not annihilate affine functions."], "CertiField"]];
  <|"A" -> am, "Rows" -> el3Rows[am[[$el3Keep, $el3Keep]]]|>];

(* the vectors w with Integrate[s.S.s_b] = Sum[X_k (w_k1 + xc1 w_k2 + xc2 w_k3 + xc3 w_k4)] over the local unknowns,
   xc the lower corner of the cell, s_b = (-b1 x1, -b2 x2, -b3 x3, 0, 0, 0); 384 x 4 *)
el3StressW[h_List, sm_List, b_List] := el3StressW[h, sm, b] = Module[{jt = el3JTab /@ h, w = ConstantArray[0, {4, 6, 64}]},
  Do[
    If[sm[[m, n]] != 0 && b[[n]] != 0,
      With[{f = t1[[1]], al = t1[[2]], cc = t1[[3]] sm[[m, n]] (-b[[n]])},
        w[[n + 1, f]] += cc Flatten[Outer[Times, jt[[1, 1, al[[1]] + 1]], jt[[2, 1, al[[2]] + 1]], jt[[3, 1, al[[3]] + 1]]]];
        w[[1, f]] += cc h[[n]] Flatten[Outer[Times, jt[[1, If[n == 1, 2, 1], al[[1]] + 1]], jt[[2, If[n == 2, 2, 1], al[[2]] + 1]], jt[[3, If[n == 3, 2, 1], al[[3]] + 1]]]]]],
    {m, 6}, {t1, $el3Beltrami[[m]]}, {n, 3}];
  w = Transpose[Flatten /@ w];
  If[! AllTrue[Flatten[Normal[el3Affine[h] . w]], # == 0 &],
    Throw[failure["Internal", "The load vectors of the stress functions do not annihilate affine functions."], "CertiField"]];
  w];

(* the nodal values of the affine functions 1, x1 - xc1, x2 - xc2, x3 - xc3 in each of the six functions (24 x 384) *)
el3Affine[h_List] := SparseArray[Flatten[Table[
    Table[If[d[[1]] != f, 0,
        With[{e = Quotient[Rest[d], 2], dd = Mod[Rest[d], 2]},
          Which[dd === {0, 0, 0}, If[r == 0, 1, h[[r]] e[[r]]], r > 0 && dd === UnitVector[3, r], 1, True, 0]]], {d, $el3Dofs}],
    {f, 6}, {r, 0, 3}], 1]];

(* exact element matrix of Integrate[e.D.e] over a cell and exact load vector for the local unknowns {a, i1, i2, i3}:
   component a at the local node (i1, i2, i3) of the Q_p element *)
el3StiffnessLocal[p_Integer, h_List, dm_List, b_List] := el3StiffnessLocal[p, h, dm, b] = Module[
  {lag = elLagrange[p], il, nl = (p + 1)^3, vol = Times @@ h, blocks, km, dofs, intL, keep},
  il = Table[elIMat[lag, da, lag, db], {da, 0, 1}, {db, 0, 1}];
  blocks = Table[ConstantArray[0, {nl, nl}], {3}, {3}];
  Do[
    If[dm[[s1[[1]], s2[[1]]]] != 0,
      blocks[[s1[[2]], s2[[2]]]] += dm[[s1[[1]], s2[[1]]]] vol/(h[[s1[[3]]]] h[[s2[[3]]]]) *
        KroneckerProduct @@ Table[il[[Boole[s1[[3]] == t] + 1, Boole[s2[[3]] == t] + 1]], {t, 3}]],
    {s1, $el3Strain}, {s2, $el3Strain}];
  km = ArrayFlatten[blocks];
  dofs = Flatten[Table[{a, i1, i2, i3}, {a, 3}, {i1, 0, p}, {i2, 0, p}, {i3, 0, p}], 3];
  intL = elPInt /@ lag;
  (* symmetric, and a constant displacement has no strain: checked, because the evaluation works with differences *)
  If[km =!= Transpose[km] || ! AllTrue[Flatten[Table[km . Boole[# == a & /@ dofs[[All, 1]]], {a, 3}]], # == 0 &],
    Throw[failure["Internal", "The element matrix is not symmetric or does not annihilate constants."], "CertiField"]];
  keep = Flatten[Position[dofs, {_, i1_, i2_, i3_} /; {i1, i2, i3} =!= {0, 0, 0}, {1}, Heads -> False]];
  <|"Dofs" -> dofs, "K" -> km, "Keep" -> keep, "Rows" -> el3Rows[km[[keep, keep]]],
    "f" -> Table[b[[d[[1]]]] vol intL[[d[[2]] + 1]] intL[[d[[3]] + 1]] intL[[d[[4]] + 1]], {d, dofs}]|>];

(* ------------------------------------------------------------------------------------------------ *)
(* rigorous evaluation of quadratic forms, cell by cell                                               *)
(* ------------------------------------------------------------------------------------------------ *)

(* enclosures {lo, hi} of a list of exact numbers, as two lists of machine numbers: the machine number nearest to each,
   moved outward by rDn and rUp, and checked by exact comparison; an entry that fails the check gets fLo and fHi *)
el3Enclose[q_List] := Module[{x = Quiet@Developer`ToPackedArray[N[q]], lo, hi, okl, okh},
  If[q === {}, Return[{{}, {}}, Module]];
  If[! Developer`PackedArrayQ[x, Real], Return[Transpose[ivC /@ q], Module]];
  lo = rDn[x]; hi = rUp[x];
  okl = Thread[SetPrecision[lo, Infinity] <= q]; okh = Thread[SetPrecision[hi, Infinity] >= q];
  {MapThread[If[#1, #2, fLo[#3]] &, {okl, lo, q}], MapThread[If[#1, #2, fHi[#3]] &, {okh, hi, q}]}];

(* the quadratic form Sum[m[[k, l]] d_k d_l] of an exact symmetric matrix m, prepared row by row: {k, enclosure of
   m[[k, k]] or None, the columns l > k with m[[k, l]] != 0, lower and upper ends of the enclosures of 2 m[[k, l]]} *)
el3Rows[m_] := Module[{n = Length[m], vals, enc},
  vals = DeleteCases[DeleteDuplicates[Join[Diagonal[m], 2 Flatten[Table[m[[k, k + 1 ;;]], {k, n - 1}]]]], 0];
  enc = AssociationThread[vals -> Transpose[el3Enclose[vals]]];
  Table[
    With[{cols = k + Flatten[Position[m[[k, k + 1 ;;]], x_ /; x != 0, {1}, Heads -> False]]},
      If[m[[k, k]] == 0 && cols === {}, Nothing,
        {k, If[m[[k, k]] == 0, None, enc[m[[k, k]]]], cols,
         Developer`ToPackedArray[N[enc[2 m[[k, #]]][[1]] & /@ cols]], Developer`ToPackedArray[N[enc[2 m[[k, #]]][[2]] & /@ cols]]}]],
    {k, n}]];

(* enclosure of the sum over the first dimension of an interval array {lo, hi}: pairwise, every addition rounded
   outward and checked *)
el3TreeSum[{lo_, hi_}] := Module[{l = lo, u = hi, k},
  While[Length[l] > 1,
    k = Quotient[Length[l], 2];
    l = Join[ndChk1[rDn[l[[1 ;; 2 k ;; 2]] + l[[2 ;; 2 k ;; 2]]]], l[[2 k + 1 ;;]]];
    u = Join[ndChk1[rUp[u[[1 ;; 2 k ;; 2]] + u[[2 ;; 2 k ;; 2]]]], u[[2 k + 1 ;;]]]];
  ndChk[{l[[1]], u[[1]]}]];

(* RIGOROUS enclosure, cell by cell, of Sum[m[[k, l]] d_k d_l] (rows from el3Rows) for interval values d = {lo, hi},
   lo and hi packed arrays of dimensions {Length[m], number of cells}.  Row by row: the products with the coefficients
   of a row are formed together, summed pairwise and multiplied by d_k.  Every product and sum is rounded outward,
   and every result is checked to consist of machine numbers (an overflow is detected, not rounded away). *)
el3Quadratic[rows_, {lo_, hi_}] := Module[{tl = {}, th = {}, a, b, p1, p2, p3, p4, s, t},
  Do[
    With[{k = r[[1]], dg = r[[2]], cols = r[[3]], cl = r[[4]], ch = r[[5]]},
      If[dg =!= None, t = ndMul[ivSqr[{lo[[k]], hi[[k]]}], dg]; AppendTo[tl, t[[1]]]; AppendTo[th, t[[2]]]];
      If[cols =!= {},
        a = lo[[cols]]; b = hi[[cols]];
        p1 = ndChk1[cl a]; p2 = ndChk1[ch a]; p3 = ndChk1[cl b]; p4 = ndChk1[ch b];
        s = el3TreeSum[{ndChk1[rDn[vMin[vMin[p1, p2], vMin[p3, p4]]]], ndChk1[rUp[vMax[vMax[p1, p2], vMax[p3, p4]]]]}];
        t = ndMul[{lo[[k]], hi[[k]]}, s]; AppendTo[tl, t[[1]]]; AppendTo[th, t[[2]]]]],
    {r, rows}];
  (* the entries of tl and th are checked packed vectors of machine numbers; stacking them packs them as a matrix *)
  el3TreeSum[{Developer`ToPackedArray[tl], Developer`ToPackedArray[th]}]];

(* ------------------------------------------------------------------------------------------------ *)
(* model normalisation                                                                               *)
(* ------------------------------------------------------------------------------------------------ *)

$el3DefaultCells = {4, 4, 4};

el3DomainQ[params_Association] := MatchQ[exactify[Lookup[params, "Domain", Missing[]]], {{_, _}, {_, _}, {_, _}}];

(* an array of triples over the nodes as a triple of arrays over the nodes *)
el3Components[v_, mm_List] := If[ArrayQ[v, 4] && Dimensions[v] === Append[mm, 3], Transpose[v, {2, 3, 4, 1}], $Failed];

linearElasticitySpec3D[params_Association] := Module[
  {dom, cells, lame, dmat, b, p, nat, dir, h, mm, uaxes, vnodes, ind, fixed, g, gval},
  dom = exactify[Lookup[params, "Domain", Missing[]]];
  If[! (MatchQ[dom, {{_, _}, {_, _}, {_, _}}] && AllTrue[Flatten[dom], MatchQ[#, _Integer | _Rational] &] && AllTrue[dom, #[[1]] < #[[2]] &]),
    Throw[failure["BadModel", "\"Domain\" must be {{a1, b1}, {a2, b2}, {a3, b3}} with rational a < b."], "CertiField"]];
  lame = Lookup[params, "Lame", Missing[]]; dmat = Lookup[params, "StiffnessMatrix", Missing[]];
  Which[
    ! MissingQ[lame] && ! MissingQ[dmat],
      Throw[failure["BadModel", "Give either \"Lame\" or \"StiffnessMatrix\", not both."], "CertiField"],
    ! MissingQ[lame],
      lame = exactify[lame];
      If[! MatchQ[lame, {_Integer | _Rational, _Integer | _Rational}],
        Throw[failure["BadModel", "\"Lame\" must be {lambda, mu}, two rational numbers."], "CertiField"]];
      dmat = Table[Which[i <= 3 && j <= 3, lame[[1]] + If[i == j, 2 lame[[2]], 0], i == j, lame[[2]], True, 0], {i, 6}, {j, 6}],
    ! MissingQ[dmat], dmat = exactify[dmat],
    True, Throw[failure["BadModel", "The material is missing: give \"Lame\" -> {lambda, mu} or \"StiffnessMatrix\" -> a 6x6 matrix."], "CertiField"]];
  If[! (MatrixQ[dmat, MatchQ[#, _Integer | _Rational] &] && Dimensions[dmat] === {6, 6} && dmat === Transpose[dmat]),
    Throw[failure["BadModel", "\"StiffnessMatrix\" must be a symmetric 6x6 matrix of rational numbers (strains in the order d1 u1, d2 u2, d3 u3, d2 u3 + d3 u2, d1 u3 + d3 u1, d1 u2 + d2 u1)."], "CertiField"]];
  (* positive definite: the leading principal minors are positive (exact arithmetic) *)
  If[! AllTrue[Range[6], Det[dmat[[1 ;; #, 1 ;; #]]] > 0 &],
    Throw[failure["BadModel", "The stiffness is not positive definite (for \"Lame\" -> {lambda, mu}: mu > 0 and 3 lambda + 2 mu > 0 are needed)."], "CertiField"]];
  b = exactify[Lookup[params, "BodyForce", {0, 0, 0}]];
  If[! MatchQ[b, {_Integer | _Rational, _Integer | _Rational, _Integer | _Rational}],
    Throw[failure["BadModel", "\"BodyForce\" must be {b1, b2, b3}, three rational numbers (a constant body force)."], "CertiField"]];
  p = Lookup[params, "Degree", 1];
  If[! MatchQ[p, 1 | 2], Throw[failure["BadModel", "\"Degree\" must be 1 or 2."], "CertiField"]];
  cells = Replace[Lookup[params, "Cells", Automatic], Automatic -> $el3DefaultCells];
  If[IntegerQ[cells], cells = {cells, cells, cells}];
  If[! (VectorQ[cells, IntegerQ] && Length[cells] == 3 && Min[cells] >= 1),
    Throw[failure["BadModel", "\"Cells\" must be a positive integer or a list of three positive integers."], "CertiField"]];
  nat = Lookup[params, "TractionFreeFaces", {}];
  If[! MatchQ[nat, {{_Integer, "Lower" | "Upper"} ...}] || ! AllTrue[nat, 1 <= #[[1]] <= 3 &],
    Throw[failure["BadModel", "\"TractionFreeFaces\" must be a list of {direction, \"Lower\" | \"Upper\"}."], "CertiField"]];
  dir = Table[{! MemberQ[nat, {i, "Lower"}], ! MemberQ[nat, {i, "Upper"}]}, {i, 3}];          (* True: Dirichlet *)
  If[! MemberQ[Flatten[dir], True],
    Throw[failure["BadModel", "At least one face must be a Dirichlet face."], "CertiField"]];
  h = (dom[[All, 2]] - dom[[All, 1]])/cells;
  mm = p cells + 1;                                                             (* displacement nodes per direction *)
  uaxes = Table[N[dom[[i, 1]] + j h[[i]]/p], {i, 3}, {j, 0, p cells[[i]]}];
  vnodes = Table[dom[[i, 1]] + j h[[i]], {i, 3}, {j, 0, cells[[i]]}];           (* exact nodes of the grid of the stress functions *)
  ind = Table[ReplacePart[ConstantArray[1, mm[[i]]], {If[dir[[i, 1]], 1 -> 0, Nothing], If[dir[[i, 2]], -1 -> 0, Nothing]}], {i, 3}];
  fixed = 1 - Outer[Times, ind[[1]], ind[[2]], ind[[3]]];
  (* boundary data: three numbers, a function of the coordinates returning three numbers, or nodal arrays {g1, g2, g3} *)
  g = Lookup[params, "BoundaryData", {0, 0, 0}];
  gval = Which[
    MatchQ[g, {_?NumericQ, _?NumericQ, _?NumericQ}], Table[ConstantArray[g[[a]], mm], {a, 3}],
    ArrayQ[g, 4] && Dimensions[g] === Prepend[mm, 3], g,
    True, el3Components[Quiet@Apply[g, Outer[List, uaxes[[1]], uaxes[[2]], uaxes[[3]]], {3}], mm]];
  (* only the values at the Dirichlet nodes are used: a function may be undefined or singular elsewhere *)
  gval = If[ArrayQ[gval, 4] && Dimensions[gval] === Prepend[mm, 3],
    With[{fx = Flatten[fixed]},
      ndRealArray[Table[ArrayReshape[ReplacePart[ConstantArray[0, Length[fx]], Thread[Flatten[Position[fx, 1]] -> Pick[Flatten[gval[[a]]], fx, 1]]], mm], {a, 3}], Prepend[mm, 3]]],
    $Failed];
  If[gval === $Failed,
    Throw[failure["BadModel", "\"BoundaryData\" must be three real numbers, a function of the coordinates returning three real numbers, or nodal arrays of dimensions {3, M1, M2, M3} with M = Degree Cells + 1."], "CertiField"]];
  <|"n" -> 3, "Domain" -> dom, "Cells" -> cells, "h" -> h, "p" -> p, "D" -> dmat, "S" -> Inverse[dmat], "b" -> b, "Dirichlet" -> dir,
    "M" -> mm, "UAxes" -> uaxes, "VNodes" -> vnodes, "Fixed" -> fixed, "G" -> gval|>];

(* ------------------------------------------------------------------------------------------------ *)
(* the traction-free faces: the polynomials psi_b, psi_1, ..., psi_r (docs/ELASTICITY_3D.md, 3.3-3.4) *)
(* ------------------------------------------------------------------------------------------------ *)

(* the pairs {f, m}: d_i^m Phi_f enters the traction on a plane normal to e_i *)
el3Relevant[i_] := Union[Flatten[Table[{#[[1]], #[[2, i]]} & /@ $el3Beltrami[[$el3Voigt[[i, j]]]], {j, 3}], 1]];
(* the coefficient of x1^g1 x2^g2 x3^g3 in function f *)
el3Col[f_, {g1_, g2_, g3_}] := 64 (f - 1) + 16 g1 + 4 g2 + g3 + 1;
(* the axis a if the two faces normal to e_a are Dirichlet faces and the four others are traction-free, else None *)
el3RingAxis[dir_] := Module[{c = Select[Range[3], dir[[#]] === {True, True} && AllTrue[Complement[Range[3], {#}], Function[t, dir[[t]] === {False, False}]] &]},
  If[c === {}, None, First[c]]];

(* force and moment about the origin of the traction s e_axis of the sextuple of polynomials with coefficient vector v,
   on the cross-section at the lower end of the axis; exact *)
el3Resultants[v_, axis_, dom_] := Module[{x1, x2, x3, xs, cf, polys, s, trac, tt, c0, int},
  xs = {x1, x2, x3};
  cf = ArrayReshape[v, {6, 4, 4, 4}];
  polys = Table[Sum[cf[[f, g1 + 1, g2 + 1, g3 + 1]] x1^g1 x2^g2 x3^g3, {g1, 0, 3}, {g2, 0, 3}, {g3, 0, 3}], {f, 6}];
  s = Table[Total[(#[[3]] D[polys[[#[[1]]]], {x1, #[[2, 1]]}, {x2, #[[2, 2]]}, {x3, #[[2, 3]]}]) & /@ $el3Beltrami[[m]]], {m, 6}];
  c0 = dom[[axis, 1]]; tt = Complement[Range[3], {axis}];
  int[e_] := Total[(#2 Product[If[t == axis, 1, (dom[[t, 2]]^(#1[[t]] + 1) - dom[[t, 1]]^(#1[[t]] + 1))/(#1[[t]] + 1)], {t, 3}]) & @@@
    CoefficientRules[Expand[e /. xs[[axis]] -> c0], xs]];
  trac = Table[s[[$el3Voigt[[axis, j]]]], {j, 3}];
  Join[int /@ trac, Table[Sum[Signature[{i, j, k}] int[xs[[j]] trac[[k]]], {j, 3}, {k, 3}], {i, 3}]]];

(* <|"Positions" -> the nodal values {f, t + 1, i1, i2, i3} fixed by the traction-free faces (lexicographic order),
     "Tb" -> their exact values for psi_b, "T" -> the same for psi_1, ..., psi_r|>.
   The traction of psi = Sum[c[[el3Col[f, g]]] x^g Phi_f] on a face {x_i = a} is a polynomial in the two other
   variables; "it equals minus the traction of s_b" is linear in c.  The reduced row echelon form is unique: psi_b is
   the solution with the free unknowns zero, and there is one null vector for each free unknown, in increasing order.
   For a ring of four traction-free faces psi_1, ..., psi_6 are the first null vectors whose resultants (force and
   moment through the cross-section) are independent; otherwise r = 0.  This is the rule of the Python reference, so
   the parameters of a field mean the same in both. *)
el3States[spec_] := el3StatesData[spec["Domain"], spec["Cells"], spec["Dirichlet"], spec["b"]];

el3StatesData[dom_, cells_, dir_, b_] := el3StatesData[dom, cells, dir, b] = Module[
  {faces, rows = {}, rhs = {}, rr, piv, free, part, axis, kept = {}, states = {}, mask, pos, vnodes, pw, vals},
  faces = Select[Tuples[{Range[3], {1, 2}}], ! dir[[#[[1]], #[[2]]]] &];
  If[faces === {}, Return[<|"Positions" -> {}, "Tb" -> {}, "T" -> {}|>, Module]];
  Do[
    With[{i = fc[[1]], a = dom[[fc[[1]], fc[[2]]]], tt = Complement[Range[3], {fc[[1]]}]},
      Do[
        Module[{row = ConstantArray[0, 384]},
          Do[
            With[{f = term[[1]], al = term[[2]], c = term[[3]], g1 = r1 + term[[2, tt[[1]]]], g2 = r2 + term[[2, tt[[2]]]]},
              If[g1 <= 3 && g2 <= 3,
                Do[
                  row[[el3Col[f, ReplacePart[{0, 0, 0}, {i -> gi, tt[[1]] -> g1, tt[[2]] -> g2}]]]] +=
                    c FactorialPower[g1, al[[tt[[1]]]]] FactorialPower[g2, al[[tt[[2]]]]] FactorialPower[gi, al[[i]]] elPow[a, gi - al[[i]]],
                  {gi, al[[i]], 3}]]],
            {term, $el3Beltrami[[$el3Voigt[[i, j]]]]}];
          AppendTo[rows, row]; AppendTo[rhs, If[j == i && r1 == 0 && r2 == 0, b[[i]] a, 0]]],
        {j, 3}, {r1, 0, 3}, {r2, 0, 3}]],
    {fc, faces}];
  rr = Select[RowReduce[MapThread[Append, {rows, rhs}]], ! AllTrue[#, # == 0 &] &];
  piv = First[FirstPosition[#, x_ /; x != 0]] & /@ rr;
  If[MemberQ[piv, 385],
    Throw[failure["BadModel", "The body force admits no polynomial equilibrium field with these traction-free faces."], "CertiField"]];
  free = Complement[Range[384], piv];
  part = ConstantArray[0, 384];
  Do[part[[piv[[k]]]] = rr[[k, 385]], {k, Length[piv]}];
  axis = el3RingAxis[dir];
  If[axis =!= None,
    Do[
      If[Length[kept] < 6,
        Module[{v = UnitVector[384, fr], res},
          Do[v[[piv[[k]]]] = -rr[[k, fr]], {k, Length[piv]}];
          res = el3Resultants[v, axis, dom];
          If[MatrixRank[Append[kept, res]] > Length[kept], AppendTo[kept, res]; AppendTo[states, v]]]],
      {fr, free}]];
  (* the nodal values fixed by each face: the pairs {f, m} of el3Relevant, all types t with d_i = m *)
  mask = ConstantArray[0, Join[{6, 8}, cells + 1]];
  Do[
    With[{i = fc[[1]], sl = If[fc[[2]] == 1, 1, cells[[fc[[1]]]] + 1]},
      Do[
        If[BitGet[t, i - 1] == fm[[2]],
          Switch[i,
            1, mask[[fm[[1]], t + 1, sl, All, All]] = 1,
            2, mask[[fm[[1]], t + 1, All, sl, All]] = 1,
            3, mask[[fm[[1]], t + 1, All, All, sl]] = 1]],
        {fm, el3Relevant[i]}, {t, 0, 7}]],
    {fc, faces}];
  pos = Position[mask, 1];
  (* nodal values of a sextuple: pw[[i, d + 1]] is the matrix (node, g) of d^d/dx_i^d x_i^g at the nodes of direction i *)
  vnodes = Table[dom[[i, 1]] + j (dom[[i, 2]] - dom[[i, 1]])/cells[[i]], {i, 3}, {j, 0, cells[[i]]}];
  pw = Table[If[g >= d, FactorialPower[g, d] elPow[x, g - d], 0], {i, 3}, {d, 0, 1}, {x, vnodes[[i]]}, {g, 0, 3}];
  vals[v_] := With[{cf = ArrayReshape[v, {6, 4, 4, 4}]},
    Extract[Table[
      With[{p1 = pw[[1, BitGet[t, 0] + 1]], p2 = pw[[2, BitGet[t, 1] + 1]], p3 = pw[[3, BitGet[t, 2] + 1]]},
        Transpose[Transpose[(p1 . cf[[f]]) . Transpose[p3], {1, 3, 2}] . Transpose[p2], {1, 3, 2}]],
      {f, 6}, {t, 0, 7}], pos]];
  <|"Positions" -> pos, "Tb" -> vals[part], "T" -> (vals /@ states)|>];

(* ------------------------------------------------------------------------------------------------ *)
(* stress fields: <|"X" -> array {6, 8, N1 + 1, N2 + 1, N3 + 1} of nodal values (the entries fixed by the  *)
(* traction-free faces are ignored), "Parameters" -> {P_1, ..., P_r}|>                               *)
(* ------------------------------------------------------------------------------------------------ *)

el3PackField[spec_, f_] := Module[{dims = Join[{6, 8}, spec["Cells"] + 1], r = Length[el3States[spec]["T"]], x, par},
  If[! AssociationQ[f], Return[$Failed, Module]];
  x = ndRealArray[Lookup[f, "X", $Failed], dims];
  par = Lookup[f, "Parameters", {}];
  par = If[r == 0, If[ListQ[par] && Length[par] == 0, {}, $Failed], ndRealArray[par, {r}]];
  If[x === $Failed || par === $Failed, $Failed, <|"X" -> x, "Parameters" -> par|>]];

el3LinearIndex[pos_, dims_] := (pos - 1) . Reverse[FoldList[Times, 1, Reverse[Rest[dims]]]] + 1;

(* RIGOROUS: the nodal values as a pair {lo, hi} of packed arrays.  Free entries: the stored machine numbers (exact).
   Entries fixed by the traction-free faces: enclosures of the exact rational values of psi_b + Sum[P_k psi_k],
   computed in rational arithmetic from the machine numbers P_k. *)
el3NodalEnclosure[spec_, f_] := Module[{st = el3States[spec], dims = Dimensions[f["X"]], lo, hi, par, vals, uv, enc, idx, e},
  lo = Flatten[f["X"]]; hi = lo;
  If[st["Positions"] =!= {},
    par = exactNumber /@ Normal[f["Parameters"]];
    vals = If[par === {}, st["Tb"], st["Tb"] + par . st["T"]];
    uv = DeleteDuplicates[vals];
    enc = AssociationThread[uv -> Transpose[el3Enclose[uv]]];
    e = Lookup[enc, vals];
    If[! VectorQ[Flatten[e], MachineNumberQ], ndOverflow[]];
    idx = el3LinearIndex[st["Positions"], dims];
    (* enclosures of exact numbers (inputs, not results of machine arithmetic): packed so that the arrays stay packed *)
    lo[[idx]] = Developer`ToPackedArray[e[[All, 1]]]; hi[[idx]] = Developer`ToPackedArray[e[[All, 2]]]];
  ndChk[{ArrayReshape[lo, dims], ArrayReshape[hi, dims]}]];

(* ------------------------------------------------------------------------------------------------ *)
(* the two certified bounds                                                                          *)
(* ------------------------------------------------------------------------------------------------ *)

(* RIGOROUS enclosure {lo, hi} (exact rationals) of the energy of the Q_p function with the nodal values u0
   (dimensions {3, M1, M2, M3}).  As in two dimensions: the element matrix annihilates constants, so it is applied to
   the differences from the value at the lower corner of the cell, component by component. *)
el3UpperBound[spec_, u0_] := ndQuiet@Module[{p = spec["p"], n1, n2, n3, u, loc, val, keep, dl, dh, en},
  {n1, n2, n3} = spec["Cells"];
  u = ndRealArray[u0, Prepend[spec["M"], 3]];
  If[u === $Failed, Throw[failure["BadCandidate", "The nodal values must be a real array of dimensions {3, M1, M2, M3} with M = Degree Cells + 1."], "CertiField"]];
  If[! AllTrue[Range[3], Max[Abs[spec["Fixed"] (u[[#]] - spec["G"][[#]])]] == 0. &],
    Throw[failure["Internal", "The Dirichlet nodal values are not imposed exactly."], "CertiField"]];
  loc = el3StiffnessLocal[p, spec["h"], spec["D"], spec["b"]];
  val[{a_, i1_, i2_, i3_}] := ivPt[Flatten[u[[a, i1 + 1 ;; i1 + 1 + p (n1 - 1) ;; p, i2 + 1 ;; i2 + 1 + p (n2 - 1) ;; p, i3 + 1 ;; i3 + 1 + p (n3 - 1) ;; p]]]];
  keep = loc["Dofs"][[loc["Keep"]]];
  {dl, dh} = Developer`ToPackedArray /@ Transpose[Table[ndSub[val[d], val[{d[[1]], 0, 0, 0}]], {d, keep}]];
  en = ivK[el3Quadratic[loc["Rows"], {dl, dh}], 1/2];
  (* the load: f.u = f.d + (sum of f over the nodes of the component) (value at the lower corner) *)
  Do[With[{fj = loc["f"][[loc["Keep"][[j]]]]}, If[fj != 0, en = ndSub[en, ivK[{dl[[j]], dh[[j]]}, fj]]]], {j, Length[keep]}];
  Do[With[{fa = Total[Pick[loc["f"], loc["Dofs"][[All, 1]], a]]}, If[fa != 0, en = ndSub[en, ivK[val[{a, 0, 0, 0}], fa]]]], {a, 3}];
  {exactTotal[en[[1]]], exactTotal[en[[2]]]}];

(* the face cells of the face {x_i = (pos-th node)} of an array over the nodes of a grid: a matrix over the two other
   directions, in increasing order *)
el3FaceSlice[arr_, i_, pos_] := Switch[i, 1, arr[[pos, All, All]], 2, arr[[All, pos, All]], 3, arr[[All, All, pos]]];

(* RIGOROUS lower bound  Integrate[(s nu).g over the Dirichlet faces] - 1/2 Integrate[s.S.s]  for the stress field of
   the given nodal values and parameters.  Returns <|"Lower" -> the bound, "Boundary" -> enclosure of the first term,
   "IntQ" -> enclosure of Integrate[s.S.s]|> (exact rationals).  Sums over cells and cell faces are exact. *)
el3LowerBound[spec_, field_] := ndQuiet@Module[
  {f, n1, n2, n3, nc, h = spec["h"], sm = spec["S"], b = spec["b"], dom = spec["Domain"], p = spec["p"], g = spec["G"], vn = spec["VNodes"],
   sa, w, xlo, xhi, cellv, val, keep, dl, dh, q, lin = {}, xc, qb, intq, bnd = {0, 0}, hl, li},
  f = el3PackField[spec, field];
  If[f === $Failed,
    Throw[failure["BadField", "The stress field must be <|\"X\" -> nodal values of dimensions {6, 8, N1 + 1, N2 + 1, N3 + 1}, \"Parameters\" -> the parameters of the traction-free faces|> for the grid of the model."], "CertiField"]];
  {n1, n2, n3} = spec["Cells"]; nc = n1 n2 n3;
  sa = el3StressA[h, sm];
  {xlo, xhi} = el3NodalEnclosure[spec, f];
  (* the values of the local unknown {f, k1, k2, k3} on all cells (C order of the cells) *)
  cellv[arr_, {ff_, k1_, k2_, k3_}] := Flatten[arr[[ff, Mod[k1, 2] + 2 Mod[k2, 2] + 4 Mod[k3, 2] + 1,
      Quotient[k1, 2] + 1 ;; Quotient[k1, 2] + n1, Quotient[k2, 2] + 1 ;; Quotient[k2, 2] + n2, Quotient[k3, 2] + 1 ;; Quotient[k3, 2] + n3]]];
  val[d_] := {cellv[xlo, d], cellv[xhi, d]};
  (* differences from the affine function with the value and the first derivatives of the lower corner of the cell,
     function by function.  Its stress is zero, so X.A.X = d.A.d and w.X = w.d (checked exactly in el3StressA and
     el3StressW).  The value at the corner is subtracted first: each subtraction is then rounded relative to a small
     number. *)
  keep = $el3Dofs[[$el3Keep]];
  {dl, dh} = Developer`ToPackedArray /@ Transpose[Table[
    With[{ff = d[[1]], e = Quotient[Rest[d], 2], dd = Mod[Rest[d], 2]},
      Which[
        dd === {0, 0, 0},
          Fold[If[e[[#2]] == 1, ndSub[#1, ivK[val[Join[{ff}, UnitVector[3, #2]]], h[[#2]]]], #1] &, ndSub[val[d], val[{ff, 0, 0, 0}]], {1, 2, 3}],
        Total[dd] == 1, ndSub[val[d], val[Join[{ff}, dd]]],
        True, val[d]]],
    {d, keep}]];
  (* the cell integrals: Integrate[s.S.s] = (stress functions) + 2 (cross term with s_b) + (s_b part) *)
  q = el3Quadratic[sa["Rows"], {dl, dh}];
  If[b =!= {0, 0, 0},
    w = el3StressW[h, sm, b];
    xc = Table[With[{e = ivC /@ vn[[i, 1 ;; spec["Cells"][[i]]]]},
        Table[Developer`ToPackedArray[N[Flatten[Switch[i,
            1, Table[e[[c1, s]], {c1, n1}, {c2, n2}, {c3, n3}],
            2, Table[e[[c2, s]], {c1, n1}, {c2, n2}, {c3, n3}],
            3, Table[e[[c3, s]], {c1, n1}, {c2, n2}, {c3, n3}]]]]], {s, 2}]], {i, 3}];
    Do[
      With[{wj = w[[$el3Keep[[j]]]]},
        If[wj =!= {0, 0, 0, 0},
          Module[{coef = ndChk[ivPt[ConstantArray[0., nc]] + ivC[wj[[1]]]]},
            Do[coef = ndAdd[coef, ndMul[xc[[i]], ivC[wj[[i + 1]]]]], {i, 3}];
            AppendTo[lin, ivK[ndMul[{dl[[j]], dh[[j]]}, coef], 2]]]]],
      {j, Length[keep]}];
    If[lin =!= {}, q = ndAdd[q, el3TreeSum[{Developer`ToPackedArray[lin[[All, 1]]], Developer`ToPackedArray[lin[[All, 2]]]}]]]];
  qb = Sum[sm[[i, j]] b[[i]] b[[j]] If[i == j,
      (dom[[i, 2]]^3 - dom[[i, 1]]^3)/3 Product[dom[[t, 2]] - dom[[t, 1]], {t, Complement[Range[3], {i}]}],
      (dom[[i, 2]]^2 - dom[[i, 1]]^2)/2 (dom[[j, 2]]^2 - dom[[j, 1]]^2)/2 Product[dom[[t, 2]] - dom[[t, 1]], {t, Complement[Range[3], {i, j}]}]],
    {i, 3}, {j, 3}];
  intq = {exactTotal[q[[1]]], exactTotal[q[[2]]]} + qb;
  (* the boundary term.  On a Dirichlet face {x_i = a}: (s nu).g = +-Sum[s_ij g_j], and every term of s_ij is a
     derivative along the face of the trace of d_i^m Phi_f; the traces are bicubic and g_j has degree p in each
     variable on every face cell. *)
  hl = Table[elIMat[$elHermite, a, elLagrange[p], 0], {a, 0, 2}]; li = elPInt /@ elLagrange[p];
  Do[
    If[spec["Dirichlet"][[i, side]],
      Module[{tt = Complement[Range[3], {i}], nt, ht, pos, posU, a = dom[[i, side]], xface, gface, terms = {}, acc},
        nt = spec["Cells"][[tt]]; ht = h[[tt]];
        pos = If[side == 1, 1, spec["Cells"][[i]] + 1]; posU = If[side == 1, 1, spec["M"][[i]]];
        xface[arr_, ff_, m_, k1_, k2_] := With[{d = ReplacePart[{0, 0, 0}, {i -> m, tt[[1]] -> Mod[k1, 2], tt[[2]] -> Mod[k2, 2]}]},
          Flatten[el3FaceSlice[arr[[ff, d[[1]] + 2 d[[2]] + 4 d[[3]] + 1]], i, pos][[Quotient[k1, 2] + 1 ;; Quotient[k1, 2] + nt[[1]], Quotient[k2, 2] + 1 ;; Quotient[k2, 2] + nt[[2]]]]]];
        gface[c_, q1_, q2_] := ivPt[Flatten[el3FaceSlice[g[[c]], i, posU][[q1 + 1 ;; q1 + 1 + p (nt[[1]] - 1) ;; p, q2 + 1 ;; q2 + 1 + p (nt[[2]] - 1) ;; p]]]];
        Do[
          Do[
            With[{ff = term[[1]], al = term[[2]], c = term[[3]]},
              Do[
                Module[{comb = None},
                  Do[
                    With[{cf = c ht[[1]]^(1 + Mod[k1, 2] - al[[tt[[1]]]]) ht[[2]]^(1 + Mod[k2, 2] - al[[tt[[2]]]]) hl[[al[[tt[[1]]]] + 1, k1 + 1, q1 + 1]] hl[[al[[tt[[2]]]] + 1, k2 + 1, q2 + 1]]},
                      If[cf != 0, With[{t = ivK[gface[j, q1, q2], cf]}, comb = If[comb === None, t, ndAdd[comb, t]]]]],
                    {q1, 0, p}, {q2, 0, p}];
                  If[comb =!= None, AppendTo[terms, ndMul[{xface[xlo, ff, al[[i]], k1, k2], xface[xhi, ff, al[[i]], k1, k2]}, comb]]]],
                {k1, 0, 3}, {k2, 0, 3}]],
            {term, $el3Beltrami[[$el3Voigt[[i, j]]]]}];
          If[j == i && b[[i]] != 0 && a != 0,
            Do[AppendTo[terms, ivK[gface[i, q1, q2], -b[[i]] a ht[[1]] ht[[2]] li[[q1 + 1]] li[[q2 + 1]]]], {q1, 0, p}, {q2, 0, p}]],
          {j, 3}];
        If[terms =!= {},
          acc = el3TreeSum[{Developer`ToPackedArray[terms[[All, 1]]], Developer`ToPackedArray[terms[[All, 2]]]}];
          bnd = bnd + If[side == 2, {exactTotal[acc[[1]]], exactTotal[acc[[2]]]}, {-exactTotal[acc[[2]]], -exactTotal[acc[[1]]]}]]]],
    {i, 3}, {side, 2}];
  <|"Lower" -> bnd[[1]] - intq[[2]]/2, "Boundary" -> bnd, "IntQ" -> intq|>];

(* ------------------------------------------------------------------------------------------------ *)
(* heuristics: the two solves (machine arithmetic, never trusted)                                    *)
(* ------------------------------------------------------------------------------------------------ *)

(* the finite-element solution in Q_p (any nodal array is a valid candidate) *)
el3SolveCandidate[spec_] := Module[{g1 = Table[elGlobal1D[spec, i], {i, 3}], dm = N[spec["D"]], b = N[spec["b"]], kg, load, u, free, sol},
  kg = ArrayFlatten[Table[
    Total[Join @@ Table[
      If[s1[[2]] == a && s2[[2]] == c && spec["D"][[s1[[1]], s2[[1]]]] != 0,
        dm[[s1[[1]], s2[[1]]]] KroneckerProduct @@ Table[g1[[t]]["GL"][[Boole[s1[[3]] == t] + 1, Boole[s2[[3]] == t] + 1]], {t, 3}],
        Nothing],
      {s1, $el3Strain}, {s2, $el3Strain}]],
    {a, 3}, {c, 3}]];
  load = With[{v = Flatten[Outer[Times, g1[[1]]["VL"], g1[[2]]["VL"], g1[[3]]["VL"]]]}, Join[b[[1]] v, b[[2]] v, b[[3]] v]];
  u = Flatten[spec["G"]];
  free = Flatten[Position[Flatten[{spec["Fixed"], spec["Fixed"], spec["Fixed"]}], 0]];
  If[free =!= {},
    sol = elSolveSparse[kg[[free, free]], (load - kg . u)[[free]]];
    If[sol =!= $Failed, u[[free]] = sol]];
  Developer`ToPackedArray[ArrayReshape[u, Prepend[spec["M"], 3]]]];

(* preconditioned conjugate gradients for a symmetric positive (semi)definite operator *)
el3PCG[op_, rhs_, dg_, tol_, maxit_] := Module[{x = 0. rhs, r = rhs, z, p, rz, rzNew, ap, it = 0, stop = tol Norm[rhs]},
  If[stop == 0., Return[x, Module]];
  z = r/dg; p = z; rz = r . z;
  While[it < maxit && Norm[r] > stop,
    ap = op[p];
    With[{alpha = rz/(p . ap)}, x += alpha p; r -= alpha ap];
    z = r/dg; rzNew = r . z; p = z + (rzNew/rz) p; rz = rzNew; it++];
  x];

(* a maximiser of the lower bound over the free nodal values and the parameters.  The bound is
   beta.x - w.x - x.A.x/2 + const with x = x0 + (free nodal values, scaled by h^d) + Sum[P_k T_k].  The system is
   singular (fields that change the stress functions and not the stress); conjugate gradients with the diagonal as
   preconditioner, after adding 10^-10 times the diagonal. *)
el3SolveStress[spec_] := Module[
  {n1, n2, n3, nn, nc, nd, h = spec["h"], hn = N[spec["h"]], b = spec["b"], p = spec["p"], g = spec["G"], dom = spec["Domain"], sm = spec["S"],
   af, cf, gt, apply, st, pidx, r, x0, tm, free, nfree, sc, wg, beta, hl, lin, tnorm, toX, fromX, rhs, da, diag, eps = 10.^-10, z, x},
  {n1, n2, n3} = spec["Cells"]; nn = spec["Cells"] + 1; nc = n1 n2 n3; nd = 48 Times @@ nn;
  af = N[el3StressA[h, sm]["A"]];
  (* the global index (C order of {6, 8, N1 + 1, N2 + 1, N3 + 1}) of every local unknown on every cell *)
  cf = Flatten[Transpose[Table[
    With[{t = Mod[d[[2]], 2] + 2 Mod[d[[3]], 2] + 4 Mod[d[[4]], 2], e = Quotient[Rest[d], 2]},
      Flatten[Table[((((d[[1]] - 1) 8 + t) nn[[1]] + c1 + e[[1]]) nn[[2]] + c2 + e[[2]]) nn[[3]] + c3 + e[[3]] + 1,
        {c1, 0, n1 - 1}, {c2, 0, n2 - 1}, {c3, 0, n3 - 1}]]],
    {d, $el3Dofs}]]];
  gt = SparseArray[Transpose[{cf, Range[384 nc]}] -> 1., {nd, 384 nc}];
  apply[v_] := gt . Flatten[Partition[v[[cf]], 384] . af];
  st = el3States[spec]; r = Length[st["T"]];
  pidx = If[st["Positions"] === {}, {}, el3LinearIndex[st["Positions"], Join[{6, 8}, nn]]];
  x0 = ConstantArray[0., nd]; If[pidx =!= {}, x0[[pidx]] = Developer`ToPackedArray[N[st["Tb"]]]];
  tm = Table[Module[{v = ConstantArray[0., nd]}, v[[pidx]] = Developer`ToPackedArray[N[st["T"][[k]]]]; v], {k, r}];
  free = Complement[Range[nd], pidx]; nfree = Length[free];
  sc = Flatten[Table[ConstantArray[hn[[1]]^BitGet[t, 0] hn[[2]]^BitGet[t, 1] hn[[3]]^BitGet[t, 2], nn], {f, 6}, {t, 0, 7}]];
  (* the linear term of the cell integrals *)
  wg = ConstantArray[0., nd];
  If[b =!= {0, 0, 0},
    Module[{w = N[el3StressW[h, sm, b]], xcs},
      xcs = N[Flatten[Table[{1, spec["VNodes"][[1, c1]], spec["VNodes"][[2, c2]], spec["VNodes"][[3, c3]]}, {c1, n1}, {c2, n2}, {c3, n3}], 2]];
      wg = gt . Flatten[xcs . Transpose[w]]]];
  (* the boundary term beta.x *)
  beta = ConstantArray[0., nd];
  hl = Table[elIMat[$elHermite, a, elLagrange[p], 0], {a, 0, 2}];
  Do[
    If[spec["Dirichlet"][[i, side]],
      Module[{tt = Complement[Range[3], {i}], nt, ht, pos, posU, sign = If[side == 1, -1., 1.], gf},
        nt = spec["Cells"][[tt]]; ht = h[[tt]];
        pos = If[side == 1, 0, spec["Cells"][[i]]]; posU = If[side == 1, 1, spec["M"][[i]]];
        gf[c_, q1_, q2_] := Flatten[el3FaceSlice[g[[c]], i, posU][[q1 + 1 ;; q1 + 1 + p (nt[[1]] - 1) ;; p, q2 + 1 ;; q2 + 1 + p (nt[[2]] - 1) ;; p]]];
        Do[
          With[{ff = term[[1]], al = term[[2]], c = term[[3]]},
            Do[
              Module[{comb = 0., d, node, idx},
                Do[comb += N[c ht[[1]]^(1 + Mod[k1, 2] - al[[tt[[1]]]]) ht[[2]]^(1 + Mod[k2, 2] - al[[tt[[2]]]]) hl[[al[[tt[[1]]]] + 1, k1 + 1, q1 + 1]] hl[[al[[tt[[2]]]] + 1, k2 + 1, q2 + 1]]] gf[j, q1, q2],
                  {q1, 0, p}, {q2, 0, p}];
                d = ReplacePart[{0, 0, 0}, {i -> al[[i]], tt[[1]] -> Mod[k1, 2], tt[[2]] -> Mod[k2, 2]}];
                idx = Flatten[Table[
                  node = ReplacePart[{0, 0, 0}, {i -> pos, tt[[1]] -> c1 + Quotient[k1, 2], tt[[2]] -> c2 + Quotient[k2, 2]}];
                  ((((ff - 1) 8 + d[[1]] + 2 d[[2]] + 4 d[[3]]) nn[[1]] + node[[1]]) nn[[2]] + node[[2]]) nn[[3]] + node[[3]] + 1,
                  {c1, 0, nt[[1]] - 1}, {c2, 0, nt[[2]] - 1}]];
                If[VectorQ[comb], beta[[idx]] += sign comb]],
              {k1, 0, 3}, {k2, 0, 3}]],
          {j, 3}, {term, $el3Beltrami[[$el3Voigt[[i, j]]]]}]]],
    {i, 3}, {side, 2}];
  lin = beta - wg - apply[x0];
  tnorm = Table[With[{m = Max[Abs[tm[[k]]]]}, If[m > 0., m, 1.]], {k, r}];
  toX[zz_] := Module[{xx = ConstantArray[0., nd]}, xx[[free]] = sc[[free]] zz[[1 ;; nfree]]; If[r > 0, xx += (zz[[nfree + 1 ;;]]/tnorm) . tm]; xx];
  fromX[v_] := Join[sc[[free]] v[[free]], If[r > 0, (tm . v)/tnorm, {}]];
  rhs = fromX[lin];
  da = gt . Flatten[ConstantArray[Diagonal[af], nc]];
  diag = Join[sc[[free]]^2 da[[free]], Table[(tm[[k]]/tnorm[[k]]) . apply[tm[[k]]/tnorm[[k]]], {k, r}]];
  diag = Map[If[# > 0., #, 1.] &, diag];
  z = el3PCG[fromX[apply[toX[#]]] + eps diag # &, rhs, diag, 10.^-9, 20000];
  x = ConstantArray[0., nd]; x[[free]] = sc[[free]] z[[1 ;; nfree]];
  <|"X" -> Developer`ToPackedArray[ArrayReshape[x, Join[{6, 8}, nn]]], "Parameters" -> If[r > 0, z[[nfree + 1 ;;]]/tnorm, {}]|>];

(* ------------------------------------------------------------------------------------------------ *)
(* the competitor, reference data, the certificate                                                   *)
(* ------------------------------------------------------------------------------------------------ *)

el3ImposeDirichlet[spec_, u_] := Developer`ToPackedArray[N[Table[(1 - spec["Fixed"]) u[[a]] + spec["Fixed"] spec["G"][[a]], {a, 3}]]];

(* cand: Automatic (the finite-element solution), a function of the coordinates returning {u1, u2, u3}, or nodal arrays
   of dimensions {3, M1, M2, M3}; refine: replace it by the finite-element solution *)
el3Competitor[spec_, cand_, refine_ : False] := Module[{mm = spec["M"], u},
  If[cand === Automatic, Return[el3SolveCandidate[spec], Module]];
  u = If[ArrayQ[cand, 4] && Dimensions[cand] === Prepend[mm, 3], cand,
    el3Components[Quiet@Apply[cand, Outer[List, spec["UAxes"][[1]], spec["UAxes"][[2]], spec["UAxes"][[3]]], {3}], mm]];
  If[u =!= $Failed, u = ndRealArray[u, Prepend[mm, 3]]];
  If[u === $Failed,
    Throw[failure["BadCandidate", "The candidate must be a function of the coordinates returning three real numbers, or nodal arrays of dimensions {3, M1, M2, M3} with M = Degree Cells + 1."], "CertiField"]];
  If[TrueQ[refine], el3SolveCandidate[spec], el3ImposeDirichlet[spec, u]]];

(* a directory written by export_dir of Reference/Elasticity/elasticity3d.py *)
el3ImportField[spec_, dir_String] := Module[{meta, dims = Join[{6, 8}, spec["Cells"] + 1], x, par},
  meta = Quiet@Import[FileNameJoin[{dir, "field.json"}], "RawJSON"];
  If[! (AssociationQ[meta] && meta["format"] === "CertiField-elasticity3d-field-1" && FileExistsQ[FileNameJoin[{dir, "X.f64"}]]),
    Throw[failure["MissingData", "No stress field of the three-dimensional elasticity reference in " <> dir], "CertiField"]];
  x = BinaryReadList[FileNameJoin[{dir, "X.f64"}], "Real64", ByteOrdering -> -1];
  par = If[FileExistsQ[FileNameJoin[{dir, "par.f64"}]], BinaryReadList[FileNameJoin[{dir, "par.f64"}], "Real64", ByteOrdering -> -1], {}];
  If[Length[x] =!= Times @@ dims, Throw[failure["MissingData", "Wrong size: X.f64"], "CertiField"]];
  <|"X" -> ArrayReshape[x, dims], "Parameters" -> par|>];

el3ImportArray[spec_, dir_String, name_String] := Module[{dims = Prepend[spec["M"], 3], v},
  v = BinaryReadList[FileNameJoin[{dir, name}], "Real64", ByteOrdering -> -1];
  If[Length[v] =!= Times @@ dims, Throw[failure["MissingData", "Wrong size: " <> name], "CertiField"]];
  Developer`ToPackedArray[ArrayReshape[v, dims]]];

certifyLinearElasticity3D[params_Association, cand_, opts_List] := Module[
  {spec, t0 = AbsoluteTime[], verbose, refine, given, u, ub, field, lb, lo, hi, margin},
  spec = linearElasticitySpec3D[KeyDrop[params, "StressField"]];
  verbose = TrueQ[OptionValue[CertifiedVariationalBound, opts, "Verbose"]];
  refine = TrueQ[OptionValue[CertifiedVariationalBound, opts, "RefineCandidate"]];
  If[verbose, Print["  LinearElasticity: ", ndGridString[spec["Cells"]], " cells, degree ", spec["p"], ", ", Length[el3States[spec]["T"]], " parameters on the traction-free faces"]];
  u = el3Competitor[spec, cand, refine];
  ub = el3UpperBound[spec, u];
  If[verbose, Print["  candidate: E <= ", ndNum[ub[[2]]], "   (", ndSec[AbsoluteTime[] - t0], " s)"]];
  given = Lookup[params, "StressField", None];
  field = Which[given === None, el3SolveStress[spec], StringQ[given], el3ImportField[spec, given], True, given];
  field = el3PackField[spec, field];
  If[field === $Failed,
    Throw[failure["BadField", "\"StressField\" must be <|\"X\" -> nodal values of dimensions {6, 8, N1 + 1, N2 + 1, N3 + 1}, \"Parameters\" -> the parameters of the traction-free faces|> for the grid of the model, or a directory written by the Python reference."], "CertiField"]];
  lb = el3LowerBound[spec, field];
  If[verbose, Print["  stress field: E >= ", ndNum[lb["Lower"]], "   (", ndSec[AbsoluteTime[] - t0], " s)"]];
  (* the claims are the computed bounds weakened by one part in 10^9 of the larger modulus, as in two dimensions *)
  lo = lb["Lower"]; hi = ub[[2]];
  margin = Max[Abs[lo], Abs[hi]]/10^9;
  <|"LowerBound" -> lo - margin, "UpperBound" -> hi + margin, "Seconds" -> AbsoluteTime[] - t0, "StressField" -> field,
    "Nodal" -> u, "BoundaryNodalValues" -> spec["G"], "Spec" -> spec|>];

elasticityDescription3D[params_, res_] := Module[{sp = res["Spec"], d = Max[displayDigits[res], 8], nat, grid, datum, free},
  nat = Flatten[Table[If[sp["Dirichlet"][[i, s]], Nothing, {i, If[s == 1, "Lower", "Upper"]}], {i, 3}, {s, 2}], 1];
  grid = ndGridString[sp["Cells"]];
  datum = "the continuous piecewise polynomial of degree " <> ToString[sp["p"]] <> " in each variable through the nodal boundary values on the " <> grid <> " grid";
  free = If[nat === {}, "", "; no traction on the faces " <> StringRiffle[elFaceString[sp["Domain"], #] & /@ nat, ", "]];
  <|"Units" -> "energy units of the model",
    "Theorem" -> Row[{"Every u in H^1 of the box ", exactString[sp["Domain"]], " with values in \[DoubleStruckCapitalR]^3 whose trace on the Dirichlet faces is ", datum,
      " has E[u] \[GreaterEqual] ", lowerString[res["LowerBound"], d], ", where E[u] = \[Integral] (e(u).D.e(u)/2 \[Minus] b.u) dx, e(u) = (\[PartialD]1 u1, \[PartialD]2 u2, \[PartialD]3 u3, \[PartialD]2 u3 + \[PartialD]3 u2, \[PartialD]1 u3 + \[PartialD]3 u1, \[PartialD]1 u2 + \[PartialD]2 u1), D = ",
      exactString[sp["D"]], ", b = ", exactString[sp["b"]], free, "; the competitor has E \[LessEqual] ", upperString[res["UpperBound"], d], "."}],
    "AdmissibleClass" -> Row[{"all u in H^1 of the box ", exactString[sp["Domain"]], " with values in \[DoubleStruckCapitalR]^3 whose trace on the Dirichlet faces is ", datum,
      " (machine numbers, stored in the certificate)", free}],
    "LowerBoundMethod" -> Row[{"complementary-energy bound (Prager\[Dash]Synge) with an exactly equilibrated stress field: Beltrami stress functions in the C^1 tricubic space on the ", grid, " grid, ",
      Length[res["StressField"]["Parameters"]], " parameters on the traction-free faces; exact rational coefficients, outward-rounded arithmetic"}],
    "UpperBoundMethod" -> Row[{"energy of the continuous piecewise polynomial displacement of degree ", sp["p"], " through the candidate's nodal values (exact quadratic form, outward-rounded arithmetic)"}],
    "TrustedBase" -> "IEEE-754 machine arithmetic (+, \[Minus], \[Times]), Wolfram Language exact rational arithmetic, and docs/ELASTICITY_3D.md (Theorem 1, Lemmas 2 and 3, Corollary 4; proofs checked by their author only)",
    "NotClaimed" -> "that the candidate is the minimiser; the gap bounds how far below the candidate any admissible displacement can be, and the energy-norm errors of the candidate and of the stress field",
    "References" -> "docs/ELASTICITY_3D.md; Prager\[Dash]Synge (1947); stress functions of Maxwell, Morera and Beltrami; Fraeijs de Veubeke (1964, 1965), dual analysis with equilibrium models."|>];

(* re-verification: both bounds recomputed from the stored nodal values and the stored stress field, with the Dirichlet
   datum stored in the certificate (it must agree with the model's boundary data to rounding) *)
verifyLinearElasticity3D[a_Association] := Module[{c = a["Certificate"], spec, g, lb, ub},
  spec = linearElasticitySpec3D[KeyDrop[a["Parameters"], "StressField"]];
  g = ndRealArray[Lookup[c, "BoundaryNodalValues", $Failed], Prepend[spec["M"], 3]];
  If[g === $Failed || AnyTrue[Range[3], Max[Abs[(1 - spec["Fixed"]) g[[#]]]] != 0. &] ||
      Max[Abs[g - spec["G"]]] > 10.^-12 (1. + Max[Abs[spec["G"]]]), Return[False, Module]];
  spec = Append[spec, "G" -> g];
  lb = el3LowerBound[spec, Lookup[c, "StressField", $Failed]]["Lower"];
  ub = el3UpperBound[spec, Lookup[c, "Nodal", $Failed]][[2]];
  NumericQ[lb] && NumericQ[ub] && lb >= exactNumber[a["LowerBound"]] && ub <= exactNumber[a["UpperBound"]]];

End[];
