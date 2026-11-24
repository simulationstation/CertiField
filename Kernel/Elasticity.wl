(* ::Package:: *)

(* ::Title:: *)
(*CertiField: linear elasticity on rectangles (plane strain)*)

(* ::Text:: *)
(*Wolfram Language port of Reference/Elasticity/elasticity.py.  Reference/Elasticity/PORTING.md maps the functions onto*)
(*each other and lists the test vectors that every rigorous function below reproduces.*)
(**)
(*Model.  E[u] = Integrate[e(u).D.e(u)/2 - b.u] over a rectangle, e(u) = (d1 u1, d2 u2, d1 u2 + d2 u1), D a constant*)
(*symmetric positive definite 3x3 matrix with rational entries (isotropic material in plane strain:*)
(*D = {{l + 2 m, l, 0}, {l, l + 2 m, 0}, {0, 0, m}}), b a constant rational body force, u = g on the Dirichlet faces, no*)
(*traction on the other faces.  The Dirichlet datum of the certified problem is the trace of the continuous function*)
(*that is a polynomial of degree p (1 or 2) in each variable on every cell and takes the nodal boundary values.*)
(*Admissible class: all u in H^1 with that trace.  Theorem and proofs: docs/ELASTICITY.md.*)
(**)
(*Lower bound (Prager-Synge).  For ANY stress field s = (s11, s22, s12) that is a polynomial on each cell, has a*)
(*continuous normal traction, satisfies div s + b = 0 exactly and has no traction on the traction-free faces,*)
(*    E[u] >= Integrate[(s nu).g over the Dirichlet faces] - 1/2 Integrate[s.S.s],     S = Inverse[D].*)
(*The field is s = (d22 phi, d11 phi, -d12 phi) + (-b1 x1, -b2 x2, 0) with phi in the space of C^1 piecewise bicubic*)
(*functions (nodal values phi, d1 phi, d2 phi, d12 phi).  Any machine numbers at the nodes off the traction-free faces*)
(*give such a field.  At the nodes on the traction-free faces the nodal values are those of a polynomial*)
(*psi_b + Sum[P_k psi_k], found in exact rational arithmetic (elTraces), with machine numbers P_k; this makes the*)
(*traction vanish there exactly.  elLowerBound evaluates the right-hand side: it is a polynomial of degree two in the*)
(*nodal values with exact rational coefficients.*)
(**)
(*Upper bound.  The energy of the Q_p function with the candidate's nodal values: a polynomial of degree two in the*)
(*nodal values with exact rational coefficients (elUpperBound).*)
(**)
(*Rigorous: elUpperBound, elLowerBound, elNodalEnclosure, elTraces and the exact local matrices.  The two linear*)
(*solves (elSolveCandidate, elSolveStress) are ordinary machine arithmetic and are never trusted.*)

Begin["CertiField`Private`"];

(* ------------------------------------------------------------------------------------------------ *)
(* exact one-dimensional tables; polynomials are ascending coefficient lists                          *)
(* ------------------------------------------------------------------------------------------------ *)

(* cubic Hermite basis on [0, 1]: value at 0, slope at 0, value at 1, slope at 1 *)
$elHermite = {{1, 0, -3, 2}, {0, 1, -2, 1}, {0, 0, 3, -2}, {0, 0, -1, 1}};
(* Lagrange basis of degree p on [0, 1] at the nodes j/p *)
elLagrange[1] = {{1, -1}, {0, 1}};
elLagrange[2] = {{1, -3, 2}, {0, 4, -4}, {0, -1, 2}};

elPMul[a_List, b_List] := ListConvolve[a, b, {1, -1}, 0];
elPDer[a_List, k_Integer] := Nest[polyD, a, k];
elPInt[a_List] := a . (1/Range[Length[a]]);                      (* integral over [0, 1] *)
(* exact matrix of Integrate[xi^k A_i^(da)(xi) B_j^(db)(xi), {xi, 0, 1}] *)
elIMat[A_List, da_Integer, B_List, db_Integer, k_Integer : 0] :=
  Table[elPInt[elPMul[UnitVector[k + 1, k + 1], elPMul[elPDer[a, da], elPDer[b, db]]]], {a, A}, {b, B}];
elPow[x_, m_Integer] := If[m == 0, 1, x^m];                      (* 0^0 = 1 *)

$elStrain = {{1, 1, 1}, {2, 2, 2}, {3, 1, 2}, {3, 2, 1}};        (* {m, a, i}: the strain e_m contains d_i u_a *)
$elAiry = {{1, 0, 2, 1}, {2, 2, 0, 1}, {3, 1, 1, -1}};           (* {m, alpha, beta, sign}: s_m = sign d1^alpha d2^beta phi *)

(* exact element matrix of Integrate[e.D.e] over a cell and exact load vector Integrate[b_a N], for the local unknowns
   {a, i1, i2}: component a at the local node (i1, i2) of the Q_p element *)
elStiffnessLocal[p_Integer, h_List, dm_List, b_List] := elStiffnessLocal[p, h, dm, b] = Module[{lag = elLagrange[p], il, dofs, km, intL},
  il = Table[elIMat[lag, da, lag, db], {da, 0, 1}, {db, 0, 1}];
  dofs = Flatten[Table[{a, i1, i2}, {a, 2}, {i1, 0, p}, {i2, 0, p}], 2];
  km = Table[
    Sum[If[s1[[2]] == pp[[1]] && s2[[2]] == qq[[1]],
        dm[[s1[[1]], s2[[1]]]] (h[[1]] h[[2]]/(h[[s1[[3]]]] h[[s2[[3]]]])) *
          il[[Boole[s1[[3]] == 1] + 1, Boole[s2[[3]] == 1] + 1, pp[[2]] + 1, qq[[2]] + 1]] *
          il[[Boole[s1[[3]] == 2] + 1, Boole[s2[[3]] == 2] + 1, pp[[3]] + 1, qq[[3]] + 1]], 0],
      {s1, $elStrain}, {s2, $elStrain}],
    {pp, dofs}, {qq, dofs}];
  intL = elPInt /@ lag;
  If[km =!= Transpose[km], Throw[failure["Internal", "The element matrix is not symmetric."], "CertiField"]];
  <|"Dofs" -> dofs, "K" -> km, "f" -> Table[b[[pp[[1]]]] h[[1]] h[[2]] intL[[pp[[2]] + 1]] intL[[pp[[3]] + 1]], {pp, dofs}]|>];

(* exact element matrix A of Integrate[s.S.s] over a cell for the stress of the 16 bicubic basis functions
   {k1, k2}, k_i = 2 (end of the cell in direction i) + (1 for a derivative in direction i), and the vectors w with
   Integrate[s.S.s_b] = Sum[X_k (w_k1 + xc1 w_k2 + xc2 w_k3)], xc the lower corner of the cell, s_b = (-b1 x1, -b2 x2, 0) *)
elAiryLocal[h_List, sm_List, b_List] := elAiryLocal[h, sm, b] = Module[{ih, jh, dofs, scale, am, w},
  ih = Table[elIMat[$elHermite, a, $elHermite, c], {a, 0, 2}, {c, 0, 2}];
  jh = Table[elPInt[elPMul[UnitVector[k + 1, k + 1], elPDer[hk, a]]], {k, 0, 1}, {a, 0, 2}, {hk, $elHermite}];
  dofs = Flatten[Table[{k1, k2}, {k1, 0, 3}, {k2, 0, 3}], 1];
  scale[{k1_, k2_}] := h[[1]]^Mod[k1, 2] h[[2]]^Mod[k2, 2];
  am = Table[
    Sum[sm[[s1[[1]], s2[[1]]]] s1[[4]] s2[[4]] h[[1]] h[[2]]/(h[[1]]^(s1[[2]] + s2[[2]]) h[[2]]^(s1[[3]] + s2[[3]])) scale[kk] scale[ll] *
        ih[[s1[[2]] + 1, s2[[2]] + 1, kk[[1]] + 1, ll[[1]] + 1]] ih[[s1[[3]] + 1, s2[[3]] + 1, kk[[2]] + 1, ll[[2]] + 1]],
      {s1, $elAiry}, {s2, $elAiry}],
    {kk, dofs}, {ll, dofs}];
  w = Table[
    Sum[With[{c = scale[kk] s1[[4]] h[[1]] h[[2]]/(h[[1]]^s1[[2]] h[[2]]^s1[[3]]), m = s1[[1]], ia = s1[[2]] + 1, ib = s1[[3]] + 1, k1 = kk[[1]] + 1, k2 = kk[[2]] + 1},
        c {(-b[[1]] sm[[m, 1]]) h[[1]] jh[[2, ia, k1]] jh[[1, ib, k2]] + (-b[[2]] sm[[m, 2]]) h[[2]] jh[[1, ia, k1]] jh[[2, ib, k2]],
           (-b[[1]] sm[[m, 1]]) jh[[1, ia, k1]] jh[[1, ib, k2]],
           (-b[[2]] sm[[m, 2]]) jh[[1, ia, k1]] jh[[1, ib, k2]]}],
      {s1, $elAiry}],
    {kk, dofs}];
  If[am =!= Transpose[am], Throw[failure["Internal", "The element matrix is not symmetric."], "CertiField"]];
  <|"Dofs" -> dofs, "A" -> am, "w" -> w|>];

(* exact matrices for the boundary term: Integrate[H_k'' L_q], Integrate[H_k' L_q], Integrate[L_q] *)
elFaceLocal[p_Integer] := elFaceLocal[p] = With[{lag = elLagrange[p]},
  <|"HL2" -> elIMat[$elHermite, 2, lag, 0], "HL1" -> elIMat[$elHermite, 1, lag, 0], "LI" -> elPInt /@ lag|>];

(* ------------------------------------------------------------------------------------------------ *)
(* model normalisation                                                                               *)
(* ------------------------------------------------------------------------------------------------ *)

$elDefaultCells = {32, 32};

(* an array of pairs over the nodes as a pair of arrays over the nodes *)
elComponents[v_, mm_List] := If[ArrayQ[v, 3] && Dimensions[v] === Append[mm, 2], Transpose[v, {2, 3, 1}], $Failed];

linearElasticitySpec[params_Association] := Module[
  {dom, cells, lame, dmat, b, p, nat, dir, h, mm, uaxes, vnodes, ind, fixed, natural, g, gval},
  dom = exactify[Lookup[params, "Domain", Missing[]]];
  If[MatchQ[dom, {{_, _}, {_, _}, {_, _}}],
    Throw[failure["NotImplemented", "LinearElasticity supports rectangles in this version (two dimensions, plane strain); a three-dimensional domain was given."], "CertiField"]];
  If[! (MatchQ[dom, {{_, _}, {_, _}}] && AllTrue[Flatten[dom], MatchQ[#, _Integer | _Rational] &] && AllTrue[dom, #[[1]] < #[[2]] &]),
    Throw[failure["BadModel", "\"Domain\" must be {{a1, b1}, {a2, b2}} with rational a < b."], "CertiField"]];
  lame = Lookup[params, "Lame", Missing[]]; dmat = Lookup[params, "StiffnessMatrix", Missing[]];
  Which[
    ! MissingQ[lame] && ! MissingQ[dmat],
      Throw[failure["BadModel", "Give either \"Lame\" or \"StiffnessMatrix\", not both."], "CertiField"],
    ! MissingQ[lame],
      lame = exactify[lame];
      If[! MatchQ[lame, {_Integer | _Rational, _Integer | _Rational}],
        Throw[failure["BadModel", "\"Lame\" must be {lambda, mu}, two rational numbers."], "CertiField"]];
      dmat = {{lame[[1]] + 2 lame[[2]], lame[[1]], 0}, {lame[[1]], lame[[1]] + 2 lame[[2]], 0}, {0, 0, lame[[2]]}},
    ! MissingQ[dmat], dmat = exactify[dmat],
    True, Throw[failure["BadModel", "The material is missing: give \"Lame\" -> {lambda, mu} or \"StiffnessMatrix\" -> a 3x3 matrix."], "CertiField"]];
  If[! (MatrixQ[dmat, MatchQ[#, _Integer | _Rational] &] && Dimensions[dmat] === {3, 3} && dmat === Transpose[dmat]),
    Throw[failure["BadModel", "\"StiffnessMatrix\" must be a symmetric 3x3 matrix of rational numbers (strains in the order d1 u1, d2 u2, d1 u2 + d2 u1)."], "CertiField"]];
  (* positive definite: the leading principal minors are positive (exact arithmetic) *)
  If[! (dmat[[1, 1]] > 0 && Det[dmat[[1 ;; 2, 1 ;; 2]]] > 0 && Det[dmat] > 0),
    Throw[failure["BadModel", "The stiffness is not positive definite (for \"Lame\" -> {lambda, mu}: mu > 0 and lambda + mu > 0 are needed)."], "CertiField"]];
  b = exactify[Lookup[params, "BodyForce", {0, 0}]];
  If[! MatchQ[b, {_Integer | _Rational, _Integer | _Rational}],
    Throw[failure["BadModel", "\"BodyForce\" must be {b1, b2}, two rational numbers (a constant body force)."], "CertiField"]];
  p = Lookup[params, "Degree", 1];
  If[! MatchQ[p, 1 | 2], Throw[failure["BadModel", "\"Degree\" must be 1 or 2."], "CertiField"]];
  cells = Replace[Lookup[params, "Cells", Automatic], Automatic -> $elDefaultCells];
  If[IntegerQ[cells], cells = {cells, cells}];
  If[! (VectorQ[cells, IntegerQ] && Length[cells] == 2 && Min[cells] >= 1),
    Throw[failure["BadModel", "\"Cells\" must be a positive integer or a list of two positive integers."], "CertiField"]];
  nat = Lookup[params, "TractionFreeFaces", {}];
  If[! MatchQ[nat, {{_Integer, "Lower" | "Upper"} ...}] || ! AllTrue[nat, 1 <= #[[1]] <= 2 &],
    Throw[failure["BadModel", "\"TractionFreeFaces\" must be a list of {direction, \"Lower\" | \"Upper\"}."], "CertiField"]];
  dir = Table[{! MemberQ[nat, {i, "Lower"}], ! MemberQ[nat, {i, "Upper"}]}, {i, 2}];          (* True: Dirichlet *)
  If[! MemberQ[Flatten[dir], True],
    Throw[failure["BadModel", "At least one face must be a Dirichlet face."], "CertiField"]];
  h = (dom[[All, 2]] - dom[[All, 1]])/cells;
  mm = p cells + 1;                                                             (* displacement nodes per direction *)
  uaxes = Table[N[dom[[i, 1]] + j h[[i]]/p], {i, 2}, {j, 0, p cells[[i]]}];
  vnodes = Table[dom[[i, 1]] + j h[[i]], {i, 2}, {j, 0, cells[[i]]}];           (* exact nodes of the Airy grid *)
  (* a displacement node is fixed if it lies on a Dirichlet face; an Airy node is constrained if it lies on a
     traction-free face *)
  ind = Table[ReplacePart[ConstantArray[1, mm[[i]]], {If[dir[[i, 1]], 1 -> 0, Nothing], If[dir[[i, 2]], -1 -> 0, Nothing]}], {i, 2}];
  fixed = 1 - Outer[Times, ind[[1]], ind[[2]]];
  ind = Table[ReplacePart[ConstantArray[1, cells[[i]] + 1], {If[! dir[[i, 1]], 1 -> 0, Nothing], If[! dir[[i, 2]], -1 -> 0, Nothing]}], {i, 2}];
  natural = 1 - Outer[Times, ind[[1]], ind[[2]]];
  (* boundary data: a pair of numbers, a function of the coordinates returning a pair, or nodal arrays {g1, g2} *)
  g = Lookup[params, "BoundaryData", {0, 0}];
  gval = Which[
    MatchQ[g, {_?NumericQ, _?NumericQ}], Table[ConstantArray[g[[a]], mm], {a, 2}],
    ArrayQ[g, 3] && Dimensions[g] === Prepend[mm, 2], g,
    True, elComponents[Quiet@Apply[g, Outer[List, uaxes[[1]], uaxes[[2]]], {2}], mm]];
  (* only the values at the Dirichlet nodes are used: a function may be undefined or singular elsewhere *)
  gval = If[ArrayQ[gval, 3] && Dimensions[gval] === Prepend[mm, 2],
    With[{fx = Flatten[fixed]},
      ndRealArray[Table[ArrayReshape[ReplacePart[ConstantArray[0, Length[fx]], Thread[Flatten[Position[fx, 1]] -> Pick[Flatten[gval[[a]], 1], fx, 1]]], mm], {a, 2}], Prepend[mm, 2]]],
    $Failed];
  If[gval === $Failed,
    Throw[failure["BadModel", "\"BoundaryData\" must be a pair of real numbers, a function of the coordinates returning a pair of real numbers, or nodal arrays of dimensions {2, M1, M2} with M = Degree Cells + 1."], "CertiField"]];
  <|"n" -> 2, "Domain" -> dom, "Cells" -> cells, "h" -> h, "p" -> p, "D" -> dmat, "S" -> Inverse[dmat], "b" -> b, "Dirichlet" -> dir,
    "M" -> mm, "UAxes" -> uaxes, "VNodes" -> vnodes, "Fixed" -> fixed, "Natural" -> natural, "G" -> gval|>];

(* ------------------------------------------------------------------------------------------------ *)
(* the traction-free faces: exact nodal values of the polynomials psi_b, psi_1, ..., psi_r             *)
(* (docs/ELASTICITY.md, Section 3)                                                                    *)
(* ------------------------------------------------------------------------------------------------ *)

(* Returns <|"Nodes" -> the Airy nodes {i, j} on the traction-free faces, "Tb" -> for every such node the four values
   (psi_b, d1 psi_b, d2 psi_b, d12 psi_b), "T" -> the same for psi_1, ..., psi_r|>, all exact rationals.
   psi = Sum[c[4 m + n + 1] x1^m x2^n, m, n = 0..3].  The conditions "second tangential derivative = b_i a" and
   "mixed derivative = 0" on every traction-free face {x_i = a} are linear in c.  The reduced row echelon form gives
   one particular solution (free unknowns zero) and one null vector for each free unknown; of these, the ones are
   kept whose nodal values are independent of each other and of those of the affine functions 1, x1, x2 (adding an
   affine function to the Airy function does not change the stress).  The reduced row echelon form is unique, so
   this is the same family as in the Python reference, and the parameters of a field mean the same in both. *)
elTraces[spec_] := elTracesData[spec["Domain"], spec["Cells"], spec["Dirichlet"], spec["b"]];

elTracesData[dom_, cells_, dir_, b_] := elTracesData[dom, cells, dir, b] = Module[
  {vn, ind, natural, nodes, rows = {}, rhs = {}, a, row, rr, piv, free, part, basis, pattern, kept = {}, push, t = {}, pt},
  vn = Table[dom[[i, 1]] + j (dom[[i, 2]] - dom[[i, 1]])/cells[[i]], {i, 2}, {j, 0, cells[[i]]}];
  ind = Table[ReplacePart[ConstantArray[1, cells[[i]] + 1], {If[! dir[[i, 1]], 1 -> 0, Nothing], If[! dir[[i, 2]], -1 -> 0, Nothing]}], {i, 2}];
  natural = 1 - Outer[Times, ind[[1]], ind[[2]]];
  nodes = Position[natural, 1];
  If[nodes === {}, Return[<|"Nodes" -> {}, "Tb" -> {}, "T" -> {}|>, Module]];
  Do[
    If[! dir[[i, side]],
      a = dom[[i, side]];
      Do[                                                  (* the second tangential derivative: coefficient of x_t^r *)
        row = ConstantArray[0, 16];
        Do[row[[If[i == 1, 4 q + (r + 2), 4 (r + 2) + q] + 1]] = elPow[a, q] (r + 2) (r + 1), {q, 0, 3}];
        AppendTo[rows, row]; AppendTo[rhs, If[r == 0, b[[i]] a, 0]],
        {r, 0, 1}];
      Do[                                                  (* the mixed derivative: coefficient of x_t^r *)
        row = ConstantArray[0, 16];
        Do[row[[If[i == 1, 4 q + (r + 1), 4 (r + 1) + q] + 1]] = q elPow[a, q - 1] (r + 1), {q, 1, 3}];
        AppendTo[rows, row]; AppendTo[rhs, 0],
        {r, 0, 2}]],
    {i, 2}, {side, 2}];
  rr = Select[RowReduce[MapThread[Append, {rows, rhs}]], ! AllTrue[#, # == 0 &] &];
  piv = First[FirstPosition[#, x_ /; x != 0]] & /@ rr;
  If[MemberQ[piv, 17], Throw[failure["Internal", "No polynomial particular solution on the traction-free faces."], "CertiField"]];
  free = Complement[Range[16], piv];
  part = ConstantArray[0, 16];
  Do[part[[piv[[j]]]] = rr[[j, 17]], {j, Length[piv]}];
  basis = Table[Module[{v = UnitVector[16, f]}, Do[v[[piv[[j]]]] = -rr[[j, f]], {j, Length[piv]}]; v], {f, free}];
  (* nodal values (psi, d1 psi, d2 psi, d12 psi) at the nodes of the traction-free faces *)
  pattern[c_] := With[{cm = Partition[c, 4]},
    Table[With[{x = vn[[1, nd[[1]]]], y = vn[[2, nd[[2]]]]},
        With[{px = {1, x, x^2, x^3}, dx = {0, 1, 2 x, 3 x^2}, py = {1, y, y^2, y^3}, dy = {0, 1, 2 y, 3 y^2}},
          {px . cm . py, dx . cm . py, px . cm . dy, dx . cm . dy}]], {nd, nodes}]];
  push[v_] := If[MatrixRank[Append[kept, v]] > Length[kept], AppendTo[kept, v]; True, False];
  Scan[push[Flatten[pattern[UnitVector[16, #]]]] &, {1, 5, 2}];                 (* 1, x1, x2 *)
  Do[pt = pattern[c]; If[push[Flatten[pt]], AppendTo[t, pt]], {c, basis}];
  <|"Nodes" -> nodes, "Tb" -> pattern[part], "T" -> t|>];

(* ------------------------------------------------------------------------------------------------ *)
(* stress fields                                                                                     *)
(* a field is <|"X" -> array of dimensions {4, N1 + 1, N2 + 1}: the nodal values (phi, d1 phi, d2 phi, *)
(* d12 phi) of the Airy function (the entries at the nodes on the traction-free faces are ignored),    *)
(* "Parameters" -> {P_1, ..., P_r}|>                                                                  *)
(* ------------------------------------------------------------------------------------------------ *)

(* the field with its arrays packed, or $Failed if it does not fit the grid *)
elPackField[spec_, f_] := Module[{dims = Prepend[spec["Cells"] + 1, 4], r = Length[elTraces[spec]["T"]], x, par},
  If[! AssociationQ[f], Return[$Failed, Module]];
  x = ndRealArray[Lookup[f, "X", $Failed], dims];
  par = Lookup[f, "Parameters", {}];
  par = If[r == 0, If[par === {}, {}, $Failed], ndRealArray[par, {r}]];
  If[x === $Failed || par === $Failed, $Failed, <|"X" -> x, "Parameters" -> par|>]];

(* RIGOROUS: the nodal values of the Airy function as a pair {lo, hi} of packed arrays.  Nodes off the traction-free
   faces: the stored machine numbers (exact).  Nodes on them: enclosures of the exact rational values of
   psi_b + Sum[P_k psi_k], computed in rational arithmetic from the machine numbers P_k.  With no traction-free face
   the affine part of the Airy function is fixed: phi, d1 phi, d2 phi are zero at the first node. *)
elNodalEnclosure[spec_, f_] := Module[{tr = elTraces[spec], lo = f["X"], hi = f["X"], par, vals, e},
  If[tr["Nodes"] === {},
    Do[lo[[t, 1, 1]] = 0.; hi[[t, 1, 1]] = 0., {t, 3}],
    par = exactNumber /@ Normal[f["Parameters"]];
    vals = tr["Tb"] + Sum[par[[k]] tr["T"][[k]], {k, Length[par]}];
    e = Map[ivC, vals, {2}];
    If[! VectorQ[Flatten[e], MachineNumberQ], ndOverflow[]];
    Do[
      With[{nd = tr["Nodes"][[q]]},
        Do[lo[[t, nd[[1]], nd[[2]]]] = e[[q, t, 1]]; hi[[t, nd[[1]], nd[[2]]]] = e[[q, t, 2]], {t, 4}]],
      {q, Length[tr["Nodes"]]}]];
  ndChk[{lo, hi}]];

(* enclosures of exact numbers as a pair of packed arrays (these are inputs of the interval evaluation, not results
   of machine arithmetic: every entry comes from fLo or fHi) *)
elPackPair[{a_, b_}] := With[{pa = Developer`ToPackedArray[a], pb = Developer`ToPackedArray[b]},
  If[Developer`PackedArrayQ[pa, Real] && Developer`PackedArrayQ[pb, Real], {pa, pb}, ndOverflow[]]];

(* RIGOROUS enclosure, cell by cell, of Sum[m[[i, j]] X_i X_j] for an exact symmetric matrix m and interval arrays X *)
elQuadratic[m_, dofs_, xc_] := Module[{n = Length[dofs], acc = None, term},
  Do[
    If[m[[i, j]] != 0,
      term = If[i == j, ivK[ivSqr[xc[dofs[[i]]]], m[[i, i]]], ivK[ndMul[xc[dofs[[i]]], xc[dofs[[j]]]], 2 m[[i, j]]]];
      acc = If[acc === None, term, ndAdd[acc, term]]],
    {i, n}, {j, i, n}];
  acc];

(* ------------------------------------------------------------------------------------------------ *)
(* the two certified bounds                                                                          *)
(* ------------------------------------------------------------------------------------------------ *)

(* RIGOROUS enclosure {lo, hi} (exact rationals) of the energy of the Q_p function with the nodal values u0
   (dimensions {2, M1, M2}): Sum over the cells of u.K.u/2 - f.u with the exact element matrix and load vector;
   the sum over the cells is exact. *)
elUpperBound[spec_, u0_] := ndQuiet@Module[{p = spec["p"], n1, n2, u, loc, xc, en},
  {n1, n2} = spec["Cells"];
  u = ndRealArray[u0, Prepend[spec["M"], 2]];
  If[u === $Failed, Throw[failure["BadCandidate", "The nodal values must be a real array of dimensions {2, M1, M2} with M = Degree Cells + 1."], "CertiField"]];
  If[! (Max[Abs[spec["Fixed"] (u[[1]] - spec["G"][[1]])]] == 0. && Max[Abs[spec["Fixed"] (u[[2]] - spec["G"][[2]])]] == 0.),
    Throw[failure["Internal", "The Dirichlet nodal values are not imposed exactly."], "CertiField"]];
  loc = elStiffnessLocal[p, spec["h"], spec["D"], spec["b"]];
  xc = Association[Table[d -> ivPt[u[[d[[1]], d[[2]] + 1 ;; d[[2]] + 1 + p (n1 - 1) ;; p, d[[3]] + 1 ;; d[[3]] + 1 + p (n2 - 1) ;; p]]], {d, loc["Dofs"]}]];
  en = ivK[elQuadratic[loc["K"], loc["Dofs"], xc], 1/2];
  Do[If[loc["f"][[j]] != 0, en = ndSub[en, ivK[xc[loc["Dofs"][[j]]], loc["f"][[j]]]]], {j, Length[loc["Dofs"]]}];
  {exactTotal[en[[1]]], exactTotal[en[[2]]]}];

(* RIGOROUS lower bound  Integrate[(s nu).g over the Dirichlet faces] - 1/2 Integrate[s.S.s]  for the stress field of
   the given nodal values and parameters.  Returns <|"Lower" -> the bound, "Boundary" -> enclosure of the first
   term, "IntQ" -> enclosure of Integrate[s.S.s]|> (exact rationals).  The coefficients are exact rationals, the nodal
   values exact machine numbers or enclosures (elNodalEnclosure); sums over cells and cell faces are exact. *)
elLowerBound[spec_, field_] := ndQuiet@Module[
  {f, n1, n2, loc, face, x, xc, q, b = spec["b"], sm = spec["S"], dom = spec["Domain"], h = spec["h"], p = spec["p"], g = spec["G"],
   vn = spec["VNodes"], x1, x2, w, coef, qb, intq, bnd = {0, 0}},
  f = elPackField[spec, field];
  If[f === $Failed,
    Throw[failure["BadField", "The stress field must be <|\"X\" -> nodal values of dimensions {4, N1 + 1, N2 + 1}, \"Parameters\" -> the parameters of the traction-free faces|> for the grid of the model."], "CertiField"]];
  {n1, n2} = spec["Cells"];
  loc = elAiryLocal[h, sm, b]; face = elFaceLocal[p];
  x = elNodalEnclosure[spec, f];
  (* for every cell: the 16 nodal values of its four corners *)
  xc = Association[Table[
    With[{t = Mod[kk[[1]], 2] + 2 Mod[kk[[2]], 2] + 1, s1 = Quotient[kk[[1]], 2] + 1 ;; Quotient[kk[[1]], 2] + n1, s2 = Quotient[kk[[2]], 2] + 1 ;; Quotient[kk[[2]], 2] + n2},
      kk -> {x[[1, t, s1, s2]], x[[2, t, s1, s2]]}], {kk, loc["Dofs"]}]];
  (* the cell integrals: Integrate[s.S.s] = (Airy part) + 2 (cross term with s_b) + (s_b part) *)
  q = elQuadratic[loc["A"], loc["Dofs"], xc];
  If[b =!= {0, 0},
    x1 = With[{e = ivC /@ vn[[1, 1 ;; n1]]}, elPackPair[{Transpose[ConstantArray[e[[All, 1]], n2]], Transpose[ConstantArray[e[[All, 2]], n2]]}]];
    x2 = With[{e = ivC /@ vn[[2, 1 ;; n2]]}, elPackPair[{ConstantArray[e[[All, 1]], n1], ConstantArray[e[[All, 2]], n1]}]];
    Do[
      w = loc["w"][[j]];
      If[w =!= {0, 0, 0},
        coef = ndAdd[ndAdd[ivC[w[[1]]], ndMul[x1, ivC[w[[2]]]]], ndMul[x2, ivC[w[[3]]]]];
        q = ndAdd[q, ivK[ndMul[xc[loc["Dofs"][[j]]], coef], 2]]],
      {j, Length[loc["Dofs"]]}]];
  qb = With[{a1 = dom[[1, 1]], c1 = dom[[1, 2]], a2 = dom[[2, 1]], c2 = dom[[2, 2]]},
    sm[[1, 1]] b[[1]]^2 (c1^3 - a1^3)/3 (c2 - a2) + 2 sm[[1, 2]] b[[1]] b[[2]] (c1^2 - a1^2)/2 (c2^2 - a2^2)/2 +
      sm[[2, 2]] b[[2]]^2 (c1 - a1) (c2^3 - a2^3)/3];
  intq = {exactTotal[q[[1]]], exactTotal[q[[2]]]} + qb;
  (* the boundary term.  On a Dirichlet face {x_i = a}, with t the other variable: (s nu).g =
     +-((d_tt phi - b_i a) g_i - d_t(d_i phi) g_t); phi and d_i phi are cubic on every cell face, g_i and g_t have degree p *)
  Do[
    If[spec["Dirichlet"][[i, side]],
      Module[{j = 3 - i, nt, ht, pos, posU, tt, tn, line, lower, upper, cA, cB, gn, gt, a = dom[[i, side]], acc = None, add, tl, th},
        nt = spec["Cells"][[j]]; ht = h[[j]];
        pos = If[side == 1, 1, spec["Cells"][[i]] + 1]; posU = If[side == 1, 1, spec["M"][[i]]];
        {tt, tn} = If[i == 1, {3, 2}, {2, 3}];           (* where the tangential and the normal derivative are among the four nodal values *)
        line[t_] := If[i == 1, {x[[1, t, pos, All]], x[[2, t, pos, All]]}, {x[[1, t, All, pos]], x[[2, t, All, pos]]}];
        lower[iv_] := {iv[[1, 1 ;; nt]], iv[[2, 1 ;; nt]]};
        upper[iv_] := {iv[[1, 2 ;; nt + 1]], iv[[2, 2 ;; nt + 1]]};
        cA = {lower[line[1]], ivK[lower[line[tt]], ht], upper[line[1]], ivK[upper[line[tt]], ht]};
        cB = {lower[line[tn]], ivK[lower[line[4]], ht], upper[line[tn]], ivK[upper[line[4]], ht]};
        gn = Table[ivPt[If[i == 1, g[[i, posU, qq + 1 ;; qq + 1 + p (nt - 1) ;; p]], g[[i, qq + 1 ;; qq + 1 + p (nt - 1) ;; p, posU]]]], {qq, 0, p}];
        gt = Table[ivPt[If[i == 1, g[[j, posU, qq + 1 ;; qq + 1 + p (nt - 1) ;; p]], g[[j, qq + 1 ;; qq + 1 + p (nt - 1) ;; p, posU]]]], {qq, 0, p}];
        add[term_] := (acc = If[acc === None, term, ndAdd[acc, term]]);
        Do[
          Do[
            If[face["HL2"][[k, qq]] != 0, add[ivK[ndMul[cA[[k]], gn[[qq]]], face["HL2"][[k, qq]]/ht]]];
            If[face["HL1"][[k, qq]] != 0, add[ivK[ndMul[cB[[k]], gt[[qq]]], -face["HL1"][[k, qq]]]]],
            {k, 4}];
          If[b[[i]] != 0 && a != 0, add[ivK[gn[[qq]], -b[[i]] a ht face["LI"][[qq]]]]],
          {qq, p + 1}];
        {tl, th} = {exactTotal[acc[[1]]], exactTotal[acc[[2]]]};
        bnd = bnd + If[side == 2, {tl, th}, {-th, -tl}]]],
    {i, 2}, {side, 2}];
  <|"Lower" -> bnd[[1]] - intq[[2]]/2, "Boundary" -> bnd, "IntQ" -> intq|>];

(* ------------------------------------------------------------------------------------------------ *)
(* heuristics: the two linear solves (machine arithmetic, never trusted)                              *)
(* ------------------------------------------------------------------------------------------------ *)

(* assemble a one-dimensional matrix from the same element matrix on n cells; consecutive cells are shifted by sr
   rows and sc columns *)
elAssemble1D[loc_, n_, sr_, sc_] := Module[{nr, nc, m},
  {nr, nc} = Dimensions[loc];
  m = ConstantArray[0., {sr (n - 1) + nr, sc (n - 1) + nc}];
  Do[m[[sr c + 1 ;; sr c + nr, sc c + 1 ;; sc c + nc]] += loc, {c, 0, n - 1}];
  SparseArray[m]];

(* one-dimensional matrices and vectors of direction i (machine numbers).  Hermite unknowns: (node, value | derivative)
   with the physical derivative as the unknown; Lagrange unknowns: the nodes of the Q_p grid.
   GH[[a + 1, c + 1]] = Integrate[H^(a) H^(c)], VH[[a + 1]] = Integrate[H^(a)], UH[[a + 1]] = Integrate[x H^(a)],
   BHL[[a]] = Integrate[H^(a) L] (a = 1, 2), GL[[a + 1, c + 1]] = Integrate[L^(a) L^(c)], VL = Integrate[L] *)
elGlobal1D[spec_, i_] := Module[{n = spec["Cells"][[i]], h = spec["h"][[i]], p = spec["p"], lag, s, ih, jh, il, xc, vec},
  lag = elLagrange[p]; s = {1, h, 1, h};
  ih = Table[elIMat[$elHermite, a, $elHermite, c], {a, 0, 2}, {c, 0, 2}];
  jh = Table[elPInt[elPMul[UnitVector[k + 1, k + 1], elPDer[hk, a]]], {k, 0, 1}, {a, 0, 2}, {hk, $elHermite}];
  il = Table[elIMat[lag, a, lag, c], {a, 0, 1}, {c, 0, 1}];
  xc = N[spec["VNodes"][[i, 1 ;; n]]];
  vec[len_, stride_, locs_] := Module[{v = ConstantArray[0., len]}, Do[v[[stride c + 1 ;; stride c + Length[locs[[c + 1]]]]] += locs[[c + 1]], {c, 0, n - 1}]; v];
  <|"GH" -> Table[elAssemble1D[N[h^(1 - a - c) Outer[Times, s, s] ih[[a + 1, c + 1]]], n, 2, 2], {a, 0, 2}, {c, 0, 2}],
    "VH" -> Table[vec[2 (n + 1), 2, ConstantArray[N[h^(1 - a) s jh[[1, a + 1]]], n]], {a, 0, 2}],
    "UH" -> Table[vec[2 (n + 1), 2, Table[xc[[c]] N[h^(1 - a) s jh[[1, a + 1]]] + N[h^(2 - a) s jh[[2, a + 1]]], {c, n}]], {a, 0, 2}],
    "BHL" -> Table[elAssemble1D[N[h^(1 - a) s elIMat[$elHermite, a, lag, 0]], n, 2, p], {a, 1, 2}],
    "GL" -> Table[elAssemble1D[N[h^(1 - a - c) il[[a + 1, c + 1]]], n, p, p], {a, 0, 1}, {c, 0, 1}],
    "VL" -> vec[p n + 1, p, ConstantArray[N[h (elPInt /@ lag)], n]]|>];

elSolveSparse[m_, rhs_] := Module[{x = Quiet@Check[LinearSolve[m, rhs], $Failed]},
  If[VectorQ[x, MachineNumberQ], x, $Failed]];

(* heuristic: the finite-element solution in Q_p (any nodal array is a valid candidate) *)
elSolveCandidate[spec_] := Module[{g1 = elGlobal1D[spec, 1], g2 = elGlobal1D[spec, 2], dm = N[spec["D"]], b = N[spec["b"]], mm = spec["M"],
   kg, load, u, free, sol},
  kg = ArrayFlatten[Table[
    Total[Join @@ Table[
      If[s1[[2]] == a && s2[[2]] == c && spec["D"][[s1[[1]], s2[[1]]]] != 0,
        dm[[s1[[1]], s2[[1]]]] KroneckerProduct[g1["GL"][[Boole[s1[[3]] == 1] + 1, Boole[s2[[3]] == 1] + 1]], g2["GL"][[Boole[s1[[3]] == 2] + 1, Boole[s2[[3]] == 2] + 1]]],
        Nothing],
      {s1, $elStrain}, {s2, $elStrain}]],
    {a, 2}, {c, 2}]];
  load = With[{v = Flatten[Outer[Times, g1["VL"], g2["VL"]]]}, Join[b[[1]] v, b[[2]] v]];
  u = Flatten[spec["G"]];
  free = Flatten[Position[Flatten[{spec["Fixed"], spec["Fixed"]}], 0]];
  If[free =!= {},
    sol = elSolveSparse[kg[[free, free]], (load - kg . u)[[free]]];
    If[sol =!= $Failed, u[[free]] = sol]];
  Developer`ToPackedArray[ArrayReshape[u, Prepend[mm, 2]]]];

(* the vector of all Airy unknowns in the order (node 1, derivative 1, node 2, derivative 2) of the Kronecker
   products below, from and to the array {4, N1 + 1, N2 + 1} with first index 1 + d1 + 2 d2 *)
elToKron[x_, n1_, n2_] := Flatten[Transpose[ArrayReshape[x, {2, 2, n1, n2}], {4, 2, 1, 3}]];
elFromKron[v_, n1_, n2_] := ArrayReshape[Transpose[ArrayReshape[v, {n1, 2, n2, 2}], {3, 2, 4, 1}], {4, n1, n2}];

(* heuristic: the maximiser of the lower bound over the free nodal values and the parameters (one sparse solve).
   The bound is beta.x + const - (x.A.x + 2 w.x + const)/2 with x = T z + x0 (free nodal values and parameters z). *)
elSolveStress[spec_] := Module[
  {g1 = elGlobal1D[spec, 1], g2 = elGlobal1D[spec, 2], sm = N[spec["S"]], b = N[spec["b"]], dir = spec["Dirichlet"], g = spec["G"], n1, n2, ntot,
   ag, wg, beta, tr, r, kidx, x0, isfree, free, tm, z, x},
  {n1, n2} = spec["Cells"] + 1; ntot = 4 n1 n2;
  ag = Total[Join @@ Table[
    If[spec["S"][[s1[[1]], s2[[1]]]] != 0,
      (sm[[s1[[1]], s2[[1]]]] s1[[4]] s2[[4]]) KroneckerProduct[g1["GH"][[s1[[2]] + 1, s2[[2]] + 1]], g2["GH"][[s1[[3]] + 1, s2[[3]] + 1]]],
      Nothing],
    {s1, $elAiry}, {s2, $elAiry}]];
  wg = Total[Table[
    s1[[4]] ((-b[[1]] sm[[s1[[1]], 1]]) Flatten[Outer[Times, g1["UH"][[s1[[2]] + 1]], g2["VH"][[s1[[3]] + 1]]]] +
      (-b[[2]] sm[[s1[[1]], 2]]) Flatten[Outer[Times, g1["VH"][[s1[[2]] + 1]], g2["UH"][[s1[[3]] + 1]]]]),
    {s1, $elAiry}]];
  beta = ConstantArray[0., ntot];
  Do[                                                      (* faces {x1 = a}: unknowns (node, 0 | 1, all of direction 2) *)
    If[dir[[1, side]],
      With[{sign = If[side == 1, -1., 1.], pos = If[side == 1, 0, n1 - 1], posU = If[side == 1, 1, spec["M"][[1]]]},
        beta[[(2 pos) (2 n2) + 1 ;; (2 pos + 1) (2 n2)]] += sign (g2["BHL"][[2]] . g[[1, posU, All]]);
        beta[[(2 pos + 1) (2 n2) + 1 ;; (2 pos + 2) (2 n2)]] += -sign (g2["BHL"][[1]] . g[[2, posU, All]])]],
    {side, 2}];
  Do[                                                      (* faces {x2 = a}: unknowns (all of direction 1, node, 0 | 1) *)
    If[dir[[2, side]],
      With[{sign = If[side == 1, -1., 1.], pos = If[side == 1, 0, n2 - 1], posU = If[side == 1, 1, spec["M"][[2]]]},
        beta[[Range[0, 2 n1 - 1] (2 n2) + 2 pos + 1]] += sign (g1["BHL"][[2]] . g[[2, All, posU]]);
        beta[[Range[0, 2 n1 - 1] (2 n2) + 2 pos + 2]] += -sign (g1["BHL"][[1]] . g[[1, All, posU]])]],
    {side, 2}];
  tr = elTraces[spec]; r = Length[tr["T"]];
  kidx[i_, j_, t_] := (2 (i - 1) + Mod[t, 2]) (2 n2) + 2 (j - 1) + Quotient[t, 2] + 1;
  x0 = ConstantArray[0., ntot]; isfree = ConstantArray[1, ntot];
  If[tr["Nodes"] === {},
    Do[isfree[[kidx[1, 1, t]]] = 0, {t, 0, 2}],
    Do[
      With[{nd = tr["Nodes"][[q]]},
        Do[With[{id = kidx[nd[[1]], nd[[2]], t]}, isfree[[id]] = 0; x0[[id]] = N[tr["Tb"][[q, t + 1]]]], {t, 0, 3}]],
      {q, Length[tr["Nodes"]]}]];
  free = Flatten[Position[isfree, 1]];
  tm = SparseArray[Join[
      Thread[Transpose[{free, Range[Length[free]]}] -> 1.],
      Flatten[Table[
        If[tr["T"][[k, q, t + 1]] != 0, {kidx[tr["Nodes"][[q, 1]], tr["Nodes"][[q, 2]], t], Length[free] + k} -> N[tr["T"][[k, q, t + 1]]], Nothing],
        {k, r}, {q, Length[tr["Nodes"]]}, {t, 0, 3}]]],
    {ntot, Length[free] + r}];
  z = elSolveSparse[Transpose[tm] . ag . tm, Transpose[tm] . (beta - wg - ag . x0)];
  If[z === $Failed, z = ConstantArray[0., Length[free] + r]];
  x = ConstantArray[0., ntot];
  x[[free]] = z[[1 ;; Length[free]]];
  <|"X" -> Developer`ToPackedArray[elFromKron[x, n1, n2]], "Parameters" -> z[[Length[free] + 1 ;;]]|>];

(* ------------------------------------------------------------------------------------------------ *)
(* the competitor                                                                                    *)
(* ------------------------------------------------------------------------------------------------ *)

elImposeDirichlet[spec_, u_] := Developer`ToPackedArray[N[Table[(1 - spec["Fixed"]) u[[a]] + spec["Fixed"] spec["G"][[a]], {a, 2}]]];

(* cand: Automatic (the finite-element solution), a function of the coordinates returning {u1, u2}, or nodal arrays
   of dimensions {2, M1, M2}.  refine: replace the given candidate by the finite-element solution (the energy is
   quadratic, so the starting point does not matter). *)
elCompetitor[spec_, cand_, refine_ : False] := Module[{mm = spec["M"], u},
  If[cand === Automatic, Return[elSolveCandidate[spec], Module]];
  u = If[ArrayQ[cand, 3] && Dimensions[cand] === Prepend[mm, 2], cand,
    elComponents[Quiet@Apply[cand, Outer[List, spec["UAxes"][[1]], spec["UAxes"][[2]]], {2}], mm]];
  If[u =!= $Failed, u = ndRealArray[u, Prepend[mm, 2]]];
  If[u === $Failed,
    Throw[failure["BadCandidate", "The candidate must be a function of the coordinates returning a pair of real numbers, or nodal arrays of dimensions {2, M1, M2} with M = Degree Cells + 1."], "CertiField"]];
  If[TrueQ[refine], elSolveCandidate[spec], elImposeDirichlet[spec, u]]];

(* ------------------------------------------------------------------------------------------------ *)
(* fields and certificates written by the Python reference                                            *)
(* ------------------------------------------------------------------------------------------------ *)

(* a directory written by export_dir: field.json and raw little-endian doubles, last index fastest *)
elImportField[spec_, dir_String] := Module[{meta, dims = Prepend[spec["Cells"] + 1, 4], x, par},
  meta = Quiet@Import[FileNameJoin[{dir, "field.json"}], "RawJSON"];
  If[! (AssociationQ[meta] && meta["format"] === "CertiField-elasticity-field-1" && FileExistsQ[FileNameJoin[{dir, "X.f64"}]]),
    Throw[failure["MissingData", "No stress field of the elasticity reference in " <> dir], "CertiField"]];
  x = BinaryReadList[FileNameJoin[{dir, "X.f64"}], "Real64", ByteOrdering -> -1];
  par = If[FileExistsQ[FileNameJoin[{dir, "par.f64"}]], BinaryReadList[FileNameJoin[{dir, "par.f64"}], "Real64", ByteOrdering -> -1], {}];
  If[Length[x] =!= Times @@ dims, Throw[failure["MissingData", "Wrong size: X.f64"], "CertiField"]];
  <|"X" -> ArrayReshape[x, dims], "Parameters" -> par|>];

(* nodal values of a candidate written by the reference (U.f64), and its boundary data (gval.f64) *)
elImportArray[spec_, dir_String, name_String] := Module[{dims = Prepend[spec["M"], 2], v},
  v = BinaryReadList[FileNameJoin[{dir, name}], "Real64", ByteOrdering -> -1];
  If[Length[v] =!= Times @@ dims, Throw[failure["MissingData", "Wrong size: " <> name], "CertiField"]];
  Developer`ToPackedArray[ArrayReshape[v, dims]]];

(* a certificate stored by the reference (save_certificate): stem.json (the problem in exact form, the claimed bounds
   as exact decimals) and stem.npz (boundary data, candidate and stress field as doubles) *)
elImportCertificate[stem_String] := Module[{meta, tmp, files, z, pr, params},
  meta = Quiet@Import[stem <> ".json", "RawJSON"];
  If[! (AssociationQ[meta] && meta["format"] === "CertiField-elasticity-certificate-1" && FileExistsQ[stem <> ".npz"]),
    Throw[failure["MissingData", "No reference certificate at " <> stem], "CertiField"]];
  tmp = CreateDirectory[];
  z = Catch[
    files = ExtractArchive[stem <> ".npz", tmp];
    Association[(FileBaseName[#] -> ndReadNPY[#]) & /@ files], "CertiField"];
  Quiet@DeleteDirectory[tmp, DeleteContents -> True];
  If[! AssociationQ[z], Throw[z, "CertiField"]];
  pr = meta["problem"];
  params = <|"StiffnessMatrix" -> Map[ndRational, pr["D"], {2}], "BodyForce" -> ndRational /@ pr["b"], "Domain" -> Map[ndRational, pr["box"], {2}],
    "Cells" -> pr["N"], "Degree" -> pr["degree"], "BoundaryData" -> z["gval"],
    "TractionFreeFaces" -> Flatten[Table[If[pr["bc"][[i, s]] === "N", {{i, If[s == 1, "Lower", "Upper"]}}, {}], {i, 2}, {s, 2}], 2]|>;
  <|"Parameters" -> params, "Nodal" -> z["U"], "StressField" -> <|"X" -> z["X"], "Parameters" -> Normal[Lookup[z, "par", {}]]|>,
    "ClaimedLower" -> ndRational[meta["lower"]], "ClaimedUpper" -> ndRational[meta["upper"]],
    "Name" -> Lookup[pr, "name", ""], "Info" -> Lookup[meta, "info", <||>]|>];

(* recompute both bounds of a reference certificate; True if they are at least as strong as its claims *)
elVerifyReferenceCertificate[stem_String] := Module[{c = elImportCertificate[stem], spec, lb, ub},
  spec = linearElasticitySpec[c["Parameters"]];
  lb = elLowerBound[spec, c["StressField"]]["Lower"];
  ub = elUpperBound[spec, c["Nodal"]][[2]];
  <|"Verified" -> TrueQ[lb >= c["ClaimedLower"] && ub <= c["ClaimedUpper"] && c["ClaimedLower"] <= c["ClaimedUpper"]],
    "LowerBound" -> lb, "UpperBound" -> ub, "ClaimedLower" -> c["ClaimedLower"], "ClaimedUpper" -> c["ClaimedUpper"]|>];

(* ------------------------------------------------------------------------------------------------ *)
(* the certificate                                                                                   *)
(* ------------------------------------------------------------------------------------------------ *)

certifyLinearElasticity[params_Association, cand_, opts_List] := Module[
  {spec, t0 = AbsoluteTime[], verbose, refine, given, u, ub, field, lb, lo, hi, margin},
  spec = linearElasticitySpec[KeyDrop[params, "StressField"]];
  verbose = TrueQ[OptionValue[CertifiedVariationalBound, opts, "Verbose"]];
  refine = TrueQ[OptionValue[CertifiedVariationalBound, opts, "RefineCandidate"]];
  If[verbose, Print["  LinearElasticity: ", ndGridString[spec["Cells"]], " cells, degree ", spec["p"], ", ", Length[elTraces[spec]["T"]], " parameters on the traction-free faces"]];
  u = elCompetitor[spec, cand, refine];
  ub = elUpperBound[spec, u];
  If[verbose, Print["  candidate: E <= ", ndNum[ub[[2]]], "   (", ndSec[AbsoluteTime[] - t0], " s)"]];
  given = Lookup[params, "StressField", None];
  field = Which[given === None, elSolveStress[spec], StringQ[given], elImportField[spec, given], True, given];
  field = elPackField[spec, field];
  If[field === $Failed,
    Throw[failure["BadField", "\"StressField\" must be <|\"X\" -> nodal values of dimensions {4, N1 + 1, N2 + 1}, \"Parameters\" -> the parameters of the traction-free faces|> for the grid of the model, or a directory written by the Python reference."], "CertiField"]];
  lb = elLowerBound[spec, field];
  If[verbose, Print["  stress field: E >= ", ndNum[lb["Lower"]], "   (", ndSec[AbsoluteTime[] - t0], " s)"]];
  (* the claims are the computed bounds weakened by one part in 10^9 of the larger modulus, as in the reference
     (claim_strings) and in ScalarFieldND *)
  lo = lb["Lower"]; hi = ub[[2]];
  margin = Max[Abs[lo], Abs[hi]]/10^9;
  <|"LowerBound" -> lo - margin, "UpperBound" -> hi + margin, "Seconds" -> AbsoluteTime[] - t0, "StressField" -> field,
    "Nodal" -> u, "BoundaryNodalValues" -> spec["G"], "Spec" -> spec|>];

elFaceString[dom_, {i_, side_}] := "x" <> ToString[i] <> " = " <> exactString[dom[[i, If[side === "Lower", 1, 2]]]];

elasticityDescription[params_, res_] := Module[{sp = res["Spec"], d = Max[displayDigits[res], 8], nat, grid, datum, free},
  nat = Flatten[Table[If[sp["Dirichlet"][[i, s]], Nothing, {i, If[s == 1, "Lower", "Upper"]}], {i, 2}, {s, 2}], 1];
  grid = ndGridString[sp["Cells"]];
  datum = "the continuous function that is a polynomial of degree " <> ToString[sp["p"]] <> " in each variable on every cell of the " <> grid <>
    " grid and takes the nodal boundary values";
  free = If[nat === {}, "", "; no traction on the faces " <> StringRiffle[elFaceString[sp["Domain"], #] & /@ nat, ", "]];
  <|"Units" -> "energy units of the model",
    "Theorem" -> Row[{"Every u in H^1 of the rectangle ", exactString[sp["Domain"]], " with values in \[DoubleStruckCapitalR]^2 whose trace on the Dirichlet faces is ", datum,
      " has E[u] \[GreaterEqual] ", lowerString[res["LowerBound"], d], ", where E[u] = \[Integral] (e(u).D.e(u)/2 \[Minus] b.u) dx, e(u) = (\[PartialD]1 u1, \[PartialD]2 u2, \[PartialD]1 u2 + \[PartialD]2 u1), D = ",
      exactString[sp["D"]], ", b = ", exactString[sp["b"]], free, "; the competitor has E \[LessEqual] ", upperString[res["UpperBound"], d], "."}],
    "AdmissibleClass" -> Row[{"all u in H^1 of the rectangle ", exactString[sp["Domain"]], " with values in \[DoubleStruckCapitalR]^2 whose trace on the Dirichlet faces is ", datum,
      " (machine numbers, stored in the certificate)", free}],
    "LowerBoundMethod" -> Row[{"complementary-energy bound (Prager\[Dash]Synge) with an exactly equilibrated stress field: Airy function in the C^1 bicubic space on the ", grid, " grid, ",
      Length[res["StressField"]["Parameters"]], " parameters on the traction-free faces; exact rational coefficients, outward-rounded arithmetic"}],
    "UpperBoundMethod" -> Row[{"energy of the continuous piecewise polynomial displacement of degree ", sp["p"], " through the candidate's nodal values (exact quadratic form, outward-rounded arithmetic)"}],
    "TrustedBase" -> "IEEE-754 machine arithmetic (+, \[Minus], \[Times]), Wolfram Language exact rational arithmetic, and docs/ELASTICITY.md (Theorem 1, Lemmas 2 and 3; proofs checked by their author only)",
    "NotClaimed" -> "that the candidate is the minimiser; the gap bounds how far below the candidate any admissible displacement can be, and the energy-norm errors of the candidate and of the stress field",
    "References" -> "docs/ELASTICITY.md; Prager\[Dash]Synge (1947); Fraeijs de Veubeke (1964, 1965), dual analysis with equilibrium models."|>];

(* re-verification: both bounds recomputed from the stored nodal values and the stored stress field.  The Dirichlet
   datum is the one stored in the certificate (exact machine numbers); it must agree with the boundary data of the
   model to rounding, because a function of the coordinates need not evaluate to the same last bits on another machine. *)
verifyLinearElasticity[a_Association] := Module[{c = a["Certificate"], spec, g, lb, ub},
  spec = linearElasticitySpec[KeyDrop[a["Parameters"], "StressField"]];
  g = ndRealArray[Lookup[c, "BoundaryNodalValues", $Failed], Prepend[spec["M"], 2]];
  If[g === $Failed || Max[Abs[(1 - spec["Fixed"]) g[[1]]]] != 0. || Max[Abs[(1 - spec["Fixed"]) g[[2]]]] != 0. ||
      Max[Abs[g - spec["G"]]] > 10.^-12 (1. + Max[Abs[spec["G"]]]), Return[False, Module]];
  spec = Append[spec, "G" -> g];
  lb = elLowerBound[spec, Lookup[c, "StressField", $Failed]]["Lower"];
  ub = elUpperBound[spec, Lookup[c, "Nodal", $Failed]][[2]];
  NumericQ[lb] && NumericQ[ub] && lb >= exactNumber[a["LowerBound"]] && ub <= exactNumber[a["UpperBound"]]];

CertiFieldModelInfo["LinearElasticity"] := <|
  "Name" -> "LinearElasticity",
  "Energy" -> "E[u] = \[Integral] (e(u).D.e(u)/2 \[Minus] b.u) dx over a rectangle (plane strain), e(u) = (\[PartialD]1 u1, \[PartialD]2 u2, \[PartialD]1 u2 + \[PartialD]2 u1)",
  "Parameters" -> <|"Lame" -> "{lambda, mu}, rational: isotropic material, D = {{lambda + 2 mu, lambda, 0}, {lambda, lambda + 2 mu, 0}, {0, 0, mu}}",
    "StiffnessMatrix" -> "instead of \"Lame\": a symmetric positive definite 3x3 matrix D with rational entries",
    "BodyForce" -> "{b1, b2}, rational, constant (default {0, 0})",
    "Domain" -> "{{a1, b1}, {a2, b2}}, rational",
    "Cells" -> "grid, an integer or one integer per direction (default 32)",
    "Degree" -> "1 or 2: the degree p of the displacement and of the datum in each variable (default 1)",
    "BoundaryData" -> "a pair of numbers, a function of the coordinates returning a pair, or nodal arrays of dimensions {2, p N1 + 1, p N2 + 1}; the certified datum is the piecewise polynomial of degree p through the nodal values on the boundary grid",
    "TractionFreeFaces" -> "list of {direction, \"Lower\" | \"Upper\"}: faces with no traction (default: none; at least one face must remain a Dirichlet face)",
    "StressField" -> "optional: a stress field to certify, <|\"X\" -> nodal values of the Airy function, \"Parameters\" -> ...|> or a directory written by the Python reference; it replaces the package's own choice"|>,
  "Options" -> <|"RefineCandidate" -> "True: the given candidate is replaced by the finite-element solution on the grid (default False: it is certified as given)",
    "Verbose" -> "True prints the two bounds as they are computed"|>,
  "AdmissibleClass" -> "all u in H^1 of the rectangle with values in \[DoubleStruckCapitalR]^2 and the given trace on the Dirichlet faces",
  "LowerBound" -> "E \[GreaterEqual] \[Integral] (\[Sigma]\[Nu]).g over the Dirichlet faces \[Minus] 1/2 \[Integral] \[Sigma].S.\[Sigma] (Prager\[Dash]Synge) for any exactly equilibrated stress field \[Sigma] with no traction on the traction-free faces; \[Sigma] comes from an Airy function in the C^1 bicubic space, with exact nodal values on the traction-free faces; evaluated with exact coefficients and outward-rounded arithmetic",
  "UpperBound" -> "energy of the continuous piecewise polynomial displacement (degree p) through the candidate's nodal values (rigorous enclosure)",
  "Reference" -> "docs/ELASTICITY.md",
  "Status" -> "Supported"|>;

End[];
