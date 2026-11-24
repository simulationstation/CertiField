(* ::Package:: *)

(* ::Title:: *)
(*CertiField: scalar fields on boxes in two and three dimensions*)

(* ::Text:: *)
(*Wolfram Language port of Reference/ScalarND/scalarnd.py.  Reference/ScalarND/PORTING.md maps the functions onto*)
(*each other, lists the test vectors every rigorous function below reproduces, and logs what was changed when the*)
(*port was first executed.*)
(**)
(*Model.  E[u] = Integrate[k |grad u|^2 + W(u)] over a box, u = g on the Dirichlet faces, free on the natural faces.*)
(*W is a polynomial with rational coefficients (it need not be convex).  The Dirichlet datum of the certified*)
(*problem is the piecewise multilinear interpolant of the nodal values of g on the boundary grid.*)
(*Admissible class: all u in H^1 with that trace.  Theorem and proofs: docs/SCALAR_ND.md.*)
(**)
(*Lower bound.  For ANY field S(x, t) = t sigma0(x) + Sum[a_m(x) b_m(clamp(t))] (lowest-order Raviart-Thomas in x,*)
(*clamped C^1 quadratic B-splines in t),*)
(*    E[u] >= Integrate[S(x, g) . nu over the Dirichlet faces] - Integrate[sup_t Phi(x, t)],*)
(*    Phi = div_x S + |d_t S|^2/(4k) - W(t).*)
(*ndLowerBound evaluates the right-hand side rigorously for a given field.  Two ways of choosing the field are*)
(*provided, both heuristic (they affect tightness, never validity):*)
(*  convex mode       the affine part alone, from a cell-centred finite-volume solve (ndConvexField); tight for convex W;*)
(*  calibration mode  the spline part alone, optimised by L-BFGS on a smoothed copy of the certified functional,*)
(*                    coarse to fine (ndCalibrate); needed when W is not convex.*)
(*A field computed elsewhere (for instance by the Python reference) can be imported with ndImportField.*)
(**)
(*Upper bound.  The energy of the Q1 (piecewise multilinear) function with the candidate's nodal values, enclosed*)
(*rigorously: exact quadratic form for the gradient term, midpoint rule with remainder for the potential term.*)
(**)
(*Rigorous: ndUpperBound, ndFluxLower, ndLowerBound, conjUpper, conjCells, tailInf, splineEval.  Everything else*)
(*(candidates, finite-volume solve, optimiser, embeddings, placement of tables) is ordinary machine arithmetic and*)
(*is never trusted.*)

Begin["CertiField`Private`"];

(* ------------------------------------------------------------------------------------------------ *)
(* checked machine arithmetic on {lo, hi} pairs (machine numbers or packed arrays)                    *)
(* ------------------------------------------------------------------------------------------------ *)

(* Overflow.  The outward rounding of RigorousArithmetic.wl encloses the exact result of an IEEE-754 operation.  When
   a machine operation overflows, the kernel does not return an infinity: it silently repeats the operation in
   software arithmetic and returns a number that is not a machine number, and an array holding one is no longer a
   packed array.  The enclosure argument does not cover that case, so every result of the helpers below is tested
   (a flag test, of no cost) and the evaluation stops as "Inconclusive" if an overflow has happened.  Two rules keep
   the test meaningful: inputs are packed before they are used, and Developer`ToPackedArray is never applied to the
   result of an arithmetic operation.
   Underflow.  Exact zeros are common here (zero boundary data, zero flux on natural faces), and outward rounding
   turns them into numbers of the order of 10^-300 whose products with 2^-49 are subnormal.  IEEE gradual underflow
   rounds those correctly to within 2^-1075, far inside the absolute widening 10^-300 of rUp and rDn, so the
   enclosures stay valid; only the kernel's warning about it (issued for unpacked numbers) is suppressed. *)
ndMachQ[x_] := MachineNumberQ[x] || Developer`PackedArrayQ[x, Real];
ndOverflow[] := Throw[failure["Inconclusive", "A machine number overflowed in the certified evaluation: the data are too large for machine arithmetic."], "CertiField"];
ndChk[iv : {a_, b_}] := If[ndMachQ[a] && ndMachQ[b], iv, ndOverflow[]];
ndChk[_] := ndOverflow[];
ndChk1[x_] := If[ndMachQ[x], x, ndOverflow[]];
SetAttributes[ndQuiet, HoldFirst];
ndQuiet[expr_] := Quiet[expr, {General::munfl}];

ndAdd[x_, y_] := ndChk[ivAdd[x, y]];
ndSub[x_, y_] := ndChk[ivSub[x, y]];
ndMul[x_, y_] := ndChk[ivMul[x, y]];
ndDiv[x_, y_] := ndChk[ivDivPos[x, y]];                    (* division by an interval with positive lower end *)
ndUp[x_] := ndChk1[rUp[x]];
ndDn[x_] := ndChk1[rDn[x]];

ivC[q_] := {fLo[q], fHi[q]};                              (* enclosure of an exact number *)
ivK[iv_, q_] := ndMul[iv, ivC[q]];                        (* times an exact number *)
ivSqr[{a_, b_}] := Module[{l2 = a a, h2 = b b, straddle = UnitStep[-a] UnitStep[b]},
  ndChk[{Ramp[rDn[(1 - straddle) vMin[l2, h2]]], rUp[vMax[l2, h2]]}]];
ivPos[{a_, b_}] := {Ramp[a], Ramp[b]};                    (* enclosure of Max[x, 0] *)
ivSum[list_List] := Fold[ndAdd, list];
(* polynomial with exact ascending coefficients at an interval *)
ivHornerQ[coeffs_List, t_] := Fold[ndAdd[ndMul[#1, t], ivC[#2]] &, ivC[Last[coeffs]], Reverse[Most[coeffs]]];
(* exact rational interval Horner: {lo, hi} rationals containing the range of the polynomial on [a, b] *)
exactRange[coeffs_List, a_, b_] := Module[{iv = Fold[#1 Interval[{a, b}] + #2 &, Interval[{Last[coeffs], Last[coeffs]}], Reverse[Most[coeffs]]]},
  {Min[iv], Max[iv]}];
polyAt[coeffs_List, t_] := Fold[#1 t + #2 &, 0, Reverse[coeffs]];
polyD[coeffs_List] := If[Length[coeffs] <= 1, {0}, Rest[coeffs] Range[Length[coeffs] - 1]];
(* exact sum of machine numbers (rational arithmetic, no rounding) *)
exactTotal[x_] := Total[SetPrecision[Flatten[{x}], Infinity]];

(* a real array of the given dimensions as a packed array; $Failed for anything else (exact numbers are converted,
   complex numbers, symbols and numbers outside the machine range are refused) *)
ndRealArray[x_, dims_List] := Module[{p = Quiet@Developer`ToPackedArray[N[x]]},
  If[Developer`PackedArrayQ[p, Real] && Dimensions[p] === dims, p, $Failed]];

(* ------------------------------------------------------------------------------------------------ *)
(* array slicing; the first n dimensions of every array are spatial                                   *)
(* ------------------------------------------------------------------------------------------------ *)

sliceAt[a_, n_, i_, sp_] := a[[Sequence @@ ReplacePart[ConstantArray[All, n], i -> sp]]];
loFace[a_, n_, i_] := sliceAt[a, n, i, Span[1, -2]];              (* for every cell: its lower face in direction i *)
hiFace[a_, n_, i_] := sliceAt[a, n, i, Span[2, All]];               (* for every cell: its upper face in direction i *)
ivLo[iv_, n_, i_] := {loFace[iv[[1]], n, i], loFace[iv[[2]], n, i]};
ivHi[iv_, n_, i_] := {hiFace[iv[[1]], n, i], hiFace[iv[[2]], n, i]};
(* for every cell: its vertex v in {0, 1}^n, from a node array *)
vertexPart[u_, n_, v_List] := u[[Sequence @@ Table[If[v[[j]] == 0, Span[1, -2], Span[2, All]], {j, n}]]];
lastPart[a_, n_, m_] := a[[Sequence @@ ConstantArray[All, n], m]];
(* a slab (an array over the other n - 1 directions) as an array with a singleton dimension i *)
asSlab[a_, i_] := ArrayReshape[a, Insert[Dimensions[a], 1, i]];

(* ------------------------------------------------------------------------------------------------ *)
(* model normalisation                                                                               *)
(* ------------------------------------------------------------------------------------------------ *)

$ndDefaultCells = <|2 -> {32, 32}, 3 -> {16, 16, 16}|>;

scalarFieldNDSpec[params_Association] := Module[
  {V, k, dom, n, cells, nat, g, s, vexpr, W, h, axes, dir, ind, fixed, gval, coords},
  V = Lookup[params, "Potential", Missing[]];
  k = exactify[Lookup[params, "Stiffness", 1/2]];
  dom = exactify[Lookup[params, "Domain", Missing[]]];
  If[! MatchQ[V, _Function | _Symbol], Throw[failure["BadModel", "\"Potential\" must be a pure function of one variable."], "CertiField"]];
  If[! (MatchQ[dom, {{_, _} ..}] && 2 <= Length[dom] <= 3 && AllTrue[Flatten[dom], MatchQ[#, _Integer | _Rational] &] && AllTrue[dom, #[[1]] < #[[2]] &]),
    Throw[failure["BadModel", "\"Domain\" must be {{a1, b1}, {a2, b2}} or {{a1, b1}, {a2, b2}, {a3, b3}} with rational a < b (use ScalarField1D in one dimension)."], "CertiField"]];
  If[! (MatchQ[k, _Integer | _Rational] && k > 0), Throw[failure["BadModel", "\"Stiffness\" must be a positive rational number."], "CertiField"]];
  n = Length[dom];
  cells = Replace[Lookup[params, "Cells", Automatic], Automatic -> $ndDefaultCells[n]];
  If[IntegerQ[cells], cells = ConstantArray[cells, n]];
  If[! (VectorQ[cells, IntegerQ] && Length[cells] == n && Min[cells] >= 1),
    Throw[failure["BadModel", "\"Cells\" must be a positive integer or a list of n positive integers."], "CertiField"]];
  nat = Lookup[params, "NaturalFaces", {}];
  If[! MatchQ[nat, {{_Integer, "Lower" | "Upper"} ...}] || ! AllTrue[nat, 1 <= #[[1]] <= n &],
    Throw[failure["BadModel", "\"NaturalFaces\" must be a list of {direction, \"Lower\" | \"Upper\"}."], "CertiField"]];
  vexpr = Expand[exactify[V[s]]];
  If[! PolynomialQ[vexpr, s],
    Throw[failure["NotImplemented", "ScalarFieldND supports polynomial potentials in this version."], "CertiField"]];
  W = CoefficientList[vexpr, s];
  If[W === {} , W = {0}];
  If[! VectorQ[W, MatchQ[#, _Integer | _Rational] &],
    Throw[failure["BadModel", "The potential must be a polynomial with rational coefficients and no free parameters."], "CertiField"]];
  h = (dom[[All, 2]] - dom[[All, 1]])/cells;
  axes = Table[N[dom[[i, 1]] + j h[[i]]], {i, n}, {j, 0, cells[[i]]}];
  dir = Table[{! MemberQ[nat, {i, "Lower"}], ! MemberQ[nat, {i, "Upper"}]}, {i, n}];          (* True: Dirichlet *)
  (* a node is fixed if it lies on a Dirichlet face *)
  ind = Table[ReplacePart[ConstantArray[1, cells[[i]] + 1], {If[dir[[i, 1]], 1 -> 0, Nothing], If[dir[[i, 2]], -1 -> 0, Nothing]}], {i, n}];
  fixed = 1 - Outer[Times, Sequence @@ ind];
  (* boundary data: a number, a function of the coordinates, or an array of nodal values on the grid *)
  g = Lookup[params, "BoundaryData", 0];
  gval = Which[
    NumericQ[g], ConstantArray[N[g], cells + 1],
    ArrayQ[g, n] && Dimensions[g] === cells + 1, g,
    True, coords = Outer[List, Sequence @@ axes]; Quiet@Apply[g, coords, {n}]];
  (* only the values at the Dirichlet nodes are used: a function may be undefined or singular elsewhere *)
  gval = If[ArrayQ[gval, n] && Dimensions[gval] === cells + 1,
    With[{fx = Flatten[fixed], flat = Flatten[gval, n - 1]},
      ndRealArray[ArrayReshape[ReplacePart[ConstantArray[0, Length[fx]], Thread[Flatten[Position[fx, 1]] -> Pick[flat, fx, 1]]], cells + 1], cells + 1]],
    $Failed];
  If[gval === $Failed,
    Throw[failure["BadModel", "\"BoundaryData\" must be a real number, a real function of the coordinates, or an array of nodal values of dimensions Cells + 1."], "CertiField"]];
  <|"n" -> n, "k" -> k, "W" -> W, "W1" -> polyD[W], "W2" -> polyD[polyD[W]], "Domain" -> dom, "Cells" -> cells, "h" -> h,
    "Axes" -> axes, "Dirichlet" -> dir, "Fixed" -> fixed, "G" -> Developer`ToPackedArray[N[fixed gval]],
    "CellVolume" -> Times @@ h, "Vexpr" -> vexpr, "s" -> s, "Superlinear" -> superlinearQ[W], "Convex" -> convexPolyQ[W]|>];

(* even degree >= 2 with positive leading coefficient *)
superlinearQ[W_List] := With[{d = Length[W] - 1}, d >= 2 && EvenQ[d] && Last[W] > 0];
(* W'' >= 0 on the whole line, decided exactly; only used to choose between the two modes *)
convexPolyQ[W_List] := Module[{t}, With[{e = polyAt[polyD[polyD[W]], t]}, TrueQ[Quiet@Resolve[ForAll[t, e >= 0], Reals]]]];

(* the same model on another grid (the coarse levels of calibration mode): an array of nodal boundary values is
   restricted to the coarser nodes *)
ndParamsOnGrid[params_, cells_List, cellsFine_List] := Module[{g = Lookup[params, "BoundaryData", 0], p = Append[params, "Cells" -> cells]},
  If[ArrayQ[g, Length[cells]] && Dimensions[g] === cellsFine + 1 && cells =!= cellsFine,
    p["BoundaryData"] = g[[Sequence @@ Table[Span[1, All, cellsFine[[j]]/cells[[j]]], {j, Length[cells]}]]]];
  p];

(* ------------------------------------------------------------------------------------------------ *)
(* the competitor: Q1 function through the candidate's nodal values, Dirichlet values imposed exactly  *)
(* ------------------------------------------------------------------------------------------------ *)

(* cand: Automatic (the package computes one), a function of the coordinates, or an array of nodal values.
   refine: use the given candidate as the starting point of the package's own minimisation instead of as it is. *)
ndCompetitor[spec_, cand_, refine_ : False, gfun_ : None] := Module[{n = spec["n"], dims = spec["Cells"] + 1, u, start = $Failed},
  If[cand === Automatic,
    (* start from the boundary function evaluated at all nodes when W is not convex and there is such a function,
       otherwise from the finite-volume solution (for convex W the minimiser is unique) *)
    If[! TrueQ[spec["Convex"]] && gfun =!= None,
      start = ndRealArray[Quiet@Apply[gfun, Outer[List, Sequence @@ spec["Axes"]], {n}], dims]];
    If[start === $Failed, start = ndNodalFromCells[spec, ndFVSolve[spec, ConstantArray[0., spec["Cells"]]]]];
    (* when W is not convex a symmetric start can sit on a symmetric saddle point of the energy: a small fixed
       perturbation lets the descent leave it *)
    If[! TrueQ[spec["Convex"]],
      start = start + 10.^-3 (1. + Max[Abs[start]]) BlockRandom[SeedRandom[20261123]; RandomReal[{-1, 1}, dims]]];
    Return[ndImposeDirichlet[spec, ndSolveCandidate[spec, ndImposeDirichlet[spec, start]]]]];
  u = If[ArrayQ[cand, n] && Dimensions[cand] === dims, cand, Quiet@Apply[cand, Outer[List, Sequence @@ spec["Axes"]], {n}]];
  u = ndRealArray[u, dims];
  If[u === $Failed,
    Throw[failure["BadCandidate", "The candidate must be a real function of the coordinates or an array of nodal values of dimensions Cells + 1."], "CertiField"]];
  u = ndImposeDirichlet[spec, u];
  If[TrueQ[refine], ndImposeDirichlet[spec, ndSolveCandidate[spec, u]], u]];

ndImposeDirichlet[spec_, u_] := Developer`ToPackedArray[N[(1 - spec["Fixed"]) u + spec["Fixed"] spec["G"]]];

(* nodal values from cell values: mean of the adjacent cells (heuristic; any nodal array is a valid candidate) *)
ndNodalFromCells[spec_, p_] := Module[{n = spec["n"], pad, cnt, tot},
  pad = ArrayPad[p, 1, 0.];
  cnt = ArrayPad[ConstantArray[1., Dimensions[p]], 1, 0.];
  tot = Total[vertexPart[pad, n, #] & /@ Tuples[{0, 1}, n]];
  tot/Total[vertexPart[cnt, n, #] & /@ Tuples[{0, 1}, n]]];

(* heuristic linear solves: conjugate gradients with an incomplete factorisation for large systems (the matrices are
   symmetric positive definite on the normal path of both Newton iterations below), a direct solve for small ones
   or if that fails *)
ndSolveSparse[m_, rhs_] := Module[{x = $Failed},
  If[Length[rhs] > 5000,
    x = Quiet@Check[LinearSolve[m, rhs, Method -> {"Krylov", Method -> "ConjugateGradient", "Preconditioner" -> "ILU0", Tolerance -> 10.^-12}], $Failed];
    If[! (VectorQ[x, MachineNumberQ] && Norm[m . x - rhs] <= 10.^-8 Norm[rhs]), x = $Failed]];
  If[x === $Failed, x = Quiet@Check[LinearSolve[m, rhs], $Failed]];
  If[VectorQ[x, MachineNumberQ], x, $Failed]];

(* heuristic: conjugate gradients for H d = -g with the diagonal as preconditioner, stopped at the first direction
   of non-positive curvature (Steihaug).  The result is a descent direction whether or not H is positive definite,
   so the Newton iteration that uses it decreases the energy and cannot settle on a saddle point. *)
ndTruncatedCG[h_, g_, tol_ : 10.^-10, maxit_ : 1000] := Module[{dg, d = 0. g, r = g, z, p, rz, hp, curv, al, rzn, g0 = Norm[g]},
  dg = Abs[Normal[Diagonal[h]]]; dg = dg + 10.^-12 Max[dg] + 10.^-300;
  z = r/dg; p = -z; rz = r . z;
  Do[
    hp = h . p; curv = p . hp;
    If[! (curv > 10.^-12 (p . p) Max[dg]), If[j == 1, d = p]; Break[]];
    al = rz/curv; d += al p; r += al hp;
    If[Norm[r] <= tol g0, Break[]];
    z = r/dg; rzn = r . z; p = -z + (rzn/rz) p; rz = rzn,
    {j, maxit}];
  d];

(* 1D matrices of direction i at the two Gauss points of every cell: values and derivatives of the hat functions *)
ndGauss1D[spec_, i_] := Module[{m = spec["Cells"][[i]], h = N[spec["h"][[i]]], g = {0.5 - 0.5/Sqrt[3.], 0.5 + 0.5/Sqrt[3.]}, pos},
  pos = Transpose[{Flatten[Table[2 (c - 1) + q, {c, m}, {q, 2}, {2}]], Flatten[Table[{c, c + 1}, {c, m}, {q, 2}]]}];
  {SparseArray[pos -> Flatten[Table[{1. - g[[q]], g[[q]]}, {c, m}, {q, 2}]], {2 m, m + 1}],
   SparseArray[pos -> Flatten[Table[{-1./h, 1./h}, {c, m}, {q, 2}]], {2 m, m + 1}]}];

(* heuristic: a local minimiser of the Q1 energy with two-point Gauss quadrature, by Newton's method from u0.
   For convex W the Hessian is positive definite and the Newton system is solved as it is.  Otherwise the Newton
   direction comes from conjugate gradients truncated at negative curvature (plain Newton steps converge to saddle
   points as readily as to minima).  If a direction is not a descent direction, the Hessian with the positive part
   of W'' is used.  Every step is accepted by an Armijo test, so the energy decreases.  Any nodal array is a valid
   candidate: this one only affects how tight the upper bound is. *)
ndSolveCandidate[spec_, u0_] := Module[
  {n = spec["n"], k = N[spec["k"]], wq = N[spec["CellVolume"]]/2^spec["n"], wf = N[spec["W"]], w1 = N[spec["W1"]], w2 = N[spec["W2"]],
   mats, pg, pgT, gm, stiff, free, u, ug, grad, hess, pgf, pgfT, sff, d, e0, e1, t, un, it = 0, en, slope},
  mats = Table[ndGauss1D[spec, i], {i, n}];
  pg = KroneckerProduct @@ mats[[All, 1]];
  gm = Table[KroneckerProduct @@ ReplacePart[mats[[All, 1]], i -> mats[[i, 2]]], {i, n}];
  stiff = 2 k wq Total[Transpose[#] . # & /@ gm];
  free = Flatten[Position[Flatten[spec["Fixed"]], 0]];
  u = Flatten[u0];
  If[free === {}, Return[u0]];
  pgT = Transpose[pg]; pgf = pg[[All, free]]; pgfT = Transpose[pgf]; sff = stiff[[free, free]];
  en[v_] := wq Total[polyAt[wf, pg . v]] + 0.5 v . (stiff . v);
  e0 = en[u];
  While[it++ < If[TrueQ[spec["Convex"]], 100, 400],
    ug = pg . u;
    grad = (wq pgT . polyAt[w1, ug] + stiff . u)[[free]];
    hess = sff + wq pgfT . (DiagonalMatrix[SparseArray[polyAt[w2, ug]]] . pgf);
    d = If[TrueQ[spec["Convex"]], ndSolveSparse[hess, -grad], ndTruncatedCG[hess, grad]];
    slope = If[VectorQ[d, MachineNumberQ], grad . d, 0.];
    If[! (slope < 0.),
      hess = sff + wq pgfT . (DiagonalMatrix[SparseArray[Ramp[polyAt[w2, ug]] + 10.^-12]] . pgf);
      d = ndSolveSparse[hess, -grad];
      If[d === $Failed, Break[]];
      slope = grad . d];
    If[! (slope < 0.) || Max[Abs[d]] <= 10.^-13 (1. + Max[Abs[u]]), Break[]];
    t = 1.; un = u; un[[free]] = u[[free]] + d; e1 = en[un];
    While[! (e1 <= e0 + 10.^-4 t slope) && t > 10.^-8, t *= 0.5; un[[free]] = u[[free]] + t d; e1 = en[un]];
    If[! (e1 < e0), Break[]];
    u = un;
    If[e0 - e1 <= 10.^-15 Abs[e0], Break[]];
    e0 = e1];
  ArrayReshape[u, spec["Cells"] + 1]];

(* RIGOROUS enclosure {lo, hi} (exact rationals) of the energy of the Q1 function with nodal values u.
   Gradient term: exact quadratic form of the nodal values.  Potential term: midpoint rule on each cell with the
   remainder |c|/24 Sum_i h_i^2 sup |W''(u_h)| (d_i u_h)^2  (docs/SCALAR_ND.md, Proposition 7). *)
ndUpperBound[spec_, u0_] := ndQuiet@Module[
  {n = spec["n"], u, verts, uc, invh, h2, umin, umax, w2, w2pm, kin, remLo, remHi, corners, d, acc, dmax2, ri, mean, wc, cellLo, cellHi},
  verts = Tuples[{0, 1}, n];
  u = ndRealArray[u0, spec["Cells"] + 1];
  If[u === $Failed, Throw[failure["BadCandidate", "The nodal values must be a real array of dimensions Cells + 1."], "CertiField"]];
  If[! (Max[Abs[spec["Fixed"] (u - spec["G"])]] == 0.),
    Throw[failure["Internal", "The Dirichlet nodal values are not imposed exactly."], "CertiField"]];
  uc = Association[# -> ivPt[vertexPart[u, n, #]] & /@ verts];
  invh = ivC[1/#] & /@ spec["h"];
  h2 = ivC[#^2] & /@ spec["h"];
  umin = Fold[vMin, uc[#][[1]] & /@ verts];
  umax = Fold[vMax, uc[#][[2]] & /@ verts];
  w2 = ivHornerQ[spec["W2"], {umin, umax}];
  w2pm = {vMin[w2[[1]], 0.] + 0. umin, vMax[w2[[2]], 0.] + 0. umin};
  corners = Tuples[{0, 1}, n - 1];
  kin = None; remLo = 0. umin; remHi = 0. umin;
  Do[
    d = Association[# -> ndMul[ndSub[uc[Insert[#, 1, i]], uc[Insert[#, 0, i]]], invh[[i]]] & /@ corners];
    acc = ivSum[Flatten[Table[
      ivK[If[e === e2, ivSqr[d[e]], ndMul[d[e], d[e2]]], Times @@ MapThread[If[#1 == #2, 1/3, 1/6] &, {e, e2}]],
      {e, corners}, {e2, corners}], 1]];
    kin = If[kin === None, acc, ndAdd[kin, acc]];
    dmax2 = Fold[vMax, ivSqr[d[#]][[2]] & /@ corners];
    ri = ivK[ndMul[ndMul[w2pm, {0. dmax2, dmax2}], h2[[i]]], 1/24];
    remLo = ndDn[remLo + ri[[1]]]; remHi = ndUp[remHi + ri[[2]]],
    {i, n}];
  kin = ivK[kin, spec["k"]];
  mean = ivK[ivSum[uc /@ verts], 1/2^n];
  wc = ivHornerQ[spec["W"], mean];
  cellLo = ndDn[ndDn[kin[[1]] + wc[[1]]] + remLo];
  cellHi = ndUp[ndUp[kin[[2]] + wc[[2]]] + remHi];
  spec["CellVolume"] {exactTotal[cellLo], exactTotal[cellHi]}];

(* ------------------------------------------------------------------------------------------------ *)
(* fields                                                                                            *)
(* a field is <|"TRange" -> {tlo, thi}, "K" -> K, "Affine" -> True | False,                           *)
(*              "S0" -> {array over the faces of direction i, ...}, "C" -> {faces x (K + 2), ...}|>     *)
(* K = 0 (tlo == thi): the affine part alone.                                                          *)
(* ------------------------------------------------------------------------------------------------ *)

faceShape[spec_, i_] := ReplacePart[spec["Cells"], i -> spec["Cells"][[i]] + 1];

(* the field with its arrays packed, or $Failed if it does not fit the grid *)
ndPackField[spec_, f_] := Module[{n = spec["n"], kk, tr, s0, cc},
  If[! AssociationQ[f], Return[$Failed]];
  kk = Lookup[f, "K", $Failed]; tr = Lookup[f, "TRange", $Failed];
  If[! (IntegerQ[kk] && kk >= 0 && MatchQ[tr, {_Integer | _Rational, _Integer | _Rational}] && If[kk == 0, tr[[1]] == tr[[2]], tr[[1]] < tr[[2]]]),
    Return[$Failed]];
  If[! (ListQ[Lookup[f, "S0", $Failed]] && Length[f["S0"]] == n && ListQ[Lookup[f, "C", $Failed]] && Length[f["C"]] == n), Return[$Failed]];
  s0 = Table[ndRealArray[f["S0"][[i]], faceShape[spec, i]], {i, n}];
  cc = If[kk == 0, ConstantArray[{}, n], Table[ndRealArray[f["C"][[i]], Append[faceShape[spec, i], kk + 2]], {i, n}]];
  If[MemberQ[s0, $Failed] || MemberQ[cc, $Failed], Return[$Failed]];
  <|"TRange" -> tr, "K" -> kk, "Affine" -> TrueQ[Lookup[f, "Affine", False]], "S0" -> s0, "C" -> cc|>];

ndFieldQ[spec_, f_] := AssociationQ[ndPackField[spec, f]];

(* the normal component must vanish on the natural faces *)
ndRespectsNaturalFaces[spec_, f_] := And @@ Flatten[Table[
  With[{sp = If[side == 1, 1, -1]},
    spec["Dirichlet"][[i, side]] || (Max[Abs[sliceAt[f["S0"][[i]], spec["n"], i, sp]]] == 0. &&
      (f["K"] == 0 || Max[Abs[sliceAt[f["C"][[i]], spec["n"], i, sp]]] == 0.))],
  {i, spec["n"]}, {side, 2}]];

(* ------------------------------------------------------------------------------------------------ *)
(* restricted convex conjugates  T(d) = sup { t d - W(t) : t >= t0 } (side +1) or t <= t0 (side -1)    *)
(* ------------------------------------------------------------------------------------------------ *)

(* RIGOROUS upper bounds (machine numbers) for T at the machine numbers dvals; W superlinear.
   On a subinterval of width w: t d - W(t) <= larger end value + w^2/8 sup (W'')^+.  Beyond R the function is
   non-increasing because W'(t) >= |d| there (docs/SCALAR_ND.md, Lemma 4). *)
Options[conjUpper] = {"Fine" -> 6000, "Coarse" -> 1500};
conjUpper[W0_List, dvals_List, t00_, side_, OptionsPattern[]] := ndQuiet@Module[
  {W = W0, dv, t0 = t00, deg, dmax, R, arg, a, b, W2, best, segs, w, ts, tl, th, wlo, w2hi, corr, plo, phi, f, seg},
  dv = ndRealArray[dvals, {Length[dvals]}];
  If[dv === $Failed, ndOverflow[]];
  If[side < 0, W = W (-1)^Range[0, Length[W] - 1]; dv = -dv; t0 = -t0];
  deg = Length[W] - 1;
  dmax = Max[Abs[SetPrecision[dv, Infinity]]];
  R = Max[1, t0, (Total[Range[1, deg - 1] Abs[W[[2 ;; deg]]]] + dmax)/(deg Last[W])] + 1;
  (* where the maxima are (machine arithmetic; only used to place the fine part of the partition) *)
  arg = conjArgmax[N[W], N[#], N[t0]] & /@ {Min[dv], Max[dv], Median[dv]};
  With[{span = Max[Max[arg] - Min[arg], 10.^-3 N[R - t0], 10.^-6]},
    a = Max[t0, Rationalize[Min[arg] - 0.25 span, 10^-9]];
    b = Min[R, Rationalize[Max[arg] + 0.25 span, 10^-9]]];
  If[! (t0 <= a < b <= R), a = t0; b = R];
  W2 = polyD[polyD[W]];
  best = ConstantArray[-$big, Length[dv]];
  segs = Select[{{t0, a, OptionValue["Coarse"]}, {a, b, OptionValue["Fine"]}, {b, R, OptionValue["Coarse"]}}, #[[2]] > #[[1]] &];
  Do[
    w = (sg[[2]] - sg[[1]])/sg[[3]];
    ts = sg[[1]] + w Range[0, sg[[3]]];
    (* freshly converted exact numbers: packing them is not packing a result of machine arithmetic *)
    tl = Developer`ToPackedArray[fLo /@ ts]; th = Developer`ToPackedArray[fHi /@ ts];
    wlo = Developer`ToPackedArray[fLo[polyAt[W, #]] & /@ ts];
    w2hi = Developer`ToPackedArray[Table[fHi[Max[0, Last[exactRange[W2, ts[[j]], ts[[j + 1]]]]]], {j, sg[[3]]}]];
    corr = ndUp[fHi[w^2/8] w2hi];
    plo = ndChk1[Outer[Times, dv, tl]]; phi = ndChk1[Outer[Times, dv, th]];
    f = ndUp[Transpose[Transpose[ndUp[vMax[plo, phi]]] - wlo]];
    seg = ndUp[Transpose[Transpose[vMax[f[[All, ;; -2]], f[[All, 2 ;;]]]] + corr]];
    best = vMax[best, Max /@ seg],
    {sg, segs}];
  best];

(* machine estimate of the maximiser of t d - W(t) on t >= t0 *)
conjArgmax[Wf_List, d_, t0_] := Module[{t, cands, vals},
  cands = Join[{t0}, Select[t /. NSolve[(Rest[Wf] Range[Length[Wf] - 1]) . t^Range[0, Length[Wf] - 2] == d, t], (Abs[Im[#]] < 10.^-9 (1 + Abs[Re[#]]) && Re[#] > t0) &] // Re];
  vals = (# d - polyAt[Wf, #]) & /@ cands;
  cands[[First[Ordering[vals, -1]]]]];

(* RIGOROUS upper bound, cell by cell, for T over the intervals {dlo, dhi} (arrays).  T is convex: on an interval it
   is at most the larger of its end values, and between two table nodes it is at most the chord through the node
   bounds. *)
conjCells[W_List, {dlo_, dhi_}, t0_, side_, nd_ : 513] := ndQuiet@Module[{lo = Min[dlo], hi = Max[dhi], pad, g, tg, at},
  pad = Max[10.^-9, 10.^-6 Max[Abs[lo], Abs[hi], hi - lo]];
  g = Subdivide[lo - pad, hi + pad, nd - 1];
  If[! (Developer`PackedArrayQ[g, Real] && g[[1]] < lo && g[[-1]] > hi && Min[Differences[g]] > 0), Throw[failure["Internal", "Conjugate table nodes."], "CertiField"]];
  tg = conjUpper[W, g, t0, side];
  at[p_] := Module[{dims = Dimensions[p], pf = Flatten[p], j, gj, gj1, lam},
    j = Clip[Floor[(pf - g[[1]])/(g[[2]] - g[[1]])] + 1, {1, nd - 1}];
    (* repair the bracket g[[j]] <= p <= g[[j + 1]] exactly (comparisons of machine numbers are exact) *)
    j = j - UnitStep[g[[j]] - pf] (1 - UnitStep[pf - g[[j]]]);
    j = Clip[j, {1, nd - 1}];
    j = j + (1 - UnitStep[g[[j + 1]] - pf]);
    j = Clip[j, {1, nd - 1}];
    gj = g[[j]]; gj1 = g[[j + 1]];
    If[! (Min[pf - gj] >= 0 && Min[gj1 - pf] >= 0), Throw[failure["Internal", "Conjugate table bracket."], "CertiField"]];
    lam = ndDiv[ndSub[ivPt[pf], ivPt[gj]], ndSub[ivPt[gj1], ivPt[gj]]];
    lam = {Clip[lam[[1]], {0., 1.}], Clip[lam[[2]], {0., 1.}]};
    ArrayReshape[ndAdd[ivPt[tg[[j]]], ndMul[lam, ndSub[ivPt[tg[[j + 1]]], ivPt[tg[[j]]]]]][[2]], dims]];
  vMax[at[dlo], at[dhi]]];

(* RIGOROUS lower bound for inf { W(t) : t >= t0 } (side +1) or t <= t0 (side -1); -$big if W is unbounded below.
   Remembered for the session: the value depends only on exact data, and calibration mode asks for it often. *)
tailInf[W_List, t0_, side_] := tailInf[W, t0, side] = Which[
  Length[W] == 1, fLo[W[[1]]],
  superlinearQ[W], -First[conjUpper[W, {0.}, t0, side]],
  True, -$big];

(* ------------------------------------------------------------------------------------------------ *)
(* RIGOROUS: the boundary flux                                                                        *)
(* ------------------------------------------------------------------------------------------------ *)

(* enclosure of the clamped spline at thin intervals t = {tlo, thi arrays}: cb has dimensions faces x (K + 2) *)
splineEval[cb0_, kk_, {t0_, t1_}, {tl0_, th0_}] := ndQuiet@Module[
  {dl = (t1 - t0)/kk, cb, tl, th, z, zlo, zhi, kap, dims, cf, nf, pos, piece, th1, r, over, th2, r2},
  tl = Developer`ToPackedArray[N[tl0]]; th = Developer`ToPackedArray[N[th0]]; dims = Dimensions[tl];
  cb = ndRealArray[cb0, Append[dims, kk + 2]];
  If[cb === $Failed || ! Developer`PackedArrayQ[tl, Real] || ! Developer`PackedArrayQ[th, Real], ndOverflow[]];
  z = ivK[ndSub[{tl, th}, ivC[t0]], 1/dl];
  zlo = Clip[z[[1]], {0., N[kk]}]; zhi = Clip[z[[2]], {0., N[kk]}];
  If[Max[zhi - zlo] >= 0.5, Throw[failure["Internal", "splineEval is for thin intervals."], "CertiField"]];
  kap = Clip[Floor[zlo], {0, kk - 1}];
  cf = ArrayReshape[cb, {Times @@ dims, kk + 2}]; nf = Times @@ dims;
  pos[kq_, off_] := ArrayReshape[Extract[cf, Transpose[{Range[nf], Flatten[kq] + 1 + off}]], dims];
  piece[kq_, tt_] := Module[{n0, n2, n1},
    n0 = ivK[ivSqr[ndSub[ivPt[1. + 0. tl], tt]], 1/2];
    n2 = ivK[ivSqr[tt], 1/2];
    n1 = ndSub[ndSub[ivPt[1. + 0. tl], n0], n2];
    ivSum[{ndMul[ivPt[pos[kq, 0]], n0], ndMul[ivPt[pos[kq, 1]], n1], ndMul[ivPt[pos[kq, 2]], n2]}]];
  th1 = {Clip[rDn[zlo - kap], {0., 1.}], Clip[rUp[zhi - kap], {0., 1.}]};
  r = piece[kap, th1];
  over = gtMask[rUp[zhi - kap], 1.] UnitStep[(kk - 2) - kap];     (* the last piece has no successor *)
  If[Total[over, Infinity] > 0,
    th2 = {0. tl, Clip[rUp[zhi - kap - 1.], {0., 1.}]};
    r2 = piece[Clip[kap + 1, {0, kk - 1}], th2];
    r = {vWhere[over, vMin[r[[1]], r2[[1]]], r[[1]]], vWhere[over, vMax[r[[2]], r2[[2]]], r[[2]]]}];
  r];

(* the 2^(n-1) arrays of nodal g values at the vertices of the boundary faces of one side *)
faceVertexValues[spec_, i_, side_] := Module[{n = spec["n"], slab},
  slab = sliceAt[spec["G"], n, i, If[side == 1, 1, -1]];
  vertexPart[slab, n - 1, #] & /@ Tuples[{0, 1}, n - 1]];

(* RIGOROUS lower bound (exact rational) for the flux through the Dirichlet faces (docs/SCALAR_ND.md, Lemma 5).
   The Lipschitz and curvature constants of Lemma 5 range over the knot intervals met by the datum on a sub-face.
   They are found by exact comparisons with enclosures of the knots: interval kap = [t_kap, t_kap+1] is included
   whenever (lower enclosure of t_kap) < gmax and (upper enclosure of t_kap+1) > gmin, a superset of the intervals
   whose interior meets [gmin, gmax].  (The reference takes floor((g - t_lo)/Delta) -1/+1 in floating point, a
   larger set; see PORTING.md.) *)
ndFluxLower[spec_, f_, r_ : 8] := ndQuiet@Module[
  {n = spec["n"], kk = f["K"], t0, t1, total = 0, vv, sign, area, gbar, faceLo, cb, dl, c0, c1, c2, d01, d12, s2, s2abs, s1abs,
   const, rr, corners, subLo, sv, loc, gmin, gmax, m, val, dev, knLo, knHi, inr, l1, l2, inside, pen2, pen1, pen},
  {t0, t1} = f["TRange"];
  corners = Tuples[{0, 1}, n - 1];
  If[kk > 0,
    dl = (t1 - t0)/kk;
    knLo = Developer`ToPackedArray[fLo /@ (t0 + Range[0, kk] dl)]; knHi = Developer`ToPackedArray[fHi /@ (t0 + Range[0, kk] dl)]];
  Do[
    If[spec["Dirichlet"][[i, side]],
      vv = faceVertexValues[spec, i, side];
      sign = If[side == 1, -1, 1];
      area = spec["CellVolume"]/spec["h"][[i]];
      gbar = ivK[ivSum[ivPt /@ vv], 1/Length[vv]];
      faceLo = ivK[ndMul[ivPt[sliceAt[f["S0"][[i]], n, i, If[side == 1, 1, -1]]], gbar], sign][[1]];        (* affine part *)
      If[kk > 0,
        cb = sliceAt[f["C"][[i]], n, i, If[side == 1, 1, -1]];                        (* faces x (K + 2) *)
        c0 = lastPart[cb, n - 1, Span[1, -3]]; c1 = lastPart[cb, n - 1, Span[2, -2]]; c2 = lastPart[cb, n - 1, Span[3, All]];
        d01 = ivK[ndSub[ivPt[c1], ivPt[c0]], 1/dl];
        d12 = ivK[ndSub[ivPt[c2], ivPt[c1]], 1/dl];
        s2 = ivK[ndSub[d12, d01], 1/dl];
        s2abs = vMax[Abs[s2[[1]]], Abs[s2[[2]]]];                                     (* |s''| on each knot interval *)
        s1abs = Fold[vMax, {Abs[d01[[1]]], Abs[d01[[2]]], Abs[d12[[1]]], Abs[d12[[2]]]}];
        const = Times @@ (1 - Unitize[# - First[vv]] & /@ Rest[vv]);
        rr = If[Min[const] == 1, 1, r];
        subLo = 0. First[vv];
        Do[
          sv = Table[
            loc = (sub + ev)/rr;
            ivSum[Table[
              With[{wgt = Times @@ MapThread[If[#2 == 1, #1, 1 - #1] &, {loc, corners[[b]]}]},
                If[wgt == 0, Nothing, ivK[ivPt[vv[[b]]], wgt]]], {b, Length[corners]}]],
            {ev, corners}];
          gmin = Fold[vMin, sv[[All, 1]]]; gmax = Fold[vMax, sv[[All, 2]]];
          m = ivK[ivSum[sv], 1/Length[sv]];                                         (* sub-face mean of the datum *)
          val = ivK[splineEval[cb, kk, {t0, t1}, m], sign];
          dev = vMax[ndUp[gmax - m[[1]]], ndUp[m[[2]] - gmin]];
          (* knot intervals met by the datum on this sub-face (signs of differences of machine numbers are exact) *)
          inr = (1 - UnitStep[Outer[Plus, -gmax, Most[knLo]]]) (1 - UnitStep[Outer[Plus, gmin, -Rest[knHi]]]);
          l2 = Map[Max, inr s2abs, {n - 1}];
          l1 = Map[Max, inr s1abs, {n - 1}];
          inside = geMask[gmin, fHi[t0]] leMask[gmax, fLo[t1]];
          pen2 = ivK[ndMul[ivPt[l2], ivSqr[ivPt[dev]]], 1/2][[2]];
          pen1 = ndMul[ivPt[l1], ivPt[dev]][[2]];
          pen = (1 - const) vWhere[inside, vMin[pen2, pen1], pen1];
          subLo = ndDn[subLo + ndDn[val[[1]] - pen]],
          {sub, Tuples[Range[0, rr - 1], n - 1]}];
        faceLo = ndDn[faceLo + ivK[ivPt[subLo], 1/rr^(n - 1)][[1]]]];
      total += area exactTotal[faceLo]],
    {i, n}, {side, 2}];
  total];

(* ------------------------------------------------------------------------------------------------ *)
(* RIGOROUS: the lower bound for a given field (docs/SCALAR_ND.md, Section 3)                          *)
(* ------------------------------------------------------------------------------------------------ *)

Options[ndLowerBound] = {"Subintervals" -> 16};
ndLowerBound[spec_, f0_, OptionsPattern[]] := ndQuiet[ndLowerBoundRaw[spec, f0, OptionValue["Subintervals"]]];
ndLowerBoundRaw[spec_, f0_, qs_] := Module[
  {n = spec["n"], f, kk, t0, t1, W = spec["W"], invh, inv4k, s0, d0, a0, tl, th, flux, verts,
   mup, dl, w, xi, cc, c, d01, d12, dd2, prev, cur, tt, theta, n0, n1, n2, w0, w1, dv, aa, base, ta, w2hi, corr, dend, intM, mc},
  f = ndPackField[spec, f0];
  If[f === $Failed || ! ndRespectsNaturalFaces[spec, f], Throw[failure["BadField", "The field does not fit the grid or has a normal component on a natural face."], "CertiField"]];
  If[! (IntegerQ[qs] && qs >= 1), Throw[failure["BadOption", "\"Subintervals\" must be a positive integer."], "CertiField"]];
  kk = f["K"]; {t0, t1} = f["TRange"];
  If[! f["Affine"] && Max[Abs[Flatten[f["S0"]]]] != 0., Throw[failure["BadField", "A field without affine part must have S0 = 0."], "CertiField"]];
  If[f["Affine"] && ! spec["Superlinear"], Throw[failure["Inconclusive", "The affine part needs a potential of even degree with positive leading coefficient."], "CertiField"]];
  verts = Tuples[{0, 1}, n];
  invh = ivC[1/#] & /@ spec["h"];
  inv4k = ivC[1/(4 spec["k"])];
  s0 = ivPt /@ f["S0"];
  d0 = ivSum[Table[ndMul[ndSub[ivHi[s0[[i]], n, i], ivLo[s0[[i]], n, i]], invh[[i]]], {i, n}]];
  a0 = Table[With[{q = ndMul[ivSqr[s0[[i]]], inv4k]}, {ivLo[q, n, i], ivHi[q, n, i]}], {i, n}];
  If[f["Affine"],
    tl = conjCells[W, d0, t0, -1]; th = conjCells[W, d0, t1, +1],
    tl = tailInf[W, t0, -1]; th = tailInf[W, t1, +1];
    If[tl <= -$big || th <= -$big, Throw[failure["Inconclusive", "The potential is not bounded below outside the t-range."], "CertiField"]];
    tl = -tl + 0. d0[[1]]; th = -th + 0. d0[[1]]];
  flux = ndFluxLower[spec, f];
  If[kk == 0,
    (* the affine part alone: sup_t Phi = |sigma0(x)|^2/(4k) + W*(d0); the cell integral of the first term is exact *)
    mc = ivSum[Table[
      With[{lo = ivLo[s0[[i]], n, i], hi = ivHi[s0[[i]], n, i]},
        ndMul[ivK[ivSum[{ivSqr[lo], ndMul[lo, hi], ivSqr[hi]}], 1/3], inv4k]], {i, n}]];
    intM = spec["CellVolume"] exactTotal[ndUp[mc[[2]] + vMax[tl, th]]];
    Return[flux - intM]];
  dl = (t1 - t0)/kk;
  w = dl/qs;
  xi = t0 + (Range[0, kk + 1] - 1/2) dl;                      (* Greville abscissae: Sum[xi_m b_m(t)] = t *)
  cc = f["C"];
  (* divergence of the spline part at the two ends of the t-range *)
  dend = Table[
    ivSum[Table[
      With[{ev = ivK[ndAdd[ivPt[lastPart[cc[[i]], n, m0]], ivPt[lastPart[cc[[i]], n, m0 + 1]]], 1/2]},
        ndMul[ndSub[ivHi[ev, n, i], ivLo[ev, n, i]], invh[[i]]]], {i, n}]],
    {m0, {1, kk + 1}}];
  mup = Association[# -> None & /@ verts];
  Do[
    c = Table[ndAdd[ivPt[lastPart[cc[[i]], n, kap + r - 1]], ivK[s0[[i]], xi[[kap + r - 1]]]], {i, n}, {r, 3}];
    d01 = Table[ndSub[c[[i, 2]], c[[i, 1]]], {i, n}];
    d12 = Table[ndSub[c[[i, 3]], c[[i, 2]]], {i, n}];
    dd2 = ivSum[Table[
      With[{s2 = ivK[ndSub[d12[[i]], d01[[i]]], 1/dl^2]}, ndMul[ndSub[ivHi[s2, n, i], ivLo[s2, n, i]], invh[[i]]]], {i, n}]];
    prev = None;
    Do[
      theta = j/qs;
      tt = t0 + (kap - 1 + theta) dl;
      {n0, n1, n2} = {(1 - theta)^2/2, (-2 theta^2 + 2 theta + 1)/2, theta^2/2};
      {w0, w1} = {1 - theta, theta};
      dv = ivSum[Table[
        With[{s = ivSum[{ivK[c[[i, 1]], n0], ivK[c[[i, 2]], n1], ivK[c[[i, 3]], n2]}]},
          ndMul[ndSub[ivHi[s, n, i], ivLo[s, n, i]], invh[[i]]]], {i, n}]];
      aa = Table[
        With[{q = ndMul[ivSqr[ivK[ndAdd[ivK[d01[[i]], w0], ivK[d12[[i]], w1]], 1/dl]], inv4k]}, {ivLo[q, n, i], ivHi[q, n, i]}], {i, n}];
      base = ndSub[dv, ivC[polyAt[W, tt]]];
      cur = Association[# -> ivSum[Join[{base}, Table[aa[[i, #[[i]] + 1]], {i, n}]]][[2]] & /@ verts];
      If[prev =!= None,
        ta = t0 + (kap - 1 + (j - 1)/qs) dl;
        w2hi = fHi[Last[exactRange[spec["W2"], ta, tt]]];
        corr = ivK[ivPos[ndSub[ivPt[w2hi + 0. dd2[[1]]], dd2]], w^2/8][[2]];
        Do[
          With[{cand = ndUp[vMax[prev[v], cur[v]] + corr]},
            mup[v] = If[mup[v] === None, cand, vMax[mup[v], cand]]], {v, verts}]];
      prev = cur,
      {j, 0, qs}],
    {kap, kk}];
  Do[
    With[{av = ivSum[Table[a0[[i, v[[i]] + 1]], {i, n}]]},
      mup[v] = vMax[mup[v], ivSum[{dend[[1]], av, ivPt[tl]}][[2]]];
      mup[v] = vMax[mup[v], ivSum[{dend[[2]], av, ivPt[th]}][[2]]]], {v, verts}];
  intM = spec["CellVolume"]/2^n exactTotal[ndChk1 /@ Values[mup]];
  flux - intM];

(* ------------------------------------------------------------------------------------------------ *)
(* heuristic: the affine part from a cell-centred finite-volume solve (convex mode)                   *)
(* ------------------------------------------------------------------------------------------------ *)

(* face means of the datum on one side (an array over the boundary faces) *)
faceMeans[spec_, i_, side_] := Mean[faceVertexValues[spec, i, side]];

(* 1D finite-volume matrix of direction i: 2k/h^2 between neighbours, 4k/h^2 to a Dirichlet face *)
fvMatrix1D[spec_, i_] := Module[{m = spec["Cells"][[i]], c = N[2 spec["k"]/spec["h"][[i]]^2], diag},
  diag = ConstantArray[2 c, m];
  If[m >= 1, diag[[1]] = c; diag[[m]] = c; If[m == 1, diag[[1]] = 0.]];
  If[spec["Dirichlet"][[i, 1]], diag[[1]] += 2 c];
  If[spec["Dirichlet"][[i, 2]], diag[[m]] += 2 c];
  If[m == 1, SparseArray[{{1, 1} -> diag[[1]]}, {1, 1}],
    SparseArray[{Band[{1, 1}] -> diag, Band[{1, 2}] -> -c, Band[{2, 1}] -> -c}, {m, m}]]];

(* minimise  Sum_faces k |c| (dp/h)^2 + Dirichlet faces 2k |c| ((g - p)/h)^2 + Sum |c| W(p)  by damped Newton steps
   with the positive part of W'' (always a descent method).  Any p gives a valid field; p only affects tightness. *)
ndFVSolve[spec_, p0_] := Module[
  {n = spec["n"], dims = spec["Cells"], lap, b, wf = N[spec["W"]], w1 = N[spec["W1"]], w2 = N[spec["W2"]], p, en, grad, hess, dp, t, pn, it = 0, e0, e1, slope},
  lap = Total[Table[KroneckerProduct @@ Table[If[j == i, fvMatrix1D[spec, i], IdentityMatrix[dims[[j]], SparseArray]], {j, n}], {i, n}]];
  b = ConstantArray[0., dims];
  Do[
    If[spec["Dirichlet"][[i, side]],
      b += ArrayPad[asSlab[N[4 spec["k"]/spec["h"][[i]]^2] faceMeans[spec, i, side], i],
        ReplacePart[ConstantArray[{0, 0}, n], i -> If[side == 1, {0, dims[[i]] - 1}, {dims[[i]] - 1, 0}]]]],
    {i, n}, {side, 2}];
  b = Flatten[b];
  en[q_] := 0.5 q . (lap . q) - b . q + Total[polyAt[wf, q]];
  p = Flatten[N[p0]];
  e0 = en[p];
  While[it++ < 60,
    grad = lap . p - b + polyAt[w1, p];
    hess = lap + DiagonalMatrix[SparseArray[Ramp[polyAt[w2, p]] + 10.^-12]];
    dp = ndSolveSparse[hess, -grad];
    If[dp === $Failed, Break[]];
    slope = grad . dp;
    (* stop when the step is at rounding level: the gradient cannot be reduced further *)
    If[! (slope < 0.) || Max[Abs[dp]] <= 10.^-13 (1. + Max[Abs[p]]), Break[]];
    t = 1.; e1 = en[pn = p + dp];
    While[! (e1 <= e0 + 10.^-4 t slope) && t > 10.^-8, t *= 0.5; e1 = en[pn = p + t dp]];
    If[! (e1 < e0), Break[]];
    p = pn;
    If[e0 - e1 <= 10.^-15 Abs[e0], Break[]];
    e0 = e1];
  ArrayReshape[p, dims]];

(* the two-point flux field of cell values p: sigma0 . e_i = 2k (p_R - p_L)/h_i on interior faces,
   sigma0 . nu = 4k (g - p)/h_i on Dirichlet faces, 0 on natural faces *)
ndConvexField[spec_, p_] := Module[{n = spec["n"], k = N[spec["k"]], hf = N[spec["h"]], s0},
  s0 = Table[
    Join[
      asSlab[If[spec["Dirichlet"][[i, 1]], -(4 k/hf[[i]]) (faceMeans[spec, i, 1] - sliceAt[p, n, i, 1]), 0. sliceAt[p, n, i, 1]], i],
      (2 k/hf[[i]]) (hiFace[p, n, i] - loFace[p, n, i]),
      asSlab[If[spec["Dirichlet"][[i, 2]], (4 k/hf[[i]]) (faceMeans[spec, i, 2] - sliceAt[p, n, i, -1]), 0. sliceAt[p, n, i, -1]], i],
      i],
    {i, n}];
  <|"TRange" -> {0, 0}, "K" -> 0, "Affine" -> True, "S0" -> (Developer`ToPackedArray[N[#]] & /@ s0), "C" -> ConstantArray[{}, n]|>];

(* ------------------------------------------------------------------------------------------------ *)
(* heuristic: the optimiser of calibration mode (the spline part alone)                               *)
(* ------------------------------------------------------------------------------------------------ *)

(* per knot interval, qo + 1 samples (both ends): sample points, value / first-derivative matrices of the K + 2
   B-splines (samples x coefficients), second-derivative matrix (knot intervals x coefficients) *)
ndSplineTables[kk_, t0_, t1_, qo_] := Module[{dl = N[(t1 - t0)/kk], q = kk (qo + 1), tq, b, db, b2, th, row},
  tq = Flatten[Table[N[t0] + (kap + j/qo) dl, {kap, 0, kk - 1}, {j, 0, qo}]];
  b = ConstantArray[0., {q, kk + 2}]; db = b; b2 = ConstantArray[0., {kk, kk + 2}];
  Do[
    b2[[kap + 1, kap + 1 ;; kap + 3]] = {1., -2., 1.}/dl^2;
    Do[
      th = N[j/qo]; row = kap (qo + 1) + j + 1;
      b[[row, kap + 1 ;; kap + 3]] = {(1 - th)^2/2, (-2 th^2 + 2 th + 1)/2, th^2/2};
      db[[row, kap + 1 ;; kap + 3]] = {th - 1, 1 - 2 th, th}/dl,
      {j, 0, qo}],
    {kap, 0, kk - 1}];
  {tq, b, db, b2}];

(* values of the K + 2 clamped B-splines at an array of points: (number of points) x (K + 2) *)
ndBasisRows[kk_, t0_, t1_, t_] := Module[{dl = N[(t1 - t0)/kk], tf = Flatten[{t}], z, kap, th, nf},
  nf = Length[tf];
  z = Clip[(tf - N[t0])/dl, {0., N[kk]}];
  kap = Clip[Floor[z], {0, kk - 1}]; th = z - kap;
  Normal[SparseArray[
    Join[Transpose[{Range[nf], kap + 1}], Transpose[{Range[nf], kap + 2}], Transpose[{Range[nf], kap + 3}]] ->
      Join[(1 - th)^2/2, (-2 th^2 + 2 th + 1)/2, th^2/2], {nf, kk + 2}]]];

(* everything the objective needs on one grid: index vectors of the lower and upper face of every cell (flattened,
   last index fastest), the mask of the natural faces, the spline tables, the tails and the flux weights *)
ndCalSetup[spec_, {t0_, t1_}, kk_, qo_ : 4] := Module[
  {n = spec["n"], cells = spec["Cells"], hf = N[spec["h"]], vol = N[spec["CellVolume"]], ncell, fs, nface, lo, hi, mask,
   tq, b, db, b2, q, w, fine, segmax, w2max, tlc, thc, fwc, gp, pts, corners, vv, sign, area, slab, gq, pk = Developer`ToPackedArray},
  ncell = Times @@ cells;
  fs = Table[faceShape[spec, i], {i, n}];
  nface = Times @@@ fs;
  lo = Table[1 + Flatten[Outer[Plus, Sequence @@ Table[Range[0, cells[[j]] - 1] (Times @@ fs[[i, j + 1 ;;]]), {j, n}]]], {i, n}];
  hi = Table[lo[[i]] + Times @@ fs[[i, i + 1 ;;]], {i, n}];
  mask = Table[
    Module[{m = ConstantArray[1., fs[[i]]]},
      If[! spec["Dirichlet"][[i, 1]], m[[Sequence @@ ReplacePart[ConstantArray[All, n], i -> 1]]] = 0.];
      If[! spec["Dirichlet"][[i, 2]], m[[Sequence @@ ReplacePart[ConstantArray[All, n], i -> -1]]] = 0.];
      Flatten[m]], {i, n}];
  {tq, b, db, b2} = ndSplineTables[kk, t0, t1, qo]; q = Length[tq];
  w = N[(t1 - t0)/kk/qo];
  (* sup W'' on every subinterval (sampled; the certified evaluation encloses it) and, per sample, over the adjacent ones *)
  fine = polyAt[N[spec["W2"]], Subdivide[N[t0], N[t1], 8 kk qo]];
  segmax = Table[Max[fine[[8 s + 1 ;; 8 s + 9]]], {s, 0, kk qo - 1}];
  w2max = Flatten[Table[Max[segmax[[Select[{kap qo + j, kap qo + j + 1}, kap qo + 1 <= # <= (kap + 1) qo &]]]], {kap, 0, kk - 1}, {j, 0, qo}]];
  tlc = tailInf[spec["W"], t0, -1]; thc = tailInf[spec["W"], t1, +1];
  If[tlc <= -$big || thc <= -$big, Throw[failure["Inconclusive", "The potential is not bounded below outside the t-range: calibration mode gives no finite lower bound."], "CertiField"]];
  (* the flux is linear in the coefficients: two-point Gauss rule on every Dirichlet face *)
  fwc = Table[ConstantArray[0., {nface[[i]], kk + 2}], {i, n}];
  gp = {0.5 - 0.5/Sqrt[3.], 0.5 + 0.5/Sqrt[3.]};
  corners = Tuples[{0, 1}, n - 1]; pts = Tuples[gp, n - 1];
  Do[
    If[spec["Dirichlet"][[i, side]],
      vv = faceVertexValues[spec, i, side]; sign = If[side == 1, -1., 1.]; area = vol/hf[[i]];
      slab = 1 + Flatten[Outer[Plus, Sequence @@ Table[
        If[j == i, {If[side == 1, 0, cells[[i]]]}, Range[0, cells[[j]] - 1]] (Times @@ fs[[i, j + 1 ;;]]), {j, n}]]];
      Do[
        gq = Sum[(Times @@ MapThread[If[#2 == 1, #1, 1 - #1] &, {xg, corners[[a]]}]) vv[[a]], {a, Length[corners]}];
        fwc[[i, slab]] += (sign area/Length[pts]) ndBasisRows[kk, t0, t1, gq],
        {xg, pts}]],
    {i, n}, {side, 2}];
  <|"n" -> n, "K" -> kk, "Qo" -> qo, "Q" -> q, "NC" -> ncell, "NF" -> nface, "h" -> hf, "k" -> N[spec["k"]], "wv" -> vol/2^n,
    "Lo" -> lo, "Hi" -> hi, "Mask" -> mask, "Off" -> Prepend[Accumulate[nface (kk + 2)], 0],
    "B" -> pk[b], "BT" -> pk[Transpose[b]], "DB" -> pk[db], "DBT" -> pk[Transpose[db]], "B2" -> pk[b2], "B2T" -> pk[Transpose[b2]],
    "Kq" -> pk[Flatten[Table[kap, {kap, kk}, {qo + 1}]]], "VqM" -> ConstantArray[pk[polyAt[N[spec["W"]], tq]], ncell],
    "W2maxM" -> ConstantArray[pk[N[w2max]], ncell], "cw" -> w^2/8, "Tl" -> -tlc, "Th" -> -thc, "FWC" -> fwc, "Verts" -> Tuples[{0, 1}, n]|>];

(* the coefficient matrices (faces x (K + 2), natural faces zeroed) of the packed vector x, and back *)
ndCalUnpack[st_, x_] := Table[
  ArrayReshape[x[[st["Off"][[i]] + 1 ;; st["Off"][[i + 1]]]], {st["NF"][[i]], st["K"] + 2}] st["Mask"][[i]], {i, st["n"]}];

(* minus the smoothed functional and its gradient.  The functional is
     J_tau = flux - Sum_cells |c| 2^-n Sum_vertices tau log-mean-exp { Phi~(kappa, j) } and the two tails ,
     Phi~ = Phi + (w^2/8) (sup W'' - D'')^+ ,   w = Delta/qo
   (the curvature allowance of Lemma 3 is part of it: without it the optimiser hides peaks between the samples).
   J_tau is concave in the coefficients.  curvature -> False is for stress tests only. *)
ndCalValueGradient[st_, x_, tau_, curvature_ : True] := Module[
  {n = st["n"], kk = st["K"], q = st["Q"], hf = st["h"], k = st["k"], wv = st["wv"], nc = st["NC"], verts = st["Verts"], it = 1./tau, cw, cm, dc, d2, ds, aa, pos, base,
   obj = 0., gtl, gth, tlv, thv, mx, ex, etl, eth, sm, wgt, exs, ga, gd, g2, gdc, grads, gi, gw},
  cw = If[TrueQ[curvature], st["cw"], 0.];
  cm = ndCalUnpack[st, x];
  dc = Sum[(cm[[i, st["Hi"][[i]]]] - cm[[i, st["Lo"][[i]]]])/hf[[i]], {i, n}];          (* cells x (K + 2): divergence coefficients *)
  d2 = dc . st["B2T"];
  ds = Table[cm[[i]] . st["DBT"], {i, n}];                                              (* faces x Q: d/dt of the normal components *)
  (* everything that goes into the exponential is divided by tau once *)
  aa = Table[With[{sq = ds[[i]]^2 (it/(4 k))}, {sq[[st["Lo"][[i]]]], sq[[st["Hi"][[i]]]]}], {i, n}];
  pos = Ramp[st["W2maxM"] - d2[[All, st["Kq"]]]];
  base = (dc . st["BT"] + cw pos - st["VqM"]) it;
  tlv = ((dc[[All, 1]] + dc[[All, 2]])/2 + st["Tl"]) it; thv = ((dc[[All, kk + 1]] + dc[[All, kk + 2]])/2 + st["Th"]) it;
  gtl = ConstantArray[0., nc]; gth = gtl;
  (* per vertex of the cells: the weights of the samples (minus the cell weight times the soft-max weights) *)
  exs = Table[
    ex = base + Sum[aa[[i, v[[i]] + 1]], {i, n}];
    mx = Max /@ ex; mx = mx + Ramp[tlv - mx]; mx = mx + Ramp[thv - mx];
    ex = Exp[ex - mx]; etl = Exp[tlv - mx]; eth = Exp[thv - mx];
    sm = Total[ex, {2}] + etl + eth;
    obj -= wv tau (Total[mx] + Total[Log[sm/(q + 2)]]);
    wgt = -wv/sm;
    gtl += etl wgt; gth += eth wgt;
    ex wgt,
    {v, verts}];
  ga = Table[Total[Pick[exs, verts[[All, i]], side]], {i, n}, {side, {0, 1}}];             (* by face of the cell *)
  gd = ga[[1, 1]] + ga[[1, 2]];
  g2 = -cw Total[ArrayReshape[gd Unitize[pos], {nc, kk, st["Qo"] + 1}], {3}];
  gdc = gd . st["B"] + g2 . st["B2"];
  gdc[[All, 1]] += gtl/2; gdc[[All, 2]] += gtl/2; gdc[[All, kk + 1]] += gth/2; gdc[[All, kk + 2]] += gth/2;
  grads = Table[
    gi = ConstantArray[0., {st["NF"][[i]], kk + 2}];
    gi[[st["Hi"][[i]]]] = gdc/hf[[i]]; gi[[st["Lo"][[i]]]] -= gdc/hf[[i]];
    gw = ConstantArray[0., {st["NF"][[i]], q}];                                           (* the weights, by face of the grid *)
    gw[[st["Lo"][[i]]]] = ga[[i, 1]]; gw[[st["Hi"][[i]]]] += ga[[i, 2]];
    obj += Total[st["FWC"][[i]] cm[[i]], 2];
    (gi + ((gw ds[[i]]) . st["DB"])/(2 k) + st["FWC"][[i]]) st["Mask"][[i]],
    {i, n}];
  {-obj, -Flatten[grads]}];

(* line search with the strong Wolfe conditions (bracketing and zoom, Nocedal and Wright, Algorithms 3.5 and 3.6).
   Returns {step, f, gradient, number of evaluations}; step 0 if no acceptable point was found. *)
ndWolfe[fg_, x_, f0_, g0_, d_, gd0_, t1_, budget_] := Module[
  {c1 = 10.^-4, c2 = 0.9, tp = 0., fp = f0, gp = gd0, t = t1, ft, gt, gtd, nev = 0, lo, hi, flo, glo, fhi, best = {0., f0, g0},
   done = False, zoom = False, den},
  Do[
    {ft, gt} = fg[x + t d]; nev++; gtd = gt . d;
    If[! (ft <= f0 + c1 t gd0) || (i > 1 && ft >= fp), {lo, flo, glo, hi, fhi} = {tp, fp, gp, t, ft}; zoom = True; Break[]];
    best = {t, ft, gt};
    If[Abs[gtd] <= -c2 gd0, done = True; Break[]];
    If[gtd >= 0., {lo, flo, glo, hi, fhi} = {t, ft, gtd, tp, fp}; zoom = True; Break[]];
    tp = t; fp = ft; gp = gtd; t = 2. t,
    {i, 1, Max[1, Min[10, budget]]}];
  If[done || ! zoom, Return[Append[best, nev]]];
  Do[
    (* minimiser of the quadratic through (lo, flo) with slope glo and (hi, fhi), kept away from the ends *)
    den = 2. (fhi - flo - glo (hi - lo));
    t = If[den > 0., lo - glo (hi - lo)^2/den, 0.5 (lo + hi)];
    t = Clip[t, Sort[{lo + 0.1 (hi - lo), lo + 0.9 (hi - lo)}]];
    {ft, gt} = fg[x + t d]; nev++; gtd = gt . d;
    If[! (ft <= f0 + c1 t gd0) || ft >= flo,
      hi = t; fhi = ft,
      best = {t, ft, gt};
      If[Abs[gtd] <= -c2 gd0, Break[]];
      If[gtd (hi - lo) >= 0., hi = lo; fhi = flo];
      lo = t; flo = ft; glo = gtd];
    If[Abs[hi - lo] <= 10.^-14 Max[Abs[lo], Abs[hi]], Break[]],
    {j, 1, 20}];
  Append[best, nev]];

(* L-BFGS (two-loop recursion, mem correction pairs) for a function fg returning {value, gradient}.  The limits are
   those of the reference (SciPy L-BFGS-B without bounds): maxit iterations, 3 maxit evaluations. *)
ndLBFGS[fg_, x0_, maxit_Integer, mem_Integer : 20] := Module[
  {x = x0, f, g, d, ss = {}, ys = {}, rhos = {}, it = 0, nev = 1, maxev = 3 maxit, q, r, al, be, gamma = 1., t, fn, gn, k, sk, yk, sy, gd0, xn},
  {f, g} = fg[x];
  While[it < maxit && nev < maxev,
    If[Max[Abs[g]] <= 10.^-12, Break[]];
    q = g; al = ConstantArray[0., Length[ss]];
    Do[al[[j]] = rhos[[j]] (ss[[j]] . q); q -= al[[j]] ys[[j]], {j, Length[ss], 1, -1}];
    r = gamma q;
    Do[be = rhos[[j]] (ys[[j]] . r); r += (al[[j]] - be) ss[[j]], {j, 1, Length[ss]}];
    d = -r; gd0 = g . d;
    If[! (gd0 < 0.), d = -g; gd0 = -(g . g); ss = {}; ys = {}; rhos = {}];
    t = If[ss === {}, Min[1., 1./Norm[d]], 1.];
    {t, fn, gn, k} = ndWolfe[fg, x, f, g, d, gd0, t, maxev - nev]; nev += k;
    If[t == 0.,
      If[ss === {}, Break[]];                                  (* steepest descent failed too: stop *)
      ss = {}; ys = {}; rhos = {}; Continue[]];               (* forget the curvature pairs and try again *)
    xn = x + t d; sk = xn - x; yk = gn - g; sy = sk . yk;
    If[sy > 10.^-10 Norm[sk] Norm[yk],
      AppendTo[ss, sk]; AppendTo[ys, yk]; AppendTo[rhos, 1./sy]; gamma = sy/(yk . yk);
      If[Length[ss] > mem, ss = Rest[ss]; ys = Rest[ys]; rhos = Rest[rhos]]];
    it++;
    If[f - fn <= 10.^-15 Max[Abs[f], Abs[fn], 1.], x = xn; f = fn; g = gn; Break[]];
    x = xn; f = fn; g = gn];
  {x, f, it}];

(* the zero field of calibration mode, and one optimisation stage at smoothing tau *)
ndZeroField[spec_, trange_, kk_] := <|"TRange" -> trange, "K" -> kk, "Affine" -> False,
  "S0" -> Table[ConstantArray[0., faceShape[spec, i]], {i, spec["n"]}],
  "C" -> Table[ConstantArray[0., Append[faceShape[spec, i], kk + 2]], {i, spec["n"]}]|>;

ndCalStage[spec_, st_, field_, tau_, iters_] := Module[{x, res, cm},
  x = Flatten[Flatten /@ field["C"]];
  res = Quiet[ndLBFGS[ndCalValueGradient[st, #, tau] &, x, iters], {General::munfl}];
  cm = ndCalUnpack[st, res[[1]]];
  {Append[field, "C" -> Table[Developer`ToPackedArray[ArrayReshape[cm[[i]], Append[faceShape[spec, i], st["K"] + 2]]], {i, spec["n"]}]], -res[[2]], res[[3]]}];

(* exact embeddings of a field into finer spaces.  refine in t: the same spline with 2K knot intervals (Chaikin
   subdivision of the coefficients).  prolong in x: the same Raviart-Thomas field on the grid with twice as many
   cells in every direction (copy on the old faces, mean of the two neighbours on the new faces across the normal
   direction, repeat in the tangential directions). *)
ndRefineT[f_] := <|"TRange" -> f["TRange"], "K" -> 2 f["K"], "Affine" -> f["Affine"], "S0" -> f["S0"],
  "C" -> (Function[c, Module[{dims = Dimensions[c], m, a, b},
      m = ArrayReshape[c, {Times @@ Most[dims], Last[dims]}];
      a = (3 m[[All, ;; -2]] + m[[All, 2 ;;]])/4; b = (m[[All, ;; -2]] + 3 m[[All, 2 ;;]])/4;
      Developer`ToPackedArray[ArrayReshape[Transpose[Riffle[Transpose[a], Transpose[b]]], Append[Most[dims], 2 Last[dims] - 2]]]]] /@ f["C"])|>;

ndProlongX[f_, specFine_] := Module[{n = specFine["n"], up, zero},
  up[a_, i_] := Module[{x = a},
    Do[x = If[j == i, Map[Riffle[#, (Most[#] + Rest[#])/2] &, x, {j - 1}], Map[Riffle[#, #] &, x, {j - 1}]], {j, n}];
    Developer`ToPackedArray[N[x]]];
  (* natural faces carry no flux (multiplied by zero, not set to zero: the coefficient arrays have one more dimension) *)
  zero[a_, i_] := Module[{x = a},
    If[! specFine["Dirichlet"][[i, 1]], x[[Sequence @@ ReplacePart[ConstantArray[All, n], i -> 1]]] *= 0.];
    If[! specFine["Dirichlet"][[i, 2]], x[[Sequence @@ ReplacePart[ConstantArray[All, n], i -> -1]]] *= 0.];
    x];
  <|"TRange" -> f["TRange"], "K" -> f["K"], "Affine" -> f["Affine"],
    "S0" -> Table[zero[up[f["S0"][[i]], i], i], {i, n}], "C" -> If[f["K"] == 0, f["C"], Table[zero[up[f["C"][[i]], i], i], {i, n}]]|>];

(* coarse-to-fine schedule ending at the given grid and number of knot intervals: the grid is halved up to three
   times (while every direction stays even and the longest has at least 8 cells), K is halved down to 8, and the
   two are then doubled alternately *)
ndAutoLevels[cells_List, kfin_Integer] := Module[{p = 0, q = 0, c = cells, kk = kfin, seq},
  While[p < 3 && AllTrue[c, EvenQ] && Max[c]/2 >= 8, c = c/2; p++];
  While[EvenQ[kk] && kk/2 >= 8, kk = kk/2; q++];
  seq = {{c, kk}};
  While[p > 0 || q > 0,
    If[p > 0, c = 2 c; p--; AppendTo[seq, {c, kk}]];
    If[q > 0, kk = 2 kk; q--; AppendTo[seq, {c, kk}]]];
  seq];

(* the t-range: the hull of the Dirichlet data and of the critical points of W, rounded outward to a binary
   fraction of the span.  (W decreases to the left of its first critical point and increases to the right of its
   last one when it is bounded below, so truncating a competitor to this range does not increase its energy.)
   Heuristic: any range gives a valid bound. *)
ndAutoTRange[spec_] := Module[{t, roots, vals, lo, hi, span, qn},
  roots = Quiet@Check[t /. Solve[polyAt[spec["W1"], t] == 0, t, Reals], {}];
  roots = Select[Flatten[{roots}], NumericQ[#] && Im[N[#]] == 0 &];
  vals = Join[SetPrecision[Pick[Flatten[spec["G"]], Flatten[spec["Fixed"]], 1], Infinity], roots];
  If[vals === {}, Return[{-1, 1}]];
  lo = Min[vals]; hi = Max[vals]; span = N[hi - lo];
  If[! (span > 0.), span = 1.; lo = lo - 1/2; hi = hi + 1/2];
  qn = 2^Floor[Log2[span/32]];
  {Floor[lo/qn] qn, Ceiling[hi/qn] qn}];

(* calibration mode, coarse to fine.  levels: {{cells, K}, ...}; consecutive levels double the cells in every
   direction and/or double K.  At every level the field with the best CERTIFIED value so far is kept, so the
   reported bound never decreases along a level.  Returns {final spec, field, rows, certified value of the field}. *)
$ndSchedule = <|"First" -> {{1.*^-2, 300}, {1.*^-3, 400}, {1.*^-4, 600}}, "Later" -> {{1.*^-4, 300}, {3.*^-5, 400}},
  "Last" -> {{1.*^-4, 300}, {3.*^-5, 600}, {1.*^-5, 600}}|>;
Options[ndCalibrate] = {"Samples" -> 4, "Subintervals" -> 16, "Scale" -> 1., "Schedule" -> Automatic, "Monitor" -> None};
ndCalibrate[params_, levels_List, trange_, OptionsPattern[]] := Module[
  {qo = OptionValue["Samples"], qs = OptionValue["Subintervals"], scale = N[OptionValue["Scale"]], sched = OptionValue["Schedule"],
   monitor = OptionValue["Monitor"], cellsFine = levels[[-1, 1]], spec, prev = None, field = None, best, bestLB, lb, st, stages, rows = {}, t0 = AbsoluteTime[], res},
  sched = If[AssociationQ[sched], Join[$ndSchedule, sched], $ndSchedule];
  Do[
    With[{cells = levels[[li, 1]], kk = levels[[li, 2]]},
      spec = scalarFieldNDSpec[ndParamsOnGrid[params, cells, cellsFine]];
      If[field === None,
        field = ndZeroField[spec, trange, kk],
        If[cells =!= prev["Cells"], field = ndProlongX[field, spec]];
        If[kk =!= field["K"], field = ndRefineT[field]]];
      best = field; bestLB = ndLowerBound[spec, field, "Subintervals" -> qs];
      st = ndCalSetup[spec, trange, kk, qo];
      stages = Which[li == 1, sched["First"], li == Length[levels], sched["Last"], True, sched["Later"]];
      Do[
        res = ndCalStage[spec, st, field, stage[[1]] scale, stage[[2]]];
        field = res[[1]];
        lb = ndLowerBound[spec, field, "Subintervals" -> qs];
        If[lb > bestLB, bestLB = lb; best = field];
        If[monitor =!= None, monitor[<|"Level" -> li, "Levels" -> Length[levels], "Cells" -> cells, "K" -> kk, "Tau" -> stage[[1]] scale,
          "Iterations" -> res[[3]], "Smoothed" -> res[[2]], "Certified" -> N[lb], "Best" -> N[bestLB], "Seconds" -> AbsoluteTime[] - t0|>]],
        {stage, stages}];
      field = best; prev = spec;
      AppendTo[rows, <|"Cells" -> cells, "K" -> kk, "LowerBound" -> N[bestLB], "Seconds" -> AbsoluteTime[] - t0|>]],
    {li, Length[levels]}];
  {spec, field, rows, bestLB}];

(* ------------------------------------------------------------------------------------------------ *)
(* fields and certificates written by the Python reference                                            *)
(* ------------------------------------------------------------------------------------------------ *)

(* an exact number read from a JSON file: an integer, or a string "p", "p/q" or a decimal.  The string is checked
   before it is converted, so that a file received from someone else cannot make the kernel evaluate code. *)
ndRational[x_] := Module[{v = Which[
      IntegerQ[x], x,
      MatchQ[x, _Real], exactNumber[x],
      StringQ[x] && StringMatchQ[x, RegularExpression["-?[0-9]+(/[1-9][0-9]*)?"]], ToExpression[x],
      StringQ[x] && StringMatchQ[x, RegularExpression["-?[0-9]+\\.[0-9]+"]],
        With[{p = StringSplit[StringDelete[x, "-"], "."]},
          If[StringStartsQ[x, "-"], -1, 1] (FromDigits[p[[1]]] + FromDigits[p[[2]]]/10^StringLength[p[[2]]])],
      True, $Failed]},
  If[MatchQ[v, _Integer | _Rational], v,
    Throw[failure["MissingData", "Expected an integer, a fraction p/q or a decimal number."], "CertiField"]]];

(* a field directory written by export_field_dir: field.json and raw little-endian doubles, last index fastest *)
ndImportField[spec_, dir_String] := Module[{meta, n = spec["n"], kk, rd, s0, cc},
  meta = Quiet@Import[FileNameJoin[{dir, "field.json"}], "RawJSON"];
  If[! AssociationQ[meta], Throw[failure["MissingData", "No field.json in " <> dir], "CertiField"]];
  kk = meta["K"];
  If[! (IntegerQ[kk] && kk >= 0) || ! MatchQ[meta["TRange"], {_, _}],
    Throw[failure["MissingData", "field.json needs K (a non-negative integer) and TRange (two numbers)."], "CertiField"]];
  rd[name_, dims_] := Module[{v = BinaryReadList[FileNameJoin[{dir, name}], "Real64", ByteOrdering -> -1]},
    If[Length[v] =!= Times @@ dims, Throw[failure["MissingData", "Wrong size: " <> name], "CertiField"]];
    ArrayReshape[v, dims]];
  s0 = Table[rd["S0_" <> ToString[i - 1] <> ".f64", faceShape[spec, i]], {i, n}];
  cc = If[kk == 0, ConstantArray[{}, n], Table[rd["C" <> ToString[i - 1] <> ".f64", Append[faceShape[spec, i], kk + 2]], {i, n}]];
  <|"TRange" -> ndRational /@ meta["TRange"], "K" -> kk, "Affine" -> TrueQ[meta["Affine"]], "S0" -> s0, "C" -> cc|>];

(* nodal values of a candidate written by the reference (U.f64), and its boundary data (gval.f64) *)
ndImportArray[spec_, dir_String, name_String] := Module[{dims = spec["Cells"] + 1, v},
  v = BinaryReadList[FileNameJoin[{dir, name}], "Real64", ByteOrdering -> -1];
  If[Length[v] =!= Times @@ dims, Throw[failure["MissingData", "Wrong size: " <> name], "CertiField"]];
  Developer`ToPackedArray[ArrayReshape[v, dims]]];

(* one array of a NumPy .npy file (version 1 or 2, little-endian doubles, C order) *)
ndReadNPY[file_String] := Module[{s = OpenRead[file, BinaryFormat -> True], magic, ver, hlen, header, shape, data},
  magic = BinaryReadList[s, "Byte", 6]; ver = BinaryReadList[s, "Byte", 2];
  If[magic =!= {147, 78, 85, 77, 80, 89} || ! MatchQ[ver, {1 | 2, _}], Close[s]; Throw[failure["MissingData", "Not a NumPy array file: " <> file], "CertiField"]];
  hlen = BinaryRead[s, If[ver[[1]] == 1, "UnsignedInteger16", "UnsignedInteger32"], ByteOrdering -> -1];
  header = FromCharacterCode[BinaryReadList[s, "Byte", hlen]];
  data = BinaryReadList[s, "Real64", ByteOrdering -> -1]; Close[s];
  shape = FromDigits /@ StringCases[First[StringCases[header, "'shape': (" ~~ sh : Shortest[___] ~~ ")" :> sh], ""], DigitCharacter ..];
  If[! StringContainsQ[header, "'descr': '<f8'"] || ! StringContainsQ[header, "'fortran_order': False"] || Length[data] =!= Times @@ shape,
    Throw[failure["MissingData", "Unsupported NumPy array (need little-endian float64 in C order): " <> file], "CertiField"]];
  Developer`ToPackedArray[ArrayReshape[data, shape]]];

(* a certificate stored by the reference (save_certificate): stem.json (the problem in exact form, the claimed
   bounds as exact decimals) and stem.npz (boundary data, candidate and field as doubles).
   Returns the model parameters, the candidate, the field and the two claims. *)
ndImportCertificate[stem_String] := Module[{meta, tmp, files, z, pr, n, W, params, spec, u, field, kk},
  meta = Quiet@Import[stem <> ".json", "RawJSON"];
  If[! (AssociationQ[meta] && meta["format"] === "CertiField-ND-certificate-1" && FileExistsQ[stem <> ".npz"]),
    Throw[failure["MissingData", "No reference certificate at " <> stem], "CertiField"]];
  tmp = CreateDirectory[];
  z = Catch[
    files = ExtractArchive[stem <> ".npz", tmp];
    Association[(FileBaseName[#] -> ndReadNPY[#]) & /@ files], "CertiField"];
  Quiet@DeleteDirectory[tmp, DeleteContents -> True];
  If[! AssociationQ[z], Throw[z, "CertiField"]];
  pr = meta["problem"]; n = Length[pr["N"]]; kk = meta["K"];
  W = ndRational /@ pr["W"];
  (* Function @@ {...}: injecting the polynomial with With would rename the bound variable *)
  params = <|"Potential" -> Function @@ {\[FormalU], W . \[FormalU]^Range[0, Length[W] - 1]},
    "Stiffness" -> ndRational[pr["k"]], "Domain" -> Map[ndRational, pr["box"], {2}], "Cells" -> pr["N"],
    "BoundaryData" -> z["gval"],
    "NaturalFaces" -> Flatten[Table[If[pr["bc"][[i, s]] === "N", {{i, If[s == 1, "Lower", "Upper"]}}, {}], {i, n}, {s, 2}], 2]|>;
  spec = scalarFieldNDSpec[params];
  u = ndRealArray[z["U"], spec["Cells"] + 1];
  field = ndPackField[spec, <|"TRange" -> {ndRational[meta["tlo"]], ndRational[meta["thi"]]}, "K" -> kk, "Affine" -> TrueQ[meta["affine"]],
    "S0" -> Table[z["S0_" <> ToString[i - 1]], {i, n}], "C" -> If[kk == 0, ConstantArray[{}, n], Table[z["C" <> ToString[i - 1]], {i, n}]]|>];
  If[u === $Failed || field === $Failed, Throw[failure["MissingData", "The arrays of the certificate do not fit its grid: " <> stem], "CertiField"]];
  <|"Parameters" -> params, "Nodal" -> u, "Field" -> field, "ClaimedLower" -> ndRational[meta["lower"]], "ClaimedUpper" -> ndRational[meta["upper"]],
    "Name" -> Lookup[pr, "name", ""], "Info" -> Lookup[meta, "info", <||>]|>];

(* recompute both bounds of a reference certificate; True if they are at least as strong as its claims *)
ndVerifyReferenceCertificate[stem_String] := Module[{c = ndImportCertificate[stem], spec, lb, ub},
  spec = scalarFieldNDSpec[c["Parameters"]];
  lb = ndLowerBound[spec, c["Field"]];
  ub = ndUpperBound[spec, c["Nodal"]][[2]];
  <|"Verified" -> TrueQ[lb >= c["ClaimedLower"] && ub <= c["ClaimedUpper"] && c["ClaimedLower"] <= c["ClaimedUpper"]],
    "LowerBound" -> lb, "UpperBound" -> ub, "ClaimedLower" -> c["ClaimedLower"], "ClaimedUpper" -> c["ClaimedUpper"]|>];

(* ------------------------------------------------------------------------------------------------ *)
(* the certificate                                                                                   *)
(* ------------------------------------------------------------------------------------------------ *)

$ndCalibrationOptions = {"TRange" -> Automatic, "KnotIntervals" -> Automatic, "Levels" -> Automatic, "Scale" -> Automatic,
  "Samples" -> 4, "Subintervals" -> 16, "Schedule" -> Automatic};

(* Method -> Automatic | "Convex" | "Calibration" | {"Calibration", suboptions} *)
ndMethod[m_] := Module[{name, sub},
  {name, sub} = Replace[m, {Automatic -> {Automatic, {}}, s_String :> {s, {}}, {s_String, o___?OptionQ} :> {s, Flatten[{o}]},
    _ :> Throw[failure["BadOption", "Method must be Automatic, \"Convex\", \"Calibration\" or {\"Calibration\", suboptions}."], "CertiField"]}];
  If[! MemberQ[{Automatic, "Convex", "Calibration"}, name],
    Throw[failure["BadOption", "Method must be Automatic, \"Convex\", \"Calibration\" or {\"Calibration\", suboptions}."], "CertiField"]];
  With[{bad = Complement[ToString /@ Keys[sub], ToString /@ Keys[$ndCalibrationOptions]]},
    If[bad =!= {} || (name =!= "Calibration" && sub =!= {}),
      Throw[failure["BadOption", "Unknown suboption of Method: " <> ToString[If[bad === {}, Keys[sub], bad]] <> ". Calibration mode takes " <> ToString[Keys[$ndCalibrationOptions]] <> "."], "CertiField"]]];
  {name, Association[$ndCalibrationOptions, Association[(ToString[#[[1]]] -> #[[2]]) & /@ sub]]}];

ndLevelsQ[levels_, n_] := MatchQ[levels, {{{__Integer}, _Integer} ..}] && AllTrue[levels, Length[#[[1]]] == n && Min[#[[1]]] >= 1 && #[[2]] >= 1 &] &&
  And @@ Table[(levels[[j + 1, 1]] === levels[[j, 1]] || levels[[j + 1, 1]] === 2 levels[[j, 1]]) &&
    (levels[[j + 1, 2]] == levels[[j, 2]] || levels[[j + 1, 2]] == 2 levels[[j, 2]]), {j, Length[levels] - 1}];

certifyScalarFieldND[params_Association, cand_, opts_List] := Module[
  {spec, t0 = AbsoluteTime[], mname, sub, verbose, refine, gfun, given, mode, u, ub, p, fields = {}, lbs, best, lo, hi, margin, qs = 16,
   trange, kfin, levels, scale, cal, rows = {}, wv, known = <||>},
  If[nd2ParamsQ[params], Return[certifyScalarFieldND2[params, cand, opts]]];       (* Kernel/ScalarFieldND2.wl *)
  spec = scalarFieldNDSpec[KeyDrop[params, "Field"]];
  {mname, sub} = ndMethod[OptionValue[CertifiedVariationalBound, opts, Method]];
  verbose = TrueQ[OptionValue[CertifiedVariationalBound, opts, "Verbose"]];
  refine = TrueQ[OptionValue[CertifiedVariationalBound, opts, "RefineCandidate"]];
  given = Lookup[params, "Field", None];
  If[StringQ[given], given = ndImportField[spec, given]];
  If[given =!= None,
    given = ndPackField[spec, given];
    If[given === $Failed, Throw[failure["BadField", "\"Field\" must be a calibration field on the grid of the model, or a directory written by the Python reference."], "CertiField"]]];
  (* convex mode for convex potentials, calibration mode otherwise; a supplied field replaces the optimiser *)
  mode = Which[mname =!= Automatic, mname, AssociationQ[given], "Given", spec["Superlinear"] && spec["Convex"], "Convex", True, "Calibration"];
  gfun = With[{g = Lookup[params, "BoundaryData", 0]}, If[NumericQ[g] || ArrayQ[g, spec["n"]], None, g]];
  If[verbose, Print["  ScalarFieldND: ", ndGridString[spec["Cells"]], " cells, ", If[mode === "Given", "supplied field", ToLowerCase[mode] <> " mode"]]];
  u = ndCompetitor[spec, cand, refine, gfun];
  ub = ndUpperBound[spec, u];
  If[verbose, Print["  candidate: E <= ", ndNum[ub[[2]]], "   (", ndSec[AbsoluteTime[] - t0], " s)"]];
  If[spec["Superlinear"] && mode =!= "Calibration",
    p = ndFVSolve[spec, Mean[vertexPart[u, spec["n"], #] & /@ Tuples[{0, 1}, spec["n"]]]];
    AppendTo[fields, ndConvexField[spec, p]]];
  If[AssociationQ[given], AppendTo[fields, given]];
  If[mode === "Calibration",
    qs = sub["Subintervals"];
    trange = Replace[sub["TRange"], Automatic :> ndAutoTRange[spec]];
    trange = exactify[trange];
    If[! (MatchQ[trange, {_Integer | _Rational, _Integer | _Rational}] && trange[[1]] < trange[[2]]),
      Throw[failure["BadOption", "\"TRange\" must be {tlo, thi} with rational tlo < thi."], "CertiField"]];
    kfin = Replace[sub["KnotIntervals"], Automatic -> 16];
    levels = Replace[sub["Levels"], Automatic :> If[IntegerQ[kfin] && kfin >= 1, ndAutoLevels[spec["Cells"], kfin], $Failed]];
    If[! ndLevelsQ[levels, spec["n"]] || levels[[-1, 1]] =!= spec["Cells"],
      Throw[failure["BadOption", "\"Levels\" must be a list of {cells, K} that doubles the cells and/or K from one level to the next and ends on the grid of the model; \"KnotIntervals\" must be a positive integer."], "CertiField"]];
    wv = polyAt[N[spec["W"]], Subdivide[N[trange[[1]]], N[trange[[2]]], 256]];
    scale = Replace[sub["Scale"], Automatic :> Max[1., Max[wv] - Min[wv]]];
    If[! (NumericQ[scale] && scale > 0), Throw[failure["BadOption", "\"Scale\" must be a positive number."], "CertiField"]];
    If[verbose, Print["  t-range ", trange, ", scale ", scale, ", levels ", StringRiffle[(ndGridString[#[[1]]] <> " K" <> ToString[#[[2]]]) & /@ levels, " -> "]]];
    cal = ndCalibrate[KeyDrop[params, "Field"], levels, trange, "Samples" -> sub["Samples"], "Subintervals" -> qs, "Scale" -> scale,
      "Schedule" -> sub["Schedule"],
      "Monitor" -> If[verbose, Function[s, Print["  level ", s["Level"], "/", s["Levels"], "  ", ndGridString[s["Cells"]], "  K ", s["K"], "  tau ", ToString[s["Tau"], InputForm],
        "  iterations ", s["Iterations"], "  certified lower bound ", ndNum[s["Best"]], "   (", ndSec[s["Seconds"]], " s)"]], None]];
    rows = cal[[3]];
    AppendTo[fields, cal[[2]]];
    known[Length[fields]] = cal[[4]]];                        (* certified on this grid with the same subintervals *)
  If[fields === {}, Throw[failure["Inconclusive", "No lower bound is available: the potential has no superlinear growth and no field was supplied."], "CertiField"]];
  lbs = Table[Lookup[known, j, ndLowerBound[spec, fields[[j]], "Subintervals" -> qs]], {j, Length[fields]}];
  best = First[Ordering[N[lbs], -1]];
  (* The claims are the computed bounds weakened by one part in 10^9 of the larger modulus.  A later verification
     recomputes both bounds, and the machine root finder that places the conjugate tables (conjArgmax) need not
     return the same last bits on another machine.  Any placement gives a valid bound; only its last digits move
     (the Python reference does the same: claim_strings). *)
  lo = exactNumber[lbs[[best]]]; hi = exactNumber[ub[[2]]];
  margin = Max[Abs[lo], Abs[hi]]/10^9;
  <|"LowerBound" -> lo - margin, "UpperBound" -> hi + margin, "Seconds" -> AbsoluteTime[] - t0, "Field" -> fields[[best]],
    "Nodal" -> u, "BoundaryNodalValues" -> spec["G"], "Subintervals" -> qs,
    "Mode" -> If[fields[[best]]["K"] == 0, "Convex", "Calibration"], "Levels" -> rows, "Spec" -> spec|>];

(* plain-text helpers for progress lines and descriptions (ASCII only: they are also printed in terminals) *)
ndGridString[cells_List] := StringRiffle[ToString /@ cells, "x"];
ndNum[x_] := ToString[NumberForm[N[x], 10]];
ndSec[x_] := ToString[NumberForm[N[x], {Infinity, 1}]];

scalarNDDescription[params_, res_] := Module[{sp = res["Spec"], d = Max[displayDigits[res], 8], nat, grid},
  If[TrueQ[sp["ND2"]], Return[nd2Description[params, res]]];
  nat =Flatten[Table[If[sp["Dirichlet"][[i, s]], Nothing, {i, If[s == 1, "Lower", "Upper"]}], {i, sp["n"]}, {s, 2}], 1];
  grid = ndGridString[sp["Cells"]];
  <|"Units" -> "energy units of the model",
    "Theorem" -> Row[{"Every u in H^1 of the box ", exactString[sp["Domain"]], " whose trace on the Dirichlet faces is the piecewise multilinear interpolant of the nodal boundary values on the ",
      grid, " grid has E[u] \[GreaterEqual] ", lowerString[res["LowerBound"], d], "; the competitor has E \[LessEqual] ", upperString[res["UpperBound"], d], "."}],
    "AdmissibleClass" -> Row[{"all u in H^1 of the box ", exactString[sp["Domain"]], " whose trace on the Dirichlet faces is the piecewise multilinear interpolant of the nodal boundary values (machine numbers, stored in the certificate) on the ",
      grid, " grid", If[nat === {}, "", Row[{"; natural (free) faces: ", nat}]]}],
    "LowerBoundMethod" -> If[res["Mode"] === "Convex",
      "complementary-energy bound with a Raviart\[Dash]Thomas flux from a cell-centred finite-volume solve; outward-rounded arithmetic",
      Row[{"calibration field (Raviart\[Dash]Thomas in x, quadratic B-splines in u, ", res["Field"]["K"], " knot intervals on ", exactString[res["Field"]["TRange"]], "); outward-rounded arithmetic"}]],
    "UpperBoundMethod" -> "energy of the piecewise multilinear function through the candidate's nodal values (exact gradient term, midpoint rule with remainder)",
    "TrustedBase" -> "IEEE-754 machine arithmetic (+, \[Minus], \[Times], \[Divide]), Wolfram Language exact rational arithmetic, and docs/SCALAR_ND.md (Theorem 1, Propositions 2, 6, 7, Lemmas 3\[Dash]5; proofs checked by one author only)",
    "NotClaimed" -> "that the candidate is the minimiser or that the minimiser is unique; the gap bounds how far below the candidate any admissible field can be",
    "References" -> "docs/SCALAR_ND.md; calibrations and duality: Alberti\[Dash]Bouchitt\[EAcute]\[Dash]Dal Maso (2003), Bouchitt\[EAcute]\[Dash]Fragal\[AGrave] (2018); complementary energy: Prager\[Dash]Synge (1947)."|>];

(* re-verification: both bounds recomputed from the stored nodal values and the stored field.  The Dirichlet datum
   is the one stored in the certificate (exact machine numbers); it must agree with the boundary data of the model
   to rounding, because a function of the coordinates need not evaluate to the same last bits on another machine. *)
verifyScalarFieldND[a_Association] := Module[{c = a["Certificate"], spec, g, lb, ub},
  If[nd2ParamsQ[a["Parameters"]], Return[verifyScalarFieldND2[a]]];
  spec =scalarFieldNDSpec[KeyDrop[a["Parameters"], "Field"]];
  g = ndRealArray[Lookup[c, "BoundaryNodalValues", $Failed], spec["Cells"] + 1];
  If[g === $Failed || Max[Abs[(1 - spec["Fixed"]) g]] != 0. || Max[Abs[g - spec["G"]]] > 10.^-12 (1. + Max[Abs[spec["G"]]]), Return[False]];
  spec = Append[spec, "G" -> g];
  lb = ndLowerBound[spec, c["Field"], "Subintervals" -> Lookup[c, "Subintervals", 16]];
  ub = ndUpperBound[spec, c["Nodal"]][[2]];
  NumericQ[lb] && NumericQ[ub] && lb >= exactNumber[a["LowerBound"]] && ub <= exactNumber[a["UpperBound"]]];

End[];
