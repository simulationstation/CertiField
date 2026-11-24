(* ::Package:: *)

(* ::Title:: *)
(*CertiField: stiffness of an elastic block with a spherical cavity*)

(* ::Text:: *)
(*Wolfram Language port of Reference/Curved3D/curved3d.py; the argument, the results and the run commands of the*)
(*reference are in Reference/Curved3D/NOTES.md.*)
(**)
(*Body.  The block (-1,1)^3 minus the closed ball |x| <= 1/2 (or the solid block), isotropic, Young's modulus 1,*)
(*Poisson ratio 3/10.  The faces y = +-1, z = +-1 and the sphere are free of traction; the faces x = +-1 are rigid*)
(*grips moved apart by delta, clamped (u = (+-delta/2, 0, 0)) or frictionless (u_x = +-delta/2).  The stiffness K is*)
(*defined by min E[u] = K delta^2/2, E[u] = (1/2) Integrate[eps(u):C:eps(u)].  One octant (0,1)^3 is computed; the*)
(*planes x = 0, y = 0, z = 0 are mirror planes.*)
(**)
(*Upper bound.  K <= 8 Integrate[eps(v):C:eps(v)] over Omega+, v = (x/2) e_x + a C^1 quadratic B-spline whose*)
(*coefficients satisfy the mirror and grip conditions exactly (unit grip separation).*)
(*Lower bound.  K >= A^2/Q, A = 4, Q = 8 Integrate[tau:C^-1:tau] over Omega+, tau = tau_p + inc(rho^2 Psi): inc the*)
(*Beltrami operator, rho = eta(x) + y^2 + z^2 - 1/4 (zero on the sphere), Psi a tensor of C^1 quadratic B-splines*)
(*whose coefficients satisfy exact ties, tau_p an axial stress outside the column |y|, |z| < 1/2 that carries the force.*)
(*For every such Psi, div tau = 0 and the traction vanishes on the sphere and on the free faces, the shear traction on*)
(*the mirror planes and on frictionless grips.*)
(*Omega+ is a set of octree boxes that contains the octant of the body (exact integer classification).  On every box*)
(*the integrands are polynomials of degree at most 12 (stress) and 4 (strain) in each variable, so the 7- and 3-point*)
(*Gauss rules are exact; nodes and weights are enclosed in exact rational intervals, and every value is computed in*)
(*the outward-rounded arithmetic of RigorousArithmetic.wl with the overflow checks of ScalarFieldND.wl (ndChk).*)
(**)
(*Rigorous: c3DualUpper, c3PrimalUpper, with c3Tables, c3Boxes and c3TiesHold.  Heuristic, never trusted: the*)
(*machine-arithmetic assembly and solves that find Psi and v (c3SolveDual, c3SolvePrimal).*)

Begin["CertiField`Private`"];

(* ------------------------------------------------------------------------------------------------ *)
(* index conventions (1-based): components xx yy zz xy yz xz; derivative variants of the splines      *)
(* ------------------------------------------------------------------------------------------------ *)

$c3Comps = {{1, 1}, {2, 2}, {3, 3}, {1, 2}, {2, 3}, {1, 3}};
c3Cidx[i_, j_] := First[FirstPosition[$c3Comps, Sort[{i, j}]]];
$c3Var = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0}, {0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {0, 1, 1}, {1, 0, 1}};
$c3First = Table[First[FirstPosition[$c3Var, UnitVector[3, k]]], {k, 3}];
$c3Second = Table[First[FirstPosition[$c3Var, UnitVector[3, c[[1]]] + UnitVector[3, c[[2]]]]], {c, $c3Comps}];

(* Beltrami operator tau_ij = e_ikl e_jmn d_k d_m Phi_ln: tau_s = Sum G[s, p, c] (d_c Phi_p), c a pair of derivatives *)
$c3G = Module[{g = ConstantArray[0, {6, 6, 6}]},
  Do[With[{i = $c3Comps[[s, 1]], j = $c3Comps[[s, 2]]},
      Do[With[{e = Signature[{i, k, l}] Signature[{j, m, n}]},
          If[e != 0, g[[s, c3Cidx[l, n], c3Cidx[k, m]]] += e]],
        {k, 3}, {l, 3}, {m, 3}, {n, 3}]],
    {s, 6}];
  g];
(* its non-zero entries {s, p, c, G[s, p, c]}, in the order of the reference's loops *)
$c3Triples = Append[#, $c3G[[Sequence @@ #]]] & /@ Position[$c3G, x_Integer /; x != 0, {3}];

(* engineering strain e_s = Sum GS[s, a, k] d_k v_a *)
$c3GS = Module[{g = ConstantArray[0, {6, 3, 3}]},
  g[[1, 1, 1]] = 1; g[[2, 2, 2]] = 1; g[[3, 3, 3]] = 1;
  g[[4, 1, 2]] = 1; g[[4, 2, 1]] = 1; g[[5, 2, 3]] = 1; g[[5, 3, 2]] = 1; g[[6, 1, 3]] = 1; g[[6, 3, 1]] = 1;
  g];

c3Compliance[nu_] := (1 + nu) DiagonalMatrix[{1, 1, 1, 2, 2, 2}] - nu Table[Boole[i <= 3 && j <= 3], {i, 6}, {j, 6}];
c3Stiffness[lam_, mu_] := DiagonalMatrix[{2 mu, 2 mu, 2 mu, mu, mu, mu}] + lam Table[Boole[i <= 3 && j <= 3], {i, 6}, {j, 6}];

(* ------------------------------------------------------------------------------------------------ *)
(* exact Gauss rules and one-dimensional tables                                                       *)
(* ------------------------------------------------------------------------------------------------ *)

c3RiMul[{a_, b_}, {c_, d_}] := With[{p = {a c, a d, b c, b d}}, {Min[p], Max[p]}];
c3RiHorner[coefs_List, x_] := Fold[With[{m = c3RiMul[#1, x]}, {m[[1]] + #2, m[[2]] + #2}] &, {Last[coefs], Last[coefs]}, Reverse[Most[coefs]]];

(* Gauss-Legendre rule on [0, 1] with n points: exact rational enclosures {lo, hi} of the nodes and of the weights
   (bisection of the Legendre polynomial in exact arithmetic; weights 1/((1 - x^2) P'(x)^2) in rational interval
   arithmetic), as gauss_enclosures of the reference *)
c3GaussEnclosures[n_Integer] := c3GaussEnclosures[n] = Module[{p, dp, m = 1024, xs, vs, roots = {}, lo, hi, neg, mid, fm},
  p = CoefficientList[LegendreP[n, \[FormalX]], \[FormalX]];
  dp = polyD[p];
  xs = Table[-1 + 2 k/m, {k, 0, m}];
  vs = polyAt[p, #] & /@ xs;
  Do[
    Which[
      vs[[k]] == 0, AppendTo[roots, {xs[[k]], xs[[k]]}],
      vs[[k + 1]] != 0 && Sign[vs[[k]]] != Sign[vs[[k + 1]]],
        lo = xs[[k]]; hi = xs[[k + 1]]; neg = vs[[k]] < 0;
        Do[mid = (lo + hi)/2; fm = polyAt[p, mid];
          If[fm == 0, lo = hi = mid; Break[]];
          If[(fm < 0) === neg, lo = mid, hi = mid], {120}];
        AppendTo[roots, {lo, hi}]],
    {k, m}];
  If[Length[roots] =!= n, Throw[failure["Internal", "Gauss rule: wrong number of roots."], "CertiField"]];
  Transpose[Map[Function[r, Module[{d = c3RiHorner[dp, r], xx = c3RiMul[r, r], den},
      den = c3RiMul[{1 - xx[[2]], 1 - xx[[1]]}, c3RiMul[d, d]];
      If[! (den[[1]] > 0), Throw[failure["Internal", "Gauss rule: the weight enclosure failed."], "CertiField"]];
      {{(1 + r[[1]])/2, (1 + r[[2]])/2}, {1/den[[2]], 1/den[[1]]}}]], roots]]];

(* the 1D functions as polynomials in the local coordinate t of cell e (x = (e + t)/n), ascending exact coefficients:
   one; the indicator of the column (cells e < kcol); x^2; and eta, x^2 flattened at the grip beyond kflat/n
   (x^2 - (x - b)^2/(1 - b), C^1 with eta'(1) = 0) *)
c3Funcs[n_, kcol_, kflat_] := <|
  "one" -> Function[e, {1}],
  "ind" -> Function[e, If[e < kcol, {1}, {0}]],
  "sq" -> Function[e, {e^2, 2 e, 1}/n^2],
  "sqf" -> Function[e, If[kflat === None || e < kflat, {e^2, 2 e, 1}/n^2,
    With[{g = n/(n - kflat)}, {e^2 - g (e - kflat)^2, 2 e - 2 g (e - kflat), 1 - g}/n^2]]]|>;

$c3Shapes = {{1/2, -1, 1/2}, {1/2, 1, -1}, {0, 0, 1/2}};          (* the quadratic B-spline on a cell, three pieces *)

c3MinQ[x_, y_] := (x + y - Abs[x - y])/2;                          (* exact elementwise min and max of rationals *)
c3MaxQ[x_, y_] := (x + y + Abs[x - y])/2;
(* exact enclosures {lo-list, hi-list} of a polynomial on the intervals [tl, th] (lists, tl >= 0): rational interval
   Horner, as ri_horner of the reference *)
c3HornerT[coefs_List, {tl_, th_}] := Fold[
  Function[{p, c}, {c3MinQ[p[[1]] tl, p[[1]] th] + c, c3MaxQ[p[[2]] tl, p[[2]] th] + c}],
  {Last[coefs] + 0 tl, Last[coefs] + 0 tl}, Reverse[Most[coefs]]];

(* largest (smallest) machine numbers <= (>=) the exact numbers of a list, checked by exact comparison; returned as a
   packed array (the entries are machine numbers by construction) *)
c3FLo[q_List] := Module[{x = Developer`ToPackedArray[N[q] + 0.], bad},
  bad = Flatten[Position[Sign[SetPrecision[x, Infinity] - q], 1]];
  While[bad =!= {},
    x[[bad]] = nextDown /@ x[[bad]];
    bad = Flatten[Position[Sign[SetPrecision[x, Infinity] - q], 1]]];
  If[! VectorQ[x, MachineNumberQ], Throw[failure["Internal", "Directed rounding failed."], "CertiField"]];
  Developer`ToPackedArray[x]];
c3FHi[q_List] := -c3FLo[-q];

c3Straddles[t_] := Max[(1 - UnitStep[t["Lo"]]) (1 - UnitStep[-t["Hi"]])] > 0;

(* 1D tables at the Gauss points of every dyadic sub-interval (level lev <= lmax, index j) of a cell, lj = 2^lev + j:
   "Bas" [lj, d, m, q]: d-th derivative of the B-spline piece m times n^d; "W" [lj, q]: weight/(n 2^lev);
   "Fn"[name] [e nlj + lj, d, q]: d-th derivative of the function name on cell e times n^d.  Each as "Lo"/"Hi", machine
   numbers enclosing the exact rational enclosure, and "Mid" (used only by the solver).  "Straddle": some entry has
   lo < 0 < hi (then the products below fall back to the general interval product). *)
c3Tables[n_Integer, lmax_Integer, nq_Integer, kflat_] := c3Tables[n, lmax, nq, kflat] = Module[
  {funcs = c3Funcs[n, n/2, kflat], nodes, weights, ljs, nlj, tl, th, wl, wh, derivs, pack, enc, bas, fn, w},
  {nodes, weights} = c3GaussEnclosures[nq];
  ljs = Flatten[Table[{lev, j}, {lev, 0, lmax}, {j, 0, 2^lev - 1}], 1];
  nlj = Length[ljs];
  tl = Flatten[Table[(lj[[2]] + nodes[[q, 1]])/2^lj[[1]], {lj, ljs}, {q, nq}]];
  th = Flatten[Table[(lj[[2]] + nodes[[q, 2]])/2^lj[[1]], {lj, ljs}, {q, nq}]];
  wl = Flatten[Table[weights[[q, 1]]/(n 2^lj[[1]]), {lj, ljs}, {q, nq}]];
  wh = Flatten[Table[weights[[q, 2]]/(n 2^lj[[1]]), {lj, ljs}, {q, nq}]];
  derivs[c_] := NestList[polyD, c, 2];
  pack[{lo_, hi_}, dims_] := <|"Lo" -> ArrayReshape[c3FLo[lo], dims], "Hi" -> ArrayReshape[c3FHi[hi], dims],
    "Mid" -> Developer`ToPackedArray[ArrayReshape[N[(lo + hi)/2], dims]]|>;
  w = pack[{wl, wh}, {nlj, nq}];
  enc = Table[c3HornerT[n^d derivs[$c3Shapes[[m]]][[d + 1]], {tl, th}], {d, 0, 2}, {m, 3}];
  bas = Map[Transpose[#, {2, 3, 1, 4}] &, pack[{Flatten[enc[[All, All, 1]]], Flatten[enc[[All, All, 2]]]}, {3, 3, nlj, nq}]];
  fn = Association[Table[name -> Module[{e2 = Table[c3HornerT[n^d derivs[funcs[name][e]][[d + 1]], {tl, th}], {e, 0, n - 1}, {d, 0, 2}]},
      Map[Flatten[Transpose[#, {1, 3, 2, 4}], 1] &, pack[{Flatten[e2[[All, All, 1]]], Flatten[e2[[All, All, 2]]]}, {n, 3, nlj, nq}]]],
    {name, Keys[funcs]}]];
  <|"nq" -> nq, "nlj" -> nlj, "Bas" -> bas, "W" -> w, "Fn" -> fn,
    "Straddle" -> (c3Straddles[bas] || c3Straddles[w] || AnyTrue[Values[fn], c3Straddles])|>];

(* ------------------------------------------------------------------------------------------------ *)
(* geometry: the boxes that cover the octant of the body                                               *)
(* ------------------------------------------------------------------------------------------------ *)

(* rows {ex, ey, ez, lev, jx, jy, jz, cut} (0-based): box j of level lev in cell e.  A box is dropped only if its
   farthest corner lies in the closed ball of radius^2 a2 (exact integer test); boxes in the solid are kept whole;
   boxes cut by the sphere are refined down to depth lmax and then kept whole (cut = 1).  Sorted by cell. *)
c3Boxes[n_Integer, lmax_Integer, a2_] := c3Boxes[n, lmax, a2] = Module[{cur, out = {}, sz, lo, rmin, rmax, thr, inside, void, cut, rows, base},
  cur = Flatten[Table[{ex, ey, ez, 0, 0, 0}, {ex, 0, n - 1}, {ey, 0, n - 1}, {ez, 0, n - 1}], 2];
  rows[sel_, lev_, flag_] := If[sel === {}, {},
    Join[sel[[All, 1 ;; 3]], ConstantArray[lev, {Length[sel], 1}], sel[[All, 4 ;; 6]], ConstantArray[flag, {Length[sel], 1}], 2]];
  Do[
    If[cur === {}, Break[]];
    sz = n 2^lev; thr = a2 sz^2;
    lo = cur[[All, 1 ;; 3]] 2^lev + cur[[All, 4 ;; 6]];
    rmin = Total[lo^2, {2}]; rmax = Total[(lo + 1)^2, {2}];
    inside = UnitStep[rmin - thr];                          (* nearest corner outside the open ball: in the solid *)
    void = (1 - inside) UnitStep[thr - rmax];               (* farthest corner in the closed ball: dropped *)
    cut = (1 - inside) (1 - void);
    out = Join[out, rows[Pick[cur, inside, 1], lev, 0]];
    If[lev == lmax,
      out = Join[out, rows[Pick[cur, cut, 1], lev, 1]]; cur = {},
      base = Pick[cur, cut, 1];
      cur = If[base === {}, {},
        base = Join[base[[All, 1 ;; 3]], 2 base[[All, 4 ;; 6]], 2];
        Join @@ Flatten[Table[(# + {0, 0, 0, dx, dy, dz}) & /@ base, {dx, 0, 1}, {dy, 0, 1}, {dz, 0, 1}], 2]]],
    {lev, 0, lmax}];
  Developer`ToPackedArray[Sort[out]]];

(* exact volumes of the boxes: {all, those proved to lie in the solid} *)
c3Volumes[boxes_, n_] := {Total[1/(n 2^boxes[[All, 4]])^3], Total[(1 - boxes[[All, 8]])/(n 2^boxes[[All, 4]])^3]};

(* ------------------------------------------------------------------------------------------------ *)
(* ties between spline coefficients (mirror planes, free faces, grips)                                 *)
(* ------------------------------------------------------------------------------------------------ *)

(* full coefficients 0..n+1 of one direction from the reduced ones: left mirror plane c0 = s c1; right "free",
   "zero2" (c_n = c_{n+1} = 0: value and normal derivative vanish) or a sign s' (c_{n+1} = s' c_n) *)
c3Expand1D[n_, s_, right_] := Module[{cols, pos, rules},
  cols = Which[right === "free", Range[n + 1], right === "zero2", Range[n - 1], True, Range[n]];
  pos = AssociationThread[cols -> Range[Length[cols]]];
  rules = Join[Table[{j + 1, pos[j]} -> 1, {j, cols}],
    If[KeyExistsQ[pos, 1], {{1, pos[1]} -> s}, {}],
    If[IntegerQ[right] && KeyExistsQ[pos, n], {{n + 2, pos[n]} -> right}, {}]];
  {SparseArray[rules, {n + 2, Length[cols]}], cols}];

(* the tie matrix P (full coefficients = P . reduced ones, every row has at most one entry +-1) and the position of the
   full coefficient that carries each reduced one *)
c3TieMatrix[n_, specs_] := Module[{m = n + 2, e, k, own = {}, rules = {}, off = 0},
  Do[
    e = c3Expand1D[n, #[[1]], #[[2]]] & /@ specs[[comp]];
    k = KroneckerProduct[e[[1, 1]], e[[2, 1]], e[[3, 1]]];
    rules = Join[rules, Map[({#[[1, 1]] + (comp - 1) m^3, #[[1, 2]] + off} -> #[[2]]) &, Most[ArrayRules[k]]]];
    own = Join[own, Flatten[Outer[(comp - 1) m^3 + (#1 m + #2) m + #3 + 1 &, e[[1, 2]], e[[2, 2]], e[[3, 2]]]]];
    off += Dimensions[k][[2]],
    {comp, Length[specs]}];
  {SparseArray[rules, {Length[specs] m^3, off}], own}];

(* stress potential: Psi_ab is even or odd at the mirror planes as the tensor demands; two layers vanish at the free
   faces y = 1, z = 1; at a frictionless grip x = 1 the grip is a mirror plane of the potential *)
c3PsiSpecs[clamped_] := Table[
  With[{s = (-1)^(Boole[c[[1]] == k] + Boole[c[[2]] == k])}, {s, If[k > 1, "zero2", If[clamped, "free", s]]}],
  {c, $c3Comps}, {k, 3}];
(* displacement: v_a odd in x_a at the mirror planes; the spline part vanishes at the grip (only v_x if frictionless) *)
c3VSpecs[clamped_] := Table[{If[a == k, -1, 1], If[k == 1, If[a == 1 || clamped, -1, "free"], "free"]}, {a, 3}, {k, 3}];
c3Ties[n_, kind_, clamped_] := c3Ties[n, kind, clamped] = c3TieMatrix[n, If[kind === "Psi", c3PsiSpecs, c3VSpecs][clamped]];

(* True if a coefficient array satisfies every tie exactly (it is reproduced from its reduced part) *)
c3TiesHold[{p_, own_}, c_] := With[{flat = Flatten[c]}, Max[Abs[p . flat[[own]] - flat]] == 0.];

c3RedMap[p_] := Module[{nz = p["NonzeroPositions"], vals = p["NonzeroValues"], red = ConstantArray[0, Length[p]], sg = ConstantArray[0, Length[p]]},
  red[[nz[[All, 1]]]] = nz[[All, 2]]; sg[[nz[[All, 1]]]] = Round[vals];
  {red, sg}];

(* ------------------------------------------------------------------------------------------------ *)
(* the problem                                                                                       *)
(* ------------------------------------------------------------------------------------------------ *)

c3Spec[params_Association] := Module[{grips, n, depths, cavity, clamped, nu = 3/10},
  grips = Lookup[params, "Grips", "Clamped"];
  If[! MemberQ[{"Clamped", "Frictionless"}, grips], Throw[failure["BadModel", "\"Grips\" must be \"Clamped\" or \"Frictionless\"."], "CertiField"]];
  n = Lookup[params, "Cells", 8];
  If[! (IntegerQ[n] && n >= 4 && Mod[n, 4] == 0),
    Throw[failure["BadModel", "\"Cells\" (spline cells per half-edge) must be a positive multiple of 4."], "CertiField"]];
  depths = Lookup[params, "Depths", {2, 3, 3, 5}];
  If[! (VectorQ[depths, IntegerQ] && Length[depths] == 4 && Min[depths] >= 0 && Max[depths] <= 8),
    Throw[failure["BadModel", "\"Depths\" must be four octree depths from 0 to 8: {stress solve, stress certificate, displacement solve, displacement certificate}."], "CertiField"]];
  cavity = Lookup[params, "Cavity", True];
  If[! BooleanQ[cavity], Throw[failure["BadModel", "\"Cavity\" must be True or False."], "CertiField"]];
  clamped = grips === "Clamped";
  (* the conditions of the construction (_check_admissible of the reference) hold for these choices: the column
     |y|, |z| < 1/2 ends on grid lines and contains the cavity; the weight is flattened at the grip beyond
     x = 3/4 >= 1/2; the solid block has no free surface with a normal that has an x-component *)
  <|"N" -> n, "M" -> n + 2, "Grips" -> grips, "Clamped" -> clamped, "Cavity" -> cavity, "a2" -> If[cavity, 1/4, 0],
    "nu" -> nu, "lam" -> nu/((1 + nu) (1 - 2 nu)), "mu" -> 1/(2 (1 + nu)), "kflat" -> If[clamped, None, 3 n/4],
    "Eta" -> {If[clamped, "sq", "sqf"], "sq", "sq"},
    (* the force-carrying stress {component, coefficient, fx, fy, fz}: 4/3 (1 - ind(y) ind(z)) e_x e_x, or e_x e_x *)
    "STerms" -> If[cavity, {{1, 4/3, "one", "one", "one"}, {1, -4/3, "one", "ind", "ind"}}, {{1, 1, "one", "one", "one"}}],
    "Area" -> 4, "Depths" -> depths|>];

(* ------------------------------------------------------------------------------------------------ *)
(* rigorous evaluation                                                                               *)
(* ------------------------------------------------------------------------------------------------ *)

(* Products of intervals.  Every result is rounded outward (rDn, rUp) and checked for machine overflow (ndChk).
   c3MulPt: exact machine numbers x (a packed array) times an interval {tl, th} whose arrays have leading dimensions
   of x.  c3MulTab: an interval times a table interval with no entry tl < 0 < th (the ends of the product are then
   products of ends chosen by signs).  c3MulNonneg: an interval {rl, rh} with rl >= 0 times an interval. *)
c3MulPt[x_, {tl_, th_}] := Module[{p1 = tl x, p2 = th x, m = UnitStep[x], om},
  om = 1 - m;
  ndChk[{rDn[m p1 + om p2], rUp[m p2 + om p1]}]];
c3MulTab[{xl_, xh_}, {tl_, th_}] := Module[{st = UnitStep[tl], ost, u, v, mu, mv},
  ost = 1 - st;
  u = st xl + ost xh; v = st xh + ost xl;
  mu = UnitStep[u]; mv = UnitStep[v];
  ndChk[{rDn[u (mu tl + (1 - mu) th)], rUp[v (mv th + (1 - mv) tl)]}]];
c3MulNonneg[{rl_, rh_}, {xl_, xh_}] := Module[{ml = UnitStep[xl], mh = UnitStep[xh]},
  ndChk[{rDn[xl (ml rl + (1 - ml) rh)], rUp[xh (mh rh + (1 - mh) rl)]}]];
(* times an integer of modulus 1 or 2: exact (negation, doubling) unless it overflows *)
c3Scale[iv : {lo_, hi_}, g_Integer] := Which[g == 1, iv, g == -1, {-hi, -lo}, g > 0, ndChk[{g lo, g hi}], True, ndChk[{g hi, g lo}]];
c3SumScaled[list_] := Fold[ndAdd, c3Scale[#[[1]], #[[2]]] & /@ list];

(* upper bounds of the row sums of non-negative upper bounds: pairwise sums, each rounded up *)
c3RowSumUp[x_] := Module[{y = x, k},
  While[(k = Last[Dimensions[y]]) > 1,
    If[OddQ[k], y = Join[y, ConstantArray[0., {Length[y], 1}], 2]; k++];
    y = ndUp[y[[All, 1 ;; k ;; 2]] + y[[All, 2 ;; k ;; 2]]]];
  y[[All, 1]]];

(* the 3 x 3 x 3 spline coefficients of every cell: [cell (ex, ey, ez) -> (ex n + ey) n + ez + 1, components, 3, 3, 3] *)
c3LocalCoefficients[c_, n_] := Flatten[Partition[c, {Length[c], 3, 3, 3}, {1, 1, 1, 1}], 3];

(* indices of the tables for a batch of boxes: lj per direction, and (e, lj) for the function tables *)
c3Gather[tb_, bx_] := Module[{p = 2^bx[[All, 4]], lj},
  lj = Table[p + bx[[All, 4 + k]], {k, 3}];
  {lj, Table[bx[[All, k]] tb["nlj"] + lj[[k]], {k, 3}]}];

(* the points of a box, (q3, q2, q1) with q3 slowest: the index of each point in the three directions *)
c3PointIndex[nq_] := c3PointIndex[nq] = Developer`ToPackedArray /@ {
  Flatten[Table[q1, {q3, nq}, {q2, nq}, {q1, nq}]], Flatten[Table[q2, {q3, nq}, {q2, nq}, {q1, nq}]], Flatten[Table[q3, {q3, nq}, {q2, nq}, {q1, nq}]]};

(* which spline values are needed: for every derivative variant v the components "pv"[v]; for the three stages of the
   sum factorisation the components "p2"[{dx, dy}] and "p1"[dx] *)
c3Plan[varIdx_List, needVP_List] := Module[{pv, keys2, p2, p1},
  pv = Association[Table[v -> Union[Cases[needVP, {v, p_} :> p]], {v, varIdx}]];
  keys2 = Union[$c3Var[[#, 1 ;; 2]] & /@ varIdx];
  p2 = Association[Table[k -> Union @@ Table[pv[v], {v, Select[varIdx, $c3Var[[#, 1 ;; 2]] === k &]}], {k, keys2}]];
  p1 = Association[Table[d -> Union @@ Table[p2[k], {k, Select[keys2, #[[1]] == d &]}], {d, Union[keys2[[All, 1]]]}]];
  <|"pv" -> pv, "p2" -> p2, "p1" -> p1|>];
(* inc(rho^2 Psi) needs, for every non-zero G[s, p, c] with c = (k, m): d_k d_m Psi_p, d_k Psi_p, d_m Psi_p, Psi_p *)
c3DualPlan[] := c3DualPlan[] = c3Plan[Range[10], Union[Flatten[Table[With[{c = $c3Comps[[t[[3]]]], p = t[[2]]},
    {{$c3Second[[t[[3]]]], p}, {$c3First[[c[[1]]]], p}, {$c3First[[c[[2]]]], p}, {1, p}}], {t, $c3Triples}], 1]]];
c3PrimalPlan[] := c3PrimalPlan[] = c3Plan[$c3First, Flatten[Table[{v, a}, {v, $c3First}, {a, 3}], 1]];

(* RIGOROUS spline values at the Gauss points of the boxes, sum-factorised (x, then y, then z) as _contract of the
   reference.  c: exact coefficients [nb, components, 3, 3, 3]; bas: per direction {lo, hi} [nb, 3 (d), 3 (m), nq].
   Returns for every variant v of the plan {lo, hi} [nb, nq^3, components of "pv"[v]], points as c3PointIndex. *)
c3Contract[c_, bas_, plan_, straddle_] := Module[{nq = Last[Dimensions[bas[[1, 1]]]], nb = Length[c], mulT, stack, a1 = <||>, a2 = <||>, out = <||>},
  mulT = If[straddle, ndMul, c3MulTab];
  (* intervals for every point of the new direction (each checked) stacked after the box index *)
  stack[list_] := {Transpose[Developer`ToPackedArray[list[[All, 1]]], {2, 1}], Transpose[Developer`ToPackedArray[list[[All, 2]]], {2, 1}]};
  Do[
    With[{ps = plan["p1"][d1]},
      a1[d1] = stack[Table[
        Fold[ndAdd, Table[c3MulPt[c[[All, ps, m]], {bas[[1, 1]][[All, d1 + 1, m, q]], bas[[1, 2]][[All, d1 + 1, m, q]]}], {m, 3}]],
        {q, nq}]]],                                                       (* [nb, q1, p, i2, i3] *)
    {d1, Keys[plan["p1"]]}];
  Do[
    With[{d1 = k[[1]], d2 = k[[2]], idx = Flatten[Position[plan["p1"][k[[1]]], #] & /@ plan["p2"][k]]},
      a2[k] = stack[Table[
        Fold[ndAdd, Table[mulT[{a1[d1][[1]][[All, All, idx, m]], a1[d1][[2]][[All, All, idx, m]]},
          {bas[[2, 1]][[All, d2 + 1, m, q]], bas[[2, 2]][[All, d2 + 1, m, q]]}], {m, 3}]],
        {q, nq}]]],                                                       (* [nb, q2, q1, p, i3] *)
    {k, Keys[plan["p2"]]}];
  Do[
    With[{k = $c3Var[[v, 1 ;; 2]], d3 = $c3Var[[v, 3]], idx = Flatten[Position[plan["p2"][$c3Var[[v, 1 ;; 2]]], #] & /@ plan["pv"][v]]},
      out[v] = Map[ArrayReshape[#, {nb, nq^3, Length[idx]}] &, stack[Table[
        Fold[ndAdd, Table[mulT[{a2[k][[1]][[All, All, All, idx, m]], a2[k][[2]][[All, All, All, idx, m]]},
          {bas[[3, 1]][[All, d3 + 1, m, q]], bas[[3, 2]][[All, d3 + 1, m, q]]}], {m, 3}]],
        {q, nq}]]]],                                                      (* [nb, q3 q2 q1, p] *)
    {v, Keys[plan["pv"]]}];
  out];

c3WeightsUp[tb_, lj_, iq_] := With[{w = Table[tb["W", "Hi"][[lj[[k]]]][[All, iq[[k]]]], {k, 3}]}, ndUp[ndUp[w[[1]] w[[2]]] w[[3]]]];

(* the stress component s of inc(rho^2 Psi) at the points: Sum over G[s, p, c] != 0, c = (k, m), of
   G (R d_k d_m Psi_p + dR_k d_m Psi_p + dR_m d_k Psi_p + d2_c Psi_p), grouped by the factors R, dR_j, d2_c *)
c3StressComponent[s_, rr_, dr_, d2_, at_] := Module[{ts = Select[$c3Triples, #[[1]] == s &], acc, terms},
  acc = c3MulNonneg[rr, c3SumScaled[Table[{at[$c3Second[[t[[3]]]], t[[2]]], t[[4]]}, {t, ts}]]];
  Do[
    terms = Join @@ Table[With[{k = $c3Comps[[t[[3]], 1]], m = $c3Comps[[t[[3]], 2]]},
        Join[If[k == j, {{at[$c3First[[m]], t[[2]]], t[[4]]}}, {}], If[m == j, {{at[$c3First[[k]], t[[2]]], t[[4]]}}, {}]]], {t, ts}];
    If[terms =!= {}, acc = ndAdd[acc, ndMul[dr[[j]], c3SumScaled[terms]]]],
    {j, 3}];
  Do[acc = ndAdd[acc, ndMul[d2[[t[[3]]]], c3Scale[at[1, t[[2]]], t[[4]]]]], {t, ts}];
  acc];

(* upper bounds, box by box, of the integral of tau:C^-1:tau, tau = tau_p + inc(rho^2 Psi) (dual_compliance) *)
c3DualBatch[sp_, tb_, bx_, loc_, plan_] := Module[
  {n = sp["N"], iq = c3PointIndex[7], nu = sp["nu"], lj, fi, bas, c, a3, at, ef, rho, rr, dr, d2, tau, fx, tr, s1, s2, q},
  {lj, fi} = c3Gather[tb, bx];
  bas = Table[{tb["Bas", "Lo"][[lj[[k]]]], tb["Bas", "Hi"][[lj[[k]]]]}, {k, 3}];
  c = loc[[(bx[[All, 1]] n + bx[[All, 2]]) n + bx[[All, 3]] + 1]];
  a3 = c3Contract[c, bas, plan, tb["Straddle"]];
  at[v_, p_] := With[{j = First[FirstPosition[plan["pv"][v], p]]}, {a3[v][[1]][[All, All, j]], a3[v][[2]][[All, All, j]]}];
  (* the weight rho = eta_x + eta_y + eta_z - a2 (zero on the sphere) and its derivatives *)
  ef = Table[With[{t = tb["Fn", sp["Eta"][[k]]]}, {t["Lo"][[fi[[k]], d + 1]][[All, iq[[k]]]], t["Hi"][[fi[[k]], d + 1]][[All, iq[[k]]]]}], {k, 3}, {d, 0, 2}];
  rho = ndSub[ndAdd[ndAdd[ef[[1, 1]], ef[[2, 1]]], ef[[3, 1]]], ivC[sp["a2"]]];
  rr = ivSqr[rho];
  dr = Table[c3Scale[ndMul[rho, ef[[k, 2]]], 2], {k, 3}];
  d2 = Table[With[{k = $c3Comps[[cc, 1]], m = $c3Comps[[cc, 2]]},
      If[k == m, ndAdd[c3Scale[ndMul[ef[[k, 2]], ef[[m, 2]]], 2], c3Scale[ndMul[rho, ef[[k, 3]]], 2]], c3Scale[ndMul[ef[[k, 2]], ef[[m, 2]]], 2]]], {cc, 6}];
  tau = Table[c3StressComponent[s, rr, dr, d2, at], {s, 6}];
  (* the force-carrying part *)
  Do[
    fx = Table[With[{t = tb["Fn", st[[2 + k]]]}, {t["Lo"][[fi[[k]], 1]][[All, iq[[k]]]], t["Hi"][[fi[[k]], 1]][[All, iq[[k]]]]}], {k, 3}];
    tau[[st[[1]]]] = ndAdd[tau[[st[[1]]]], ivK[ndMul[ndMul[fx[[1]], fx[[2]]], fx[[3]]], st[[2]]]],
    {st, sp["STerms"]}];
  (* tau:C^-1:tau = (1 + nu) (tau_xx^2 + tau_yy^2 + tau_zz^2 + 2 (tau_xy^2 + tau_yz^2 + tau_xz^2)) - nu (trace)^2 *)
  tr = ndAdd[ndAdd[tau[[1]], tau[[2]]], tau[[3]]];
  s1 = ndAdd[ndAdd[ivSqr[tau[[1]]], ivSqr[tau[[2]]]], ivSqr[tau[[3]]]];
  s2 = ndAdd[ndAdd[ivSqr[tau[[4]]], ivSqr[tau[[5]]]], ivSqr[tau[[6]]]];
  q = ndSub[ivK[ndAdd[s1, c3Scale[s2, 2]], 1 + nu], ivK[ivSqr[tr], nu]];
  c3RowSumUp[ndUp[Ramp[q[[2]]] c3WeightsUp[tb, lj, iq]]]];

(* upper bounds, box by box, of the integral of eps(v):C:eps(v), v = (x/2) e_x + spline (primal_energy) *)
c3PrimalBatch[sp_, tb_, bx_, loc_, plan_] := Module[{n = sp["N"], iq = c3PointIndex[3], lj, bas, c, a3, dv, e, tr, q},
  lj = First[c3Gather[tb, bx]];
  bas = Table[{tb["Bas", "Lo"][[lj[[k]]]], tb["Bas", "Hi"][[lj[[k]]]]}, {k, 3}];
  c = loc[[(bx[[All, 1]] n + bx[[All, 2]]) n + bx[[All, 3]] + 1]];
  a3 = c3Contract[c, bas, plan, tb["Straddle"]];
  dv[a_, k_] := {a3[$c3First[[k]]][[1]][[All, All, a]], a3[$c3First[[k]]][[2]][[All, All, a]]};      (* d_k v_a *)
  e = {ndAdd[dv[1, 1], ivC[1/2]], dv[2, 2], dv[3, 3], ndAdd[dv[1, 2], dv[2, 1]], ndAdd[dv[2, 3], dv[3, 2]], ndAdd[dv[1, 3], dv[3, 1]]};
  tr = ndAdd[ndAdd[e[[1]], e[[2]]], e[[3]]];
  q = ndAdd[ndAdd[ivK[ivSqr[tr], sp["lam"]], ivK[ndAdd[ndAdd[ivSqr[e[[1]]], ivSqr[e[[2]]]], ivSqr[e[[3]]]], 2 sp["mu"]]],
    ivK[ndAdd[ndAdd[ivSqr[e[[4]]], ivSqr[e[[5]]]], ivSqr[e[[6]]]], sp["mu"]]];
  c3RowSumUp[ndUp[Ramp[q[[2]]] c3WeightsUp[tb, lj, iq]]]];

(* RIGOROUS upper bounds (exact rationals) of the integrals over the boxes of one octant; the field arrays must be
   packed arrays of machine numbers (they are taken as the exact binary numbers they are) *)
c3DualUpper[sp_, boxes_, psi_, monitor_ : None] := c3Integrate[sp, boxes, psi, 7, c3DualBatch, c3DualPlan[], 128, monitor];
c3PrimalUpper[sp_, boxes_, v_, monitor_ : None] := c3Integrate[sp, boxes, v, 3, c3PrimalBatch, c3PrimalPlan[], 2048, monitor];
c3Integrate[sp_, boxes_, field_, nq_, batchFn_, plan_, batch_, monitor_] := ndQuiet@Module[
  {n = sp["N"], tb, loc, nbx = Length[boxes], parts = {}, t0 = AbsoluteTime[]},
  tb = c3Tables[n, Max[boxes[[All, 4]]], nq, sp["kflat"]];
  loc = c3LocalCoefficients[field, n];
  Do[
    AppendTo[parts, batchFn[sp, tb, boxes[[s ;; Min[s + batch - 1, nbx]]], loc, plan]];
    If[monitor =!= None, monitor[Min[s + batch - 1, nbx], nbx, AbsoluteTime[] - t0]],
    {s, 1, nbx, batch}];
  exactTotal[Join @@ parts]];

(* ------------------------------------------------------------------------------------------------ *)
(* the solver (heuristic: machine arithmetic; any result gives valid bounds)                          *)
(* ------------------------------------------------------------------------------------------------ *)

(* tensor-product spline values at the points: X, Y, Z [nb, 3 (m), nq] -> [nb, points, 27 (m1 slowest)] *)
c3TensorFloat[x_, y_, z_, iq_] := Module[{xp = x[[All, All, iq[[1]]]], yp = y[[All, All, iq[[2]]]], zp = z[[All, All, iq[[3]]]]},
  Transpose[Developer`ToPackedArray[Flatten[Table[xp[[All, m1]] yp[[All, m2]] zp[[All, m3]], {m1, 3}, {m2, 3}, {m3, 3}], 2]], {3, 1, 2}]];

c3WeightsMid[tb_, lj_, iq_] := Times @@ Table[tb["W", "Mid"][[lj[[k]]]][[All, iq[[k]]]], {k, 3}];

(* Gram matrices of the weighted stress of the local coefficients, box by box: [nb, 162, 162] and [nb, 162] *)
c3DualGramFloat[sp_, tb_, bx_, ut_] := Module[{n = sp["N"], iq = c3PointIndex[7], nb = Length[bx], np = 343, lj, fi, bm, bd, ef, rho, rr, dr, h, t, sw, tp, taup, tpp},
  {lj, fi} = c3Gather[tb, bx];
  bm = Table[tb["Bas", "Mid"][[lj[[k]]]], {k, 3}];
  bd = Table[c3TensorFloat[bm[[1]][[All, v[[1]] + 1]], bm[[2]][[All, v[[2]] + 1]], bm[[3]][[All, v[[3]] + 1]], iq], {v, $c3Var}];
  ef = Table[tb["Fn", sp["Eta"][[k]], "Mid"][[fi[[k]], d + 1]][[All, iq[[k]]]], {k, 3}, {d, 0, 2}];
  rho = ef[[1, 1]] + ef[[2, 1]] + ef[[3, 1]] - N[sp["a2"]];
  rr = rho^2; dr = Table[2 rho ef[[k, 2]], {k, 3}];
  h = Table[With[{k = $c3Comps[[cc, 1]], m = $c3Comps[[cc, 2]]},
      rr bd[[$c3Second[[cc]]]] + dr[[k]] bd[[$c3First[[m]]]] + dr[[m]] bd[[$c3First[[k]]]] +
        (2 ef[[k, 2]] ef[[m, 2]] + If[k == m, 2 rho ef[[k, 3]], 0. rho]) bd[[1]]], {cc, 6}];
  t = ConstantArray[0., {6, 6, nb, np, 27}];
  Do[t[[tr[[1]], tr[[2]]]] += tr[[4]] h[[tr[[3]]]], {tr, $c3Triples}];
  sw = Sqrt[c3WeightsMid[tb, lj, iq]];
  tp = ArrayReshape[sw Transpose[ut . t, {3, 4, 1, 2, 5}], {nb, 6 np, 162}];
  taup = Total[Table[N[st[[2]]] Times @@ Table[tb["Fn", st[[2 + k]], "Mid"][[fi[[k]], 1]][[All, iq[[k]]]], {k, 3}], {st, sp["STerms"]}]];
  tpp = ArrayReshape[Outer[Times, sw taup, ut[[All, 1]]], {nb, 6 np}];
  {Transpose[#] . # & /@ tp, MapThread[Transpose[#1] . #2 &, {tp, tpp}]}];

c3PrimalGramFloat[sp_, tb_, bx_, ut_] := Module[{iq = c3PointIndex[3], nb = Length[bx], np = 27, lj, bm, bd, bmat, sw, bp, e0p},
  lj = First[c3Gather[tb, bx]];
  bm = Table[tb["Bas", "Mid"][[lj[[k]]]], {k, 3}];
  bd = Table[c3TensorFloat[bm[[1]][[All, v[[1]] + 1]], bm[[2]][[All, v[[2]] + 1]], bm[[3]][[All, v[[3]] + 1]], iq], {v, $c3Var[[$c3First]]}];
  bmat = Developer`ToPackedArray[N[Table[Sum[$c3GS[[s, a, k]] bd[[k]], {k, 3}], {s, 6}, {a, 3}]]];
  sw = Sqrt[c3WeightsMid[tb, lj, iq]];
  bp = ArrayReshape[sw Transpose[ut . bmat, {3, 4, 1, 2, 5}], {nb, 6 np, 81}];
  e0p = ArrayReshape[Outer[Times, sw, ut . {0.5, 0., 0., 0., 0., 0.}], {nb, 6 np}];
  {Transpose[#] . # & /@ bp, MapThread[Transpose[#1] . #2 &, {bp, e0p}]}];

c3LocalToGlobal[cells_, nc_, m_] := Module[{loc = Flatten[Table[(i1 m + i2) m + i3, {i1, 0, 2}, {i2, 0, 2}, {i3, 0, 2}]], pl},
  pl = Flatten[Table[(p - 1) m^3 + loc + 1, {p, nc}]];
  Outer[Plus, (cells[[All, 1]] m + cells[[All, 2]]) m + cells[[All, 3]], pl]];

(* the cell matrices of a batch (summed over the boxes of each cell) in the reduced unknowns: {u, A_u, b_u} with
   distinct reduced indices u (full coefficients tied to the same reduced one are combined with their signs) *)
c3CellBlocks[bx_, gram_, lin_, nc_, m_, red_, sg_] := Module[{groups, cells, g},
  groups = SplitBy[Range[Length[bx]], bx[[#, 1 ;; 3]] &];
  cells = bx[[First /@ groups, 1 ;; 3]];
  g = c3LocalToGlobal[cells, nc, m];
  Table[
    Module[{grp = groups[[i]], r = red[[g[[i]]]], s = sg[[g[[i]]]], keep, u, sel},
      keep = Flatten[Position[Unitize[s], 1]];
      If[keep === {}, Nothing,
        u = DeleteDuplicates[r[[keep]]];
        sel = SparseArray[Transpose[{Range[Length[keep]], Lookup[AssociationThread[u -> Range[Length[u]]], r[[keep]]]}] -> N[s[[keep]]],
          {Length[keep], Length[u]}];
        {u, Normal[Transpose[sel] . Total[gram[[grp]]][[keep, keep]] . sel], Normal[Transpose[sel] . Total[lin[[grp]]][[keep]]]}]],
    {i, Length[groups]}]];

c3SparseSum[pos_, val_, dims_] := Module[{old = "TreatRepeatedEntries" /. ("SparseArrayOptions" /. SystemOptions["SparseArrayOptions"])},
  Internal`WithLocalSettings[SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 1}],
    SparseArray[pos -> val, dims],
    SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> old}]]];

(* minimise c.A.c/2 - rhs.c as solve_reduced of the reference: rows with a negligible diagonal are pinned, the
   diagonal gets 10^-13 of its maximum; dense LU up to 7000 unknowns, Jacobi-preconditioned CG beyond *)
c3SolveReduced[a_, rhs_] := Module[{n = Length[rhs], d, scale, dead, b = rhs, ad, y},
  d = Normal[Diagonal[a]];
  scale = Max[Append[d, 0.]]; If[! (scale > 0), scale = 1.];
  dead = Flatten[Position[UnitStep[1.*^-14 scale - d], 1]];
  b[[dead]] = 0.;
  y = If[n <= 7000,
    ad = Normal[a];
    If[dead =!= {}, ad[[dead, All]] = 0.; ad[[All, dead]] = 0.; Do[ad[[i, i]] = scale, {i, dead}]];
    Do[ad[[i, i]] += 1.*^-13 scale, {i, n}];
    Quiet[LinearSolve[ad, b]],
    c3PCG[a, b, dead, d, scale]];
  If[VectorQ[y, MachineNumberQ] && Length[y] == n, y, ConstantArray[0., n]]];

c3PCG[a_, b_, dead_, d_, scale_] := Module[{n = Length[b], keep, kd, ak, dd, x, r, z, p, rz, rzn, ap, al, it = 0, nb = Norm[b]},
  keep = ConstantArray[1., n]; keep[[dead]] = 0.;
  kd = SparseArray[Band[{1, 1}] -> keep];
  ak = kd . a . kd + SparseArray[Band[{1, 1}] -> ((1. - keep) + 1.*^-13 scale)];
  dd = d; dd[[dead]] = 1.;
  x = 0. b; r = b; z = r/dd; p = z; rz = r . z;
  While[it < 3000 && Norm[r] > 1.*^-10 nb && rz != 0.,
    ap = ak . p; al = rz/(p . ap); x += al p; r -= al ap; z = r/dd; rzn = r . z; p = z + (rzn/rz) p; rz = rzn; it++];
  x];

(* the field that minimises the discretised energy (displacement) or complementary energy (stress potential) on the
   boxes of the given depth, as a full coefficient array satisfying the ties exactly *)
c3Solve[sp_, boxes_, ties_, nq_, nc_, gramFn_, ut_, batch_] := Module[{n = sp["N"], m = sp["M"], tb, red, sg, nred, dense, amat, bvec, blocks, gl, y},
  tb = c3Tables[n, Max[boxes[[All, 4]]], nq, sp["kflat"]];
  {red, sg} = c3RedMap[ties[[1]]];
  nred = Dimensions[ties[[1]]][[2]];
  dense = nred <= 7000;
  amat = If[dense, ConstantArray[0., {nred, nred}], SparseArray[{}, {nred, nred}, 0.]];
  bvec = ConstantArray[0., nred];
  Do[
    With[{bx = boxes[[s ;; Min[s + batch - 1, Length[boxes]]]]},
      gl = gramFn[sp, tb, bx, ut];
      blocks = c3CellBlocks[bx, gl[[1]], gl[[2]], nc, m, red, sg];
      Do[bvec[[blk[[1]]]] += blk[[3]], {blk, blocks}];
      If[dense,
        Do[amat[[blk[[1]], blk[[1]]]] += blk[[2]], {blk, blocks}],
        If[blocks =!= {}, amat += c3SparseSum[Join @@ (Tuples[{#[[1]], #[[1]]}] & /@ blocks), Join @@ (Flatten[#[[2]]] & /@ blocks), {nred, nred}]]]],
    {s, 1, Length[boxes], batch}];
  y = c3SolveReduced[amat, -bvec];
  Developer`ToPackedArray[N[ArrayReshape[Normal[ties[[1]] . y], {nc, m, m, m}]]]];
c3SolveDual[sp_, boxes_, ties_] := c3Solve[sp, boxes, ties, 7, 6, c3DualGramFloat, CholeskyDecomposition[N[c3Compliance[sp["nu"]]]], 16];
c3SolvePrimal[sp_, boxes_, ties_] := c3Solve[sp, boxes, ties, 3, 3, c3PrimalGramFloat, CholeskyDecomposition[N[c3Stiffness[sp["lam"], sp["mu"]]]], 512];

(* ------------------------------------------------------------------------------------------------ *)
(* fields written by the Python reference (NumPy .npz archives)                                       *)
(* ------------------------------------------------------------------------------------------------ *)

(* one array of a .npy file (little-endian float64 or int64 in C order); a 0-dimensional array gives a number.  The
   header is parsed by string matching only. *)
c3ReadNPY[file_String] := Module[{s = OpenRead[file, BinaryFormat -> True], magic, ver, hlen, header, type, sh, shape, data},
  magic = BinaryReadList[s, "Byte", 6]; ver = BinaryReadList[s, "Byte", 2];
  If[magic =!= {147, 78, 85, 77, 80, 89} || ! MatchQ[ver, {1 | 2, _}], Close[s]; Throw[failure["MissingData", "Not a NumPy array file: " <> file], "CertiField"]];
  hlen = BinaryRead[s, If[ver[[1]] == 1, "UnsignedInteger16", "UnsignedInteger32"], ByteOrdering -> -1];
  header = FromCharacterCode[BinaryReadList[s, "Byte", hlen]];
  type = Which[StringContainsQ[header, "'descr': '<f8'"], "Real64", StringContainsQ[header, "'descr': '<i8'"], "Integer64", True, None];
  data = If[type === None, {}, BinaryReadList[s, type, ByteOrdering -> -1]]; Close[s];
  sh = StringCases[header, "'shape': (" ~~ x : Shortest[___] ~~ ")" :> x];
  shape = If[sh === {}, $Failed, FromDigits /@ StringCases[First[sh], DigitCharacter ..]];
  If[type === None || shape === $Failed || ! StringContainsQ[header, "'fortran_order': False"] || Length[data] =!= Times @@ shape,
    Throw[failure["MissingData", "Unsupported NumPy array (need little-endian float64 or int64 in C order): " <> file], "CertiField"]];
  If[shape === {}, First[data], ArrayReshape[data, shape]]];

(* the arrays of a .npz archive written by the reference (psi, v, K_low, K_up), by name *)
c3ReadNPZ[file_String] := Module[{tmp, files, z},
  If[! FileExistsQ[file], Throw[failure["MissingData", "No file " <> file], "CertiField"]];
  tmp = CreateDirectory[];
  z = Catch[
    files = Quiet@ExtractArchive[file, tmp];
    If[! ListQ[files], Throw[failure["MissingData", "Not a NumPy archive: " <> file], "CertiField"]];
    Association[(FileBaseName[#] -> c3ReadNPY[#]) & /@ Select[files, StringEndsQ[#, ".npy"] &]], "CertiField"];
  Quiet@DeleteDirectory[tmp, DeleteContents -> True];
  If[! AssociationQ[z], Throw[z, "CertiField"]];
  z];

(* a field supplied by the user: Automatic (the solver's), an array, or a .npz file with the array key; it must be a
   real array of the given dimensions that satisfies the ties exactly *)
c3FieldInput[x_, key_String, dims_, ties_, what_String] := Module[{arr = x, z},
  If[x === Automatic || x === None, Return[None, Module]];
  If[StringQ[x],
    z = c3ReadNPZ[x];
    If[! KeyExistsQ[z, key], Throw[failure["MissingData", "The file " <> x <> " has no array " <> key <> "."], "CertiField"]];
    arr = z[key]];
  arr = ndRealArray[arr, dims];
  If[arr === $Failed, Throw[failure["Bad" <> what, "The " <> ToLowerCase[what] <> " must be a real array of dimensions " <> ToString[dims] <> " (spline coefficients)."], "CertiField"]];
  If[! c3TiesHold[ties, arr],
    Throw[failure["Bad" <> what, "The " <> ToLowerCase[what] <> " violates the mirror, free-face or grip conditions of the spline coefficients."], "CertiField"]];
  arr];

(* ------------------------------------------------------------------------------------------------ *)
(* the model                                                                                         *)
(* ------------------------------------------------------------------------------------------------ *)

certifyElasticCavity3D[params_Association, cand_, opts_List] := Module[
  {sp = c3Spec[params], t0 = AbsoluteTime[], verbose, say, mon, n, m, dp, tv, tp, v, psi, bc, bcp, kup, qup, klow, margin},
  verbose = TrueQ[OptionValue[CertifiedVariationalBound, opts, "Verbose"]];
  say[x__] := If[verbose, Print["  ", x, "   (", ndSec[AbsoluteTime[] - t0], " s)"]];
  mon[what_] := If[verbose, Module[{last = 0}, Function[{k, tot, sec}, If[k == tot || k >= last + tot/5, last = k; Print["    ", what, " ", k, "/", tot, " boxes   (", ndSec[sec], " s)"]]]], None];
  n = sp["N"]; m = sp["M"]; dp = sp["Depths"];
  tv = c3Ties[n, "V", sp["Clamped"]]; tp = c3Ties[n, "Psi", sp["Clamped"]];
  v = c3FieldInput[cand, "v", {3, m, m, m}, tv, "Candidate"];
  psi = c3FieldInput[Lookup[params, "Field", Automatic], "psi", {6, m, m, m}, tp, "Field"];
  If[verbose, Print["  ElasticCavity3D: ", If[sp["Cavity"], "cavity", "solid block"], ", ", ToLowerCase[sp["Grips"]], " grips, ", n, "^3 cells per octant, depths ", dp]];
  (* upper bound *)
  If[v === None, v = c3SolvePrimal[sp, c3Boxes[n, dp[[3]], sp["a2"]], tv]; say["displacement found on the boxes of depth ", dp[[3]]]];
  bcp = c3Boxes[n, dp[[4]], sp["a2"]];
  kup = 8 c3PrimalUpper[sp, bcp, v, mon["displacement"]];
  say["K <= ", ndNum[kup], " (", Length[bcp], " boxes of depth ", dp[[4]], ")"];
  (* lower bound *)
  If[psi === None, psi = c3SolveDual[sp, c3Boxes[n, dp[[1]], sp["a2"]], tp]; say["stress potential found on the boxes of depth ", dp[[1]]]];
  bc = c3Boxes[n, dp[[2]], sp["a2"]];
  qup = 8 c3DualUpper[sp, bc, psi, mon["stress"]];
  If[! (qup > 0), Throw[failure["Inconclusive", "The complementary energy could not be bounded."], "CertiField"]];
  klow = sp["Area"]^2/qup;
  say["K >= ", ndNum[klow], " (", Length[bc], " boxes of depth ", dp[[2]], ")"];
  (* the claims are the computed bounds weakened by one part in 10^9 of the larger modulus *)
  margin = Max[Abs[klow], Abs[kup]]/10^9;
  <|"LowerBound" -> klow - margin, "UpperBound" -> kup + margin, "Seconds" -> AbsoluteTime[] - t0,
    "Psi" -> psi, "V" -> v, "Boxes" -> {Length[bc], Length[bcp]}, "Spec" -> sp|>];

elasticCavityDescription[params_, res_] := Module[{sp = res["Spec"], d = Max[displayDigits[res], 7], n, dp, body, free, grips, adm},
  n = sp["N"]; dp = sp["Depths"];
  body = If[sp["Cavity"], "the block (-1,1)^3 minus the closed ball |x| \[LessEqual] 1/2 at its centre", "the block (-1,1)^3 (no cavity)"];
  free = If[sp["Cavity"], "the faces y = \[PlusMinus]1, z = \[PlusMinus]1 and the sphere |x| = 1/2 are free of traction", "the faces y = \[PlusMinus]1, z = \[PlusMinus]1 are free of traction"];
  grips = If[sp["Clamped"], "the faces x = \[PlusMinus]1 are clamped to rigid grips moved apart by \[Delta], u = (\[PlusMinus]\[Delta]/2, 0, 0)",
    "the faces x = \[PlusMinus]1 rest on frictionless rigid grips moved apart by \[Delta], u_x = \[PlusMinus]\[Delta]/2 with the tangential displacement free"];
  adm = If[sp["Clamped"], "u = (\[PlusMinus]\[Delta]/2, 0, 0) on x = \[PlusMinus]1", "u_x = \[PlusMinus]\[Delta]/2 on x = \[PlusMinus]1"];
  <|"Units" -> "stiffness K: force on a grip per unit separation of the grips (Young's modulus 1)", "Quantity" -> "K",
    "Theorem" -> Row[{"Let \[CapitalOmega] be ", body, ", isotropic linear elastic with Young's modulus 1 and Poisson ratio 3/10; ", free, "; ", grips,
      ". The stiffness K, defined by min over admissible u of (1/2) \[Integral]_\[CapitalOmega] \[CurlyEpsilon](u):C:\[CurlyEpsilon](u) dx = K \[Delta]^2/2, satisfies ",
      lowerString[res["LowerBound"], d], " \[LessEqual] K \[LessEqual] ", upperString[res["UpperBound"], d], "."}],
    "AdmissibleClass" -> Row[{"all u in H^1(\[CapitalOmega]; \[DoubleStruckCapitalR]^3) with ", adm, "; no condition on the other surfaces"}],
    "LowerBoundMethod" -> Row[{"complementary energy, K \[GreaterEqual] 16/\[Integral] \[Tau]:C^-1:\[Tau], for the equilibrated stress \[Tau] = \[Tau]_p + inc(\[Rho]^2 \[CapitalPsi]) ",
      "(Beltrami potential of C^1 quadratic B-splines on ", n, "^3 cells per octant, \[Rho] = 0 on the sphere: zero traction there exactly), integrated over ",
      res["Boxes"][[1]], " octree boxes (depth ", dp[[2]], ") that cover the body; exact Gauss rules, outward-rounded arithmetic"}],
    "UpperBoundMethod" -> Row[{"energy of a C^1 quadratic spline displacement on ", n, "^3 cells per octant, integrated over ", res["Boxes"][[2]],
      " octree boxes (depth ", dp[[4]], ") that cover the body; exact Gauss rules, outward-rounded arithmetic"}],
    "TrustedBase" -> "IEEE-754 machine arithmetic (+, \[Minus], \[Times]) with outward rounding, Wolfram Language exact rational arithmetic, and the argument of Reference/Curved3D/NOTES.md (\"Why the bounds are true\"; checked by its author, the stored fields of the reference also by an exact rational re-check)",
    "NotClaimed" -> "that either spline field is optimal; the gap bounds the distance of both from the true stiffness",
    "References" -> "Reference/Curved3D/NOTES.md; complementary energy (Prager\[Dash]Synge 1947); Beltrami stress functions (Gurtin 1963)."|>];

(* re-verification: the ties of the stored fields are checked exactly and both bounds recomputed from them *)
verifyElasticCavity3D[a_Association] := Module[{c = a["Certificate"], sp, n, m, psi, v, kup, klow},
  sp = c3Spec[a["Parameters"]]; n = sp["N"]; m = sp["M"];
  psi = ndRealArray[Lookup[c, "Psi", $Failed], {6, m, m, m}];
  v = ndRealArray[Lookup[c, "V", $Failed], {3, m, m, m}];
  If[psi === $Failed || v === $Failed, Return[False, Module]];
  If[! (c3TiesHold[c3Ties[n, "Psi", sp["Clamped"]], psi] && c3TiesHold[c3Ties[n, "V", sp["Clamped"]], v]), Return[False, Module]];
  kup = 8 c3PrimalUpper[sp, c3Boxes[n, sp["Depths"][[4]], sp["a2"]], v];
  klow = sp["Area"]^2/(8 c3DualUpper[sp, c3Boxes[n, sp["Depths"][[2]], sp["a2"]], psi]);
  klow >= exactNumber[a["LowerBound"]] && kup <= exactNumber[a["UpperBound"]]];

CertiFieldModelInfo["ElasticCavity3D"] := <|
  "Name" -> "ElasticCavity3D",
  "Energy" -> "E[u] = (1/2) \[Integral] \[CurlyEpsilon](u):C:\[CurlyEpsilon](u) dx over the block (-1,1)^3 minus the ball |x| \[LessEqual] 1/2 (isotropic, Young's modulus 1, Poisson ratio 3/10); the faces x = \[PlusMinus]1 are rigid grips moved apart by \[Delta], all other surfaces are free. Reported: the stiffness K, min E = K \[Delta]^2/2",
  "Parameters" -> <|"Grips" -> "\"Clamped\" (u = (\[PlusMinus]\[Delta]/2, 0, 0) on the grips, default) or \"Frictionless\" (u_x = \[PlusMinus]\[Delta]/2 only)",
    "Cells" -> "spline cells per half-edge, a multiple of 4 (default 8)",
    "Depths" -> "{stress solve, stress certificate, displacement solve, displacement certificate}: depths of the octrees of boxes in each cell (default {2, 3, 3, 5})",
    "Cavity" -> "True (default), or False for the solid block, whose stiffness between frictionless grips is 2",
    "Field" -> "optional: the coefficients of the stress potential, an array 6 x (Cells+2)^3, or a .npz file written by Reference/Curved3D (array psi)"|>,
  "Candidate" -> "optional: the coefficients of the displacement spline, an array 3 x (Cells+2)^3, or a .npz file written by Reference/Curved3D (array v)",
  "AdmissibleClass" -> "all u in H^1 of the body with the grip condition on x = \[PlusMinus]1",
  "LowerBound" -> "complementary energy K \[GreaterEqual] 16/\[Integral] \[Tau]:C^-1:\[Tau] with a Beltrami stress \[Tau] = \[Tau]_p + inc(\[Rho]^2 \[CapitalPsi]) that is traction-free on the sphere exactly; integral over boxes that cover the body, outward-rounded arithmetic",
  "UpperBound" -> "energy of a spline displacement over boxes that cover the body, outward-rounded arithmetic",
  "Reference" -> "Reference/Curved3D/NOTES.md",
  "Status" -> "Supported"|>;

End[];
