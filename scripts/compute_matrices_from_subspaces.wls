#!/usr/bin/env wolframscript

ClearAll["Global`*"];

(* ====================================================== *)
(* Algebraic Rules & Simplification                       *)
(* ====================================================== *)

Clear[\[Omega]];

(* Expand conjugates and reduce all omega powers down to <= 2 *)
omegaRules = {
  Conjugate[\[Omega]] -> \[Omega]^2,
  Conjugate[\[Omega]^2] -> \[Omega],
  \[Omega]^n_Integer /; n > 2 :> \[Omega]^Mod[n, 3],
  \[Omega]^2 + \[Omega] + 1 -> 0
};

SimplifyMatrix[M_] := FullSimplify[M /. omegaRules] /. omegaRules;

(* Factor out 1/3 from matrices for clean visual presentation *)
FactorThird[M_] := If[M === ConstantArray[0, {9, 9}], M, Inactive[Times][1/3, MatrixForm[SimplifyMatrix[3 * M]]]];

(* Vectorization *)
vec[X_] := Flatten[Normal[X]];

(* Matrix generators P (Updated Shift) and Q (Clock) *)
P = {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}};
Q = {{1, 0, 0}, {0, \[Omega], 0}, {0, 0, \[Omega]^2}};

(* ====================================================== *)
(* Block Projector Construction                           *)
(* ====================================================== *)

Projector[G1_, G2_] := Module[{v1, v2},
  v1 = vec[G1];
  v2 = vec[G2];
  SimplifyMatrix[(KroneckerProduct[v1, Conjugate[v1]] + KroneckerProduct[v2, Conjugate[v2]]) / 3]
];

(* The 4 elementary 2D Fourier-invariant blocks *)
P1 = Projector[{{1,0,0},{0,1,0},{0,0,1}}, {{0,1,0},{0,0,1},{1,0,0}} + {{0,0,1},{1,0,0},{0,1,0}}];
P2 = Projector[P, P . P];
P3 = Projector[P . Q, (P . P) . (Q . Q)];
P4 = Projector[(P . P) . Q, P . (Q . Q)];

blocks = {P1, P2, P3, P4};
Psl3 = SimplifyMatrix[Total[blocks]];

(* ====================================================== *)
(* Subspace Aggregation                                   *)
(* ====================================================== *)

(* Collect all matrices in order of m = 0, 2, 4, 6, 8 *)
allSubspaces = Join[
  (* m = 0 *)
  {{"m = 0 : V_0", ConstantArray[0, {9, 9}]}},

  (* m = 2 *)
  {
    {"m = 2 : P_1", P1},
    {"m = 2 : P_2", P2},
    {"m = 2 : P_3", P3},
    {"m = 2 : P_4", P4}
  },

  (* m = 4 *)
  Table[
    {"m = 4 : W_" <> ToString[pair[[1]]] <> ToString[pair[[2]]],
     SimplifyMatrix[blocks[[pair[[1]]]] + blocks[[pair[[2]]]]]},
    {pair, Subsets[Range[4], {2}]}
  ],

  (* m = 6 *)
  Table[
    {"m = 6 : U_" <> ToString[k],
     SimplifyMatrix[Psl3 - blocks[[k]]]},
    {k, 1, 4}
  ],

  (* m = 8 *)
  {{"m = 8 : sl_3", Psl3}}
];

(* ====================================================== *)
(* Layout and PDF Export                                 *)
(* ====================================================== *)

(* Display each matrix on one single horizontal row with 1/3 factored out *)
pdfPages = Table[
  Row[{
    Text[Style[item[[1]] <> "  =  ", Bold, 12]],
    If[MatrixQ[item[[2]]], FactorThird[item[[2]]], MatrixForm[item[[2]]]]
  }, Alignment -> Center],
  {item, allSubspaces}
];

(* Export to single multi-page PDF silently *)
Export["SelfReflectiveSubspaces.pdf", pdfPages];
