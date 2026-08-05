#!/usr/bin/env wolframscript

ClearAll["Global`*"];

(* ====================================================== *)
(* Algebraic Rules & Simplification                       *)
(* ====================================================== *)

Clear[\[Omega]];

(* Replacement rules for omega reductions *)
omegaRules = {
  Conjugate[\[Omega]] -> \[Omega]^2,
  Conjugate[\[Omega]^2] -> \[Omega],
  \[Omega]^n_Integer /; n > 2 :> \[Omega]^Mod[n, 3],

  (* Exact linear binomial sum simplifications *)
  1 + \[Omega] -> -\[Omega]^2,
  \[Omega] + 1 -> -\[Omega]^2,
  1 + \[Omega]^2 -> -\[Omega],
  \[Omega]^2 + 1 -> -\[Omega],
  \[Omega] + \[Omega]^2 -> -1,
  \[Omega]^2 + \[Omega] -> -1,

  (* Full relation *)
  1 + \[Omega] + \[Omega]^2 -> 0
};

(* Full simplification forcing reduction modulo (1 + omega + omega^2) *)
SimplifyMatrix[M_] := Module[{reduced},
  reduced = M /. omegaRules;
  reduced = Map[
    Function[elem,
      Last[PolynomialReduce[elem, {1 + \[Omega] + \[Omega]^2}, {\[Omega]}]]
    ],
    reduced,
    {-1}
  ];
  FullSimplify[reduced /. omegaRules] /. omegaRules
];

(* Factor out 1/3 from matrices for clean visual presentation *)
FactorThird[M_] := If[M === ConstantArray[0, {9, 9}], M, Inactive[Times][1/3, MatrixForm[SimplifyMatrix[3 * M]]]];

(* Matrix generators P (Shift) and Q (Clock) *)
P = {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}};
Q = {{1, 0, 0}, {0, \[Omega], 0}, {0, 0, \[Omega]^2}};

(* ====================================================== *)
(* Tensor Product Block Construction                      *)
(* ====================================================== *)

(* Computes G x G* + G* x G *)
TensorBlock[G_] := Module[{GStar},
  GStar = ConjugateTranspose[G];
  SimplifyMatrix[KroneckerProduct[G, GStar] + KroneckerProduct[GStar, G]]
];

(* The 4 elementary tensor product blocks *)
P1 = TensorBlock[P];               (* P x P* + P* x P *)
P2 = TensorBlock[Q];               (* Q x Q* + Q* x Q *)
P3 = TensorBlock[P . Q];           (* PQ x (PQ)* + (PQ)* x PQ *)
P4 = TensorBlock[(P . P) . Q];    (* P^2 Q x (P^2 Q)* + (P^2 Q)* x P^2 Q *)

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
