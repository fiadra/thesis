#!/usr/bin/env wolframscript

ClearAll["Global`*"];

(* ====================================================== *)
(* Algebraic Rules & Simplification                       *)
(* ====================================================== *)

Clear[\[Omega]];

(* Algebraic replacement rules for omega *)
omegaRules = {
  Conjugate[\[Omega]] -> \[Omega]^2,
  Conjugate[\[Omega]^2] -> \[Omega],

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

(* Force omega^3 -> 1 and polynomial reduction *)
SimplifyOmegaExpr[expr_] := Module[{poly},
  (* Step 1: Replace omega^3 with 1 at the power structure level *)
  poly = expr /. Power[\[Omega], n_Integer] :> \[Omega]^Mod[n, 3];

  (* Step 2: Algebraically reduce modulo (1 + omega + omega^2) *)
  poly = Last[PolynomialReduce[poly, {1 + \[Omega] + \[Omega]^2}, {\[Omega]}]];

  (* Step 3: Final application of linear rules *)
  FullSimplify[poly /. omegaRules] /. omegaRules
];

SimplifyMatrix[M_] := Map[SimplifyOmegaExpr, M, {-1}];

(* ====================================================== *)
(* Matrix Definitions A1 through A8                       *)
(* ====================================================== *)

A8 = {
  {1, 0, 0},
  {0, \[Omega]^2, 0},
  {0, 0, \[Omega]}
};

A7 = {
  {1, 0, 0},
  {0, \[Omega], 0},
  {0, 0, \[Omega]^2}
};

A6 = {
  {0, 1, 0},
  {0, 0, \[Omega]},
  {\[Omega]^2, 0, 0}
};

A5 = {
  {0, 0, 1},
  {\[Omega]^2, 0, 0},
  {0, \[Omega], 0}
};

A4 = {
  {0, 1, 0},
  {0, 0, \[Omega]^2},
  {\[Omega], 0, 0}
};

A3 = {
  {0, 0, 1},
  {\[Omega], 0, 0},
  {0, \[Omega]^2, 0}
};

A2 = {
  {0, 1, 0},
  {0, 0, 1},
  {1, 0, 0}
};

A1 = {
  {0, 0, 1},
  {1, 0, 0},
  {0, 1, 0}
};

A[8] = A8;
A[7] = A7;
A[6] = A6;
A[5] = A5;
A[4] = A4;
A[3] = A3;
A[2] = A2;
A[1] = A1;

(* Complex conjugate helper for A_k* *)
AStar[M_] := Conjugate[M];

(* ====================================================== *)
(* Calculation of Sums sum_{k=m}^8 A_k (x) A_k*           *)
(* ====================================================== *)

ComputeTensorSum[m_Integer] := Module[{sum},
  sum = Sum[
    KroneckerProduct[A[k], AStar[A[k]]],
    {k, m, 8}
  ];
  SimplifyMatrix[sum]
];

(* ====================================================== *)
(* Execution and PDF Export                               *)
(* ====================================================== *)

mValues = {8, 7, 6, 5, 4, 3, 2, 1};

pdfPages = Table[
  Module[{sumMatrix, labelStr},
    sumMatrix = ComputeTensorSum[m];
    labelStr = "Sum_{k=" <> ToString[m] <> "}^{8}  A_k \[TensorProduct] A_k*";

    Column[{
      Text[Style[labelStr, Bold, 16]],
      Spacer[12],
      MatrixForm[sumMatrix]
    }, Alignment -> Center]
  ],
  {m, mValues}
];

(* Export to PDF silently without standard output *)
Export["MariuszMatricesCheck.pdf", pdfPages];
Print["Results for m=1..8 exported successfully to MariuszMatricesCheck.pdf"];
