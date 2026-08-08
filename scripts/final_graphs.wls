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

    (* Reduce all integer powers modulo 3 *)
    \[Omega]^n_Integer :> \[Omega]^Mod[n, 3],

    (* Convert omega^0 to 1 *)
    \[Omega]^0 -> 1,

    (* Linear identities *)
    1 + \[Omega] -> -\[Omega]^2,
    \[Omega] + 1 -> -\[Omega]^2,
    1 + \[Omega]^2 -> -\[Omega],
    \[Omega]^2 + 1 -> -\[Omega],
    \[Omega] + \[Omega]^2 -> -1,
    \[Omega]^2 + \[Omega] -> -1,

    (* Cyclotomic relation *)
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

(* ====================================================== *)
(* Vectorization & Matrix Units                           *)
(* ====================================================== *)

vec[X_] := Flatten[Normal[X]];
unvec[v_, n_] := ArrayReshape[v, {n, n}];

e[i_, j_, n_] := SparseArray[{{i, j} -> 1}, {n, n}];
basis[n_] := Flatten[Table[e[i, j, n], {i, n}, {j, n}], 1];

(* ====================================================== *)
(* Operator matrix representation                         *)
(* ====================================================== *)

OperatorMatrix[F_, n_] := Module[{b},
  b = basis[n];
  Transpose[vec /@ (F /@ b)]
];

(* ====================================================== *)
(* Comultiplication & Multiplication                      *)
(* ====================================================== *)

ComultAux[i_, j_, n_] := Sum[KroneckerProduct[e[i, k, n], e[k, j, n]], {k, 1, n}];

Comult[X_] := Module[{n = Length[X]},
  Sum[X[[i, j]] ComultAux[i, j, n], {i, 1, n}, {j, 1, n}]
];

MultOperator[n_] := Module[{b, products},
  b = basis[n];
  products = Flatten[Table[b[[i]] . b[[j]], {i, Length[b]}, {j, Length[b]}], 1];
  Transpose[vec /@ products]
];

Mult = MultOperator[3];

(* ====================================================== *)
(* X matrices                                             *)
(* ====================================================== *)

X0 = {
  {1,0,0,0,1,0,0,0,1}, {0,1,0,0,0,1,1,0,0}, {0,0,1,1,0,0,0,1,0},
  {0,0,1,1,0,0,0,1,0}, {1,0,0,0,1,0,0,0,1}, {0,1,0,0,0,1,1,0,0},
  {0,1,0,0,0,1,1,0,0}, {0,0,1,1,0,0,0,1,0}, {1,0,0,0,1,0,0,0,1}
};

X1 = {
  {1,0,0,0,1,0,0,0,1}, {0,\[Omega]^2,0,0,0,\[Omega]^2,\[Omega]^2,0,0}, {0,0,\[Omega],\[Omega],0,0,0,\[Omega],0},
  {0,0,\[Omega],\[Omega],0,0,0,\[Omega],0}, {1,0,0,0,1,0,0,0,1}, {0,\[Omega]^2,0,0,0,\[Omega]^2,\[Omega]^2,0,0},
  {0,\[Omega]^2,0,0,0,\[Omega]^2,\[Omega]^2,0,0}, {0,0,\[Omega],\[Omega],0,0,0,\[Omega],0}, {1,0,0,0,1,0,0,0,1}
};

X2 = {
  {1,0,0,0,1,0,0,0,1}, {0,\[Omega],0,0,0,\[Omega],\[Omega],0,0}, {0,0,\[Omega]^2,\[Omega]^2,0,0,0,\[Omega]^2,0},
  {0,0,\[Omega]^2,\[Omega]^2,0,0,0,\[Omega]^2,0}, {1,0,0,0,1,0,0,0,1}, {0,\[Omega],0,0,0,\[Omega],\[Omega],0,0},
  {0,\[Omega],0,0,0,\[Omega],\[Omega],0,0}, {0,0,\[Omega]^2,\[Omega]^2,0,0,0,\[Omega]^2,0}, {1,0,0,0,1,0,0,0,1}
};

(* ====================================================== *)
(* Parameterized Operator Function                        *)
(* ====================================================== *)

ComputeFinalOperator[A1_, A2_, A3_] := Module[
  {EOperator, EStarOperator, Aux0, Aux1, Aux2, auxMat0, auxMat1, auxMat2,
   tensorOp0, tensorOp1, tensorOp2, Operator, finalOp},

  EOperator[A_] := 3 Diagonal[A];
  EStarOperator[v_] := DiagonalMatrix[v];

  Aux0[A_] := EStarOperator[A1 . EOperator[A]];
  Aux1[A_] := EStarOperator[A2 . EOperator[A]];
  Aux2[A_] := EStarOperator[A3 . EOperator[A]];

  auxMat0 = OperatorMatrix[Aux0, 3];
  auxMat1 = OperatorMatrix[Aux1, 3];
  auxMat2 = OperatorMatrix[Aux2, 3];

  tensorOp0 = KroneckerProduct[X0, auxMat0];
  tensorOp1 = KroneckerProduct[X1, auxMat1];
  tensorOp2 = KroneckerProduct[X2, auxMat2];

  Operator[A_] := Module[{delta, transformed0, transformed1, transformed2, resultVec},
    delta = Comult[A];
    transformed0 = tensorOp0 . vec[delta];
    transformed1 = tensorOp1 . vec[delta];
    transformed2 = tensorOp2 . vec[delta];
    resultVec = Mult . (transformed0 + transformed1 + transformed2);
    unvec[resultVec/3, 3]
  ];

  finalOp = OperatorMatrix[Operator, 3];
  SimplifyMatrix[finalOp]
];

(* ====================================================== *)
(* Define Triples to Test                                 *)
(* ====================================================== *)

triplesToTest = {

{
  "Only 1- and 2- loops",
  {{0,0,0},
   {0,0,0},
   {0,0,0}},
  {{1,0,0},
   {0,1,0},
   {0,0,1}},
  {{1,0,0},
   {0,1,0},
   {0,0,1}}
},

{
  "1- and 2- cycles along with 1- and 2- loops",
  {{0,0,0},
   {0,0,0},
   {0,0,0}},
  {{1,1,0},
   {0,1,1},
   {1,0,1}},
  {{1,0,1},
   {1,1,0},
   {0,1,1}}
},

{
  "All 1's and all 2's",
  {{0,0,0},
   {0,0,0},
   {0,0,0}},
  {{1,1,1},
   {1,1,1},
   {1,1,1}},
  {{1,1,1},
   {1,1,1},
   {1,1,1}}
},

{
  "Everything except 0-loops",
  {{0,1,1},
   {1,0,1},
   {1,1,0}},
  {{1,1,1},
   {1,1,1},
   {1,1,1}},
  {{1,1,1},
   {1,1,1},
   {1,1,1}}
}

};

(* ====================================================== *)
(* Execution and PDF Export                               *)
(* ====================================================== *)

pdfPages = Table[
  Module[{label, A1, A2, A3, resultMatrix},
    {label, A1, A2, A3} = triple;
    resultMatrix = ComputeFinalOperator[A1, A2, A3];

    Row[{
      Text[Style[label <> " = ", Bold, 12]],
      MatrixForm[resultMatrix]
    }, Alignment -> Center]
  ],
  {triple, triplesToTest}
];

(* Export to PDF silently without standard output *)
Export["FinalGraphs.pdf", pdfPages];
