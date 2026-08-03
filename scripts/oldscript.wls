#!/usr/bin/env wolframscript

ClearAll["Global`*"];

(* ====================================================== *)
(* Vectorization                                          *)
(* ====================================================== *)

vec[X_] := Flatten[Normal[X]];

unvec[v_, n_] := ArrayReshape[v, {n, n}];

(* ====================================================== *)
(* Matrix units                                           *)
(* ====================================================== *)

e[i_, j_, n_] :=
  SparseArray[
    {{i, j} -> 1},
    {n, n}
  ];

(* ====================================================== *)
(* Basis of M_n                                           *)
(* ====================================================== *)

basis[n_] :=
  Flatten[
    Table[
      e[i, j, n],
      {i, n},
      {j, n}
    ],
    1
  ];

(* ====================================================== *)
(* Operator matrix representation                         *)
(* ====================================================== *)

OperatorMatrix[F_, n_] :=
 Module[
   {b},

   b = basis[n];

   Transpose[
     vec /@ (F /@ b)
   ]
 ];

(* ====================================================== *)
(* Comultiplication                                       *)
(* ====================================================== *)

ComultAux[i_, j_, n_] :=
  Sum[
    KroneckerProduct[
      e[i, k, n],
      e[k, j, n]
    ],
    {k, 1, n}
  ];

Comult[X_] :=
 Module[
   {n = Length[X]},

   Sum[
     X[[i, j]] ComultAux[i, j, n],
     {i, 1, n},
     {j, 1, n}
   ]
 ];

(* ====================================================== *)
(* Multiplication operator                                *)
(* ====================================================== *)

MultOperator[n_] :=
 Module[
   {b, products},

   b = basis[n];

   products =
     Flatten[
       Table[
         b[[i]] . b[[j]],
         {i, Length[b]},
         {j, Length[b]}
       ],
       1
     ];

   Transpose[
     vec /@ products
   ]
 ];

Mult = MultOperator[3];

(* ====================================================== *)
(* Cube root of unity                                     *)
(* ====================================================== *)

Clear[\[Omega]];

(* ====================================================== *)
(* X matrices                                             *)
(* ====================================================== *)

X0 = {
  {1,0,0,0,1,0,0,0,1},
  {0,1,0,0,0,1,1,0,0},
  {0,0,1,1,0,0,0,1,0},
  {0,0,1,1,0,0,0,1,0},
  {1,0,0,0,1,0,0,0,1},
  {0,1,0,0,0,1,1,0,0},
  {0,1,0,0,0,1,1,0,0},
  {0,0,1,1,0,0,0,1,0},
  {1,0,0,0,1,0,0,0,1}
};

X1 = {
  {1,0,0,0,1,0,0,0,1},
  {0,\[Omega]^2,0,0,0,\[Omega]^2,\[Omega]^2,0,0},
  {0,0,\[Omega],\[Omega],0,0,0,\[Omega],0},
  {0,0,\[Omega],\[Omega],0,0,0,\[Omega],0},
  {1,0,0,0,1,0,0,0,1},
  {0,\[Omega]^2,0,0,0,\[Omega]^2,\[Omega]^2,0,0},
  {0,\[Omega]^2,0,0,0,\[Omega]^2,\[Omega]^2,0,0},
  {0,0,\[Omega],\[Omega],0,0,0,\[Omega],0},
  {1,0,0,0,1,0,0,0,1}
};

X2 = {
  {1,0,0,0,1,0,0,0,1},
  {0,\[Omega],0,0,0,\[Omega],\[Omega],0,0},
  {0,0,\[Omega]^2,\[Omega]^2,0,0,0,\[Omega]^2,0},
  {0,0,\[Omega]^2,\[Omega]^2,0,0,0,\[Omega]^2,0},
  {1,0,0,0,1,0,0,0,1},
  {0,\[Omega],0,0,0,\[Omega],\[Omega],0,0},
  {0,\[Omega],0,0,0,\[Omega],\[Omega],0,0},
  {0,0,\[Omega]^2,\[Omega]^2,0,0,0,\[Omega]^2,0},
  {1,0,0,0,1,0,0,0,1}
};

(* ====================================================== *)
(* Circulant matrices                                     *)
(* ====================================================== *)

# A1 = {
#    {a1,c1,b1},
#    {b1,a1,c1},
#    {c1,b1,a1}
# };
A1 = {
   {0,0,0},
   {0,0,0},
   {0,0,0}
};


A2 = {
   {1,0,0},
   {0,1,0},
   {0,0,1}
};

A3 = {
   {1,0,0},
   {0,1,0},
   {0,0,1}
};

(* ====================================================== *)
(* E and E*                                               *)
(* ====================================================== *)

EOperator[A_] := 3 Diagonal[A];

EStarOperator[v_] := DiagonalMatrix[v];

(* ====================================================== *)
(* Auxiliary maps                                         *)
(* ====================================================== *)

Aux0[A_] :=
  EStarOperator[
    A1 . EOperator[A]
  ];

Aux1[A_] :=
  EStarOperator[
    A2 . EOperator[A]
  ];

Aux2[A_] :=
  EStarOperator[
    A3 . EOperator[A]
  ];

(* ====================================================== *)
(* Matrix representations of Aux maps                     *)
(* ====================================================== *)

auxMat0 = OperatorMatrix[Aux0, 3];
auxMat1 = OperatorMatrix[Aux1, 3];
auxMat2 = OperatorMatrix[Aux2, 3];

(* ====================================================== *)
(* Tensor operators                                       *)
(* ====================================================== *)

tensorOp0 = KroneckerProduct[X0, auxMat0];
tensorOp1 = KroneckerProduct[X1, auxMat1];
tensorOp2 = KroneckerProduct[X2, auxMat2];

(* ====================================================== *)
(* Main operator                                          *)
(* ====================================================== *)

Operator[A_] :=
 Module[
   {
     delta,
     transformed0,
     transformed1,
     transformed2,
     resultVec
   },

   delta = Comult[A];

   transformed0 = tensorOp0 . vec[delta];
   transformed1 = tensorOp1 . vec[delta];
   transformed2 = tensorOp2 . vec[delta];

   resultVec =
     Mult . (
       transformed0 +
       transformed1 +
       transformed2
     );

   unvec[resultVec/3, 3]
 ];

(* ====================================================== *)
(* Final operator matrix with Symbolic \[Omega] Logic     *)
(* ====================================================== *)

(* Define the polynomial relation for the cube root of unity *)
omegaRelation = {\[Omega]^2 + \[Omega] + 1 -> 0, \[Omega]^3 -> 1};

FinalOperator =
  OperatorMatrix[Operator, 3] // Expand;

(* Apply the simplification rules specifically for omega *)
FinalOperator = Simplify[FinalOperator /. omegaRelation];

(* Export the MatrixForm to the PDF for better visual formatting *)
Export["FinalOperator.pdf", MatrixForm[FinalOperator]];

Print["Final operator matrix:"];

(* Use MatrixForm for the console output *)
Print[
  MatrixForm[FinalOperator]
];
