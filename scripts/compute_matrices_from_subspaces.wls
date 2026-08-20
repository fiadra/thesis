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

(* Full simplification forcing reduction modulo
   (1 + omega + omega^2) *)
SimplifyMatrix[M_] := Module[{reduced},

    reduced = M /. omegaRules;

    reduced = Map[
        Function[elem,
            Last[
                PolynomialReduce[
                    elem,
                    {1 + \[Omega] + \[Omega]^2},
                    {\[Omega]}
                ]
            ]
        ],
        reduced,
        {-1}
    ];

    FullSimplify[reduced /. omegaRules] /. omegaRules
];


(* ====================================================== *)
(* Display                                                *)
(* ====================================================== *)

(* Display matrices in simplified form.
   The zero matrix is left unchanged. *)
FactorThird[M_] := If[
    M === ConstantArray[0, {9, 9}],
    M,
    MatrixForm[SimplifyMatrix[M]]
];


(* ====================================================== *)
(* Matrix generators P (Shift) and Q (Clock)              *)
(* ====================================================== *)

P = {
    {1, 0, 0},
    {0, \[Omega], 0},
    {0, 0, \[Omega]^2}
};

Q = {
    {0, 0, 1},
    {1, 0, 0},
    {0, 1, 0}
};


(* ====================================================== *)
(* P_k Sets                                               *)
(* ====================================================== *)

P1 = {
    P,
    P . P
};

P2 = {
    Q,
    Q . Q
};

P3 = {
    P . Q,
    (P . P) . (Q . Q)
};

P4 = {
    (P . P) . Q,
    P . (Q . Q)
};

blocksSets = {P1, P2, P3, P4};


(* ====================================================== *)
(* Tensor Product Block Construction                      *)
(* ====================================================== *)

(* MODIFIED: Now computes G x G* + H x H* to match the
   conjugate-tensor computation of the second script. *)
TensorBlock[{G_, H_}] := Module[{},
    SimplifyMatrix[
        KroneckerProduct[G, Conjugate[G]] +
        KroneckerProduct[H, Conjugate[H]]
    ]
];

blocks = TensorBlock /@ blocksSets;


(* ====================================================== *)
(* Total sl_3 Block                                       *)
(* ====================================================== *)

Psl3 = SimplifyMatrix[Total[blocks]];


(* ====================================================== *)
(* Subspace Aggregation                                   *)
(* ====================================================== *)

allSubspaces = Join[

    {
        {
            "m = 0 : V_0",
            ConstantArray[0, {9, 9}]
        }
    },

    {
        {"m = 2 : P_1", blocks[[1]]},
        {"m = 2 : P_2", blocks[[2]]},
        {"m = 2 : P_3", blocks[[3]]},
        {"m = 2 : P_4", blocks[[4]]}
    },

    Table[
        {
            "m = 4 : W_" <>
                ToString[pair[[1]]] <>
                ToString[pair[[2]]],

            SimplifyMatrix[
                blocks[[pair[[1]]]] +
                blocks[[pair[[2]]]]
            ]
        },
        {pair, Subsets[Range[4], {2}]}
    ],

    Table[
        {
            "m = 6 : U_" <> ToString[k],

            SimplifyMatrix[
                Psl3 - blocks[[k]]
            ]
        },
        {k, 1, 4}
    ],

    {
        {
            "m = 8 : sl_3",
            Psl3
        }
    }
];


(* ====================================================== *)
(* Layout and PDF Export                                  *)
(* ====================================================== *)

pdfPages = Table[
    Row[
        {
            Text[
                Style[
                    item[[1]] <> "  =  ",
                    Bold,
                    12
                ]
            ],

            If[
                MatrixQ[item[[2]]],
                FactorThird[item[[2]]],
                MatrixForm[item[[2]]]
            ]
        },
        Alignment -> Center
    ],
    {item, allSubspaces}
];

Export[
    "SelfReflectiveSubspaces.pdf",
    pdfPages
];
