#!/usr/bin/env wolframscript

(*
  verify_voltage_graphs.wl

  Checks that the four "derived graph" figures (Gamma^lambda) are correct
  lifts of their Z3-voltage base graphs (Gamma), using TWO independent
  constructions that must agree:

    (A) DeriveGraph        -- the general voltage-graph lifting formula
                               (a,i) -- (b, i+g mod n) for every base
                               edge/loop {a,b,g}.

    (B) AsDrawnGraph        -- a literal transcription of the \drawlift
                               and \drawderivedloops TikZ macros from the
                               .tex file (the exact pairs of fiber
                               coordinates that get connected on the page).

  If (A) and (B) produce the same edge set for every case, the figure is a
  faithful drawing of the voltage graph it claims to depict. On top of
  that, each derived graph is checked against the specific structural
  claims implied by the figure captions (edge counts, regularity, K9,
  disjoint triangles, etc).

  Run with:   wolframscript -file verify_voltage_graphs.wl
*)

n = 3; (* Z_3 voltages throughout *)

(* ------------------------------------------------------------------ *)
(* (A) General voltage-graph lift, from the {a, b, g} edge list       *)
(* ------------------------------------------------------------------ *)

(* baseEdges: list of {a, b, g}. a=b encodes a loop. g is the voltage,
   taken in Z_n. Returns a canonical list of undirected {v1, v2} pairs
   (v1, v2 sorted) with no self-loops and no duplicates. *)
DeriveEdgePairs[baseEdges_List, nn_Integer] := Module[{pairs},
  pairs = Flatten[
    Table[
      With[{a = e[[1]], b = e[[2]], g = e[[3]]},
        Table[
          Sort[{{a, i}, {b, Mod[i + g, nn]}}],
          {i, 0, nn - 1}
        ]
      ],
      {e, baseEdges}
    ],
    1
  ];
  pairs = Select[pairs, First[#] =!= Last[#] &];  (* g=0 loops would be degenerate; none occur here *)
  DeleteDuplicates[pairs]
]

DeriveGraph[baseEdges_List, nn_Integer] :=
  Graph[UndirectedEdge @@@ DeriveEdgePairs[baseEdges, nn],
    VertexLabels -> "Name", GraphLayout -> "SpringElectricalEmbedding"]

(* ------------------------------------------------------------------ *)
(* (B) Literal transcription of the TikZ macros                       *)
(* ------------------------------------------------------------------ *)

(* \drawderivedloops: for every base vertex a, connects
     (a,0)-(a,1), (a,1)-(a,2), (a,0)-(a,2)     -- always drawn, all 4 cases *)
DrawnTriangle[a_] :=
  Sort /@ {{{a, 0}, {a, 1}}, {{a, 1}, {a, 2}}, {{a, 0}, {a, 2}}}

(* \drawlift{a}{b}{g}: hand-transcribed from the three \ifnum branches
   in the .tex source (NOT reusing the Mod formula, so this is an
   independent check on the macro's literal wiring). *)
DrawnLift[a_, b_, 0] := Sort /@ {
   {{a, 0}, {b, 0}}, {{a, 1}, {b, 1}}, {{a, 2}, {b, 2}}}
DrawnLift[a_, b_, 1] := Sort /@ {
   {{a, 0}, {b, 1}}, {{a, 1}, {b, 2}}, {{a, 2}, {b, 0}}}
DrawnLift[a_, b_, 2] := Sort /@ {
   {{a, 0}, {b, 2}}, {{a, 1}, {b, 0}}, {{a, 2}, {b, 1}}}

AsDrawnEdgePairs[baseVerts_List, liftCalls : {{_, _, _} ...}] :=
  DeleteDuplicates[
    Join[
      Flatten[DrawnTriangle /@ baseVerts, 1],
      Flatten[DrawnLift @@@ liftCalls, 1]
    ]
  ]

AsDrawnGraph[baseVerts_List, liftCalls_List] :=
  Graph[UndirectedEdge @@@ AsDrawnEdgePairs[baseVerts, liftCalls],
    VertexLabels -> "Name"]

(* ------------------------------------------------------------------ *)
(* The four cases, exactly as labeled in the figures                  *)
(* ------------------------------------------------------------------ *)

baseVerts = {"u", "v", "w"};
loopEdges = {{"u", "u", 1}, {"u", "u", 2}, {"v", "v", 1}, {"v", "v", 2},
             {"w", "w", 1}, {"w", "w", 2}};   (* loops labeled "1,2" at every vertex, all cases *)

case1Base = loopEdges;
case2Base = Join[loopEdges, {{"u", "v", 1}, {"v", "w", 1}, {"w", "u", 1}}];
case3Base = Join[loopEdges, {{"u", "v", 1}, {"u", "v", 2}, {"v", "w", 1},
                              {"v", "w", 2}, {"w", "u", 1}, {"w", "u", 2}}];
case4Base = Join[loopEdges, {{"u", "v", 0}, {"u", "v", 1}, {"u", "v", 2},
                              {"v", "w", 0}, {"v", "w", 1}, {"v", "w", 2},
                              {"w", "u", 0}, {"w", "u", 1}, {"w", "u", 2}}];

case1Lifts = {};                                            (* Figure 2: no cross-fiber lifts *)
case2Lifts = {{"u", "v", 1}, {"v", "w", 1}, {"w", "u", 1}}; (* Figure 4 *)
case3Lifts = {{"u", "v", 1}, {"u", "v", 2}, {"v", "w", 1},
              {"v", "w", 2}, {"w", "u", 1}, {"w", "u", 2}}; (* Figure 6 *)
case4Lifts = {{"u", "v", 0}, {"u", "v", 1}, {"u", "v", 2},
              {"v", "w", 0}, {"v", "w", 1}, {"v", "w", 2},
              {"w", "u", 0}, {"w", "u", 1}, {"w", "u", 2}}; (* Figure 8 *)

cases = <|
  1 -> <|"caption" -> "Figure 2: undirected 3-cycles",
         "base" -> case1Base, "lifts" -> case1Lifts,
         "expectedEdges" -> 9, "expectedDegree" -> 2|>,
  2 -> <|"caption" -> "Figure 4: exactly 18 undirected edges",
         "base" -> case2Base, "lifts" -> case2Lifts,
         "expectedEdges" -> 18, "expectedDegree" -> 4|>,
  3 -> <|"caption" -> "Figure 6: partial K33 between fibers",
         "base" -> case3Base, "lifts" -> case3Lifts,
         "expectedEdges" -> 27, "expectedDegree" -> 6|>,
  4 -> <|"caption" -> "Figure 8: perfect K33 between all fibers",
         "base" -> case4Base, "lifts" -> case4Lifts,
         "expectedEdges" -> 36, "expectedDegree" -> 8|>
|>;

(* ------------------------------------------------------------------ *)
(* Helper: induced bipartite graph between two fibers, cross edges only *)
(* ------------------------------------------------------------------ *)

CrossFiberGraph[g_Graph, fiberA_, fiberB_] := Module[{es},
  es = Select[
    EdgeList[g],
    (First[#[[1]]] === fiberA && First[#[[2]]] === fiberB) ||
    (First[#[[1]]] === fiberB && First[#[[2]]] === fiberA) &
  ];
  Graph[es]
]

(* ------------------------------------------------------------------ *)
(* Run all checks                                                      *)
(* ------------------------------------------------------------------ *)

results = <||>;
allPassed = True;

Do[
  Module[{c, gFormula, gDrawn, pairsFormula, pairsDrawn, sameEdges,
          nVerts, nEdges, degrees, regular, checks, crossUV},

    c = cases[k];
    gFormula = DeriveGraph[c["base"], n];
    gDrawn = AsDrawnGraph[baseVerts, c["lifts"]];

    pairsFormula = Sort[DeriveEdgePairs[c["base"], n]];
    pairsDrawn = Sort[AsDrawnEdgePairs[baseVerts, c["lifts"]]];
    sameEdges = (pairsFormula === pairsDrawn);

    nVerts = VertexCount[gFormula];
    nEdges = EdgeCount[gFormula];
    degrees = DeleteDuplicates[VertexDegree[gFormula]];
    regular = (Length[degrees] == 1);

    checks = <|
      "Formula matches TikZ macros (same edge set)" -> sameEdges,
      "9 vertices (3 fibers x 3)" -> (nVerts == 9),
      "Edge count = " <> ToString[c["expectedEdges"]] -> (nEdges == c["expectedEdges"]),
      "Regular of degree " <> ToString[c["expectedDegree"]] ->
        (regular && degrees == {c["expectedDegree"]})
    |>;

    (* Case-specific structural fingerprints *)
    Switch[k,
      1,
      checks["3 disjoint triangles (C3 + C3 + C3)"] =
        IsomorphicGraphQ[gFormula,
          GraphDisjointUnion[CycleGraph[3], CycleGraph[3], CycleGraph[3]]];
      ,
      2,
      crossUV = CrossFiberGraph[gFormula, "u", "v"];
      checks["Cross-fiber (u,v) subgraph is a perfect matching"] =
        (EdgeCount[crossUV] == 3 &&
         DeleteDuplicates[VertexDegree[crossUV]] == {1});
      ,
      3,
      crossUV = CrossFiberGraph[gFormula, "u", "v"];
      checks["Cross-fiber (u,v) subgraph is K33 minus a perfect matching (\[TildeTilde] C6)"] =
        IsomorphicGraphQ[crossUV, CycleGraph[6]];
      checks["Complement of derived graph (in K9) is 3 disjoint triangles"] =
        IsomorphicGraphQ[
          GraphComplement[gFormula],
          GraphDisjointUnion[CycleGraph[3], CycleGraph[3], CycleGraph[3]]
        ];
      ,
      4,
      checks["Derived graph is the complete graph K9"] =
        IsomorphicGraphQ[gFormula, CompleteGraph[9]];
    ];

    results[k] = <|"caption" -> c["caption"], "checks" -> checks|>;
    If[Not[And @@ Values[checks]], allPassed = False];
  ],
  {k, 1, 4}
];

(* ------------------------------------------------------------------ *)
(* Report                                                              *)
(* ------------------------------------------------------------------ *)

Print["Voltage graph verification (Z_3, base triangle u-v-w)\n"];

Do[
  Module[{entry = results[k]},
    Print["Case ", k, "  (", entry["caption"], ")"];
    KeyValueMap[
      Print["    [", If[#2, "PASS", "FAIL"], "]  ", #1] &,
      entry["checks"]
    ];
    Print[];
  ],
  {k, 1, 4}
];

Print[If[allPassed,
  "ALL CHECKS PASSED \[Dash] every derived graph matches its voltage assignment, \
and every figure matches the general lifting formula.",
  "SOME CHECKS FAILED \[Dash] see above."
]];
