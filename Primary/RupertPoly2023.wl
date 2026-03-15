(* ::Package:: *)

BeginPackage["RupertPoly2023`"];

dof::usage = "dof[m,n] - Degrees of freedom for (m,n) problem";

rot::usage = "rot[t,i,j,n] - Rotation matrix by t[i,j] in (i,j) slice of n dimensional Euclidean space";

rot2::usage = "rot2[x,i,j,n] - Rotation matrix by ArcSin[Sqrt[x[i,j]]] in (i,j) slice of n dimensional Euclidean space";

fullRotMat::usage = "fullRotMat[t,m,n] - nxn Orthogonal Matrix for (m,n) problem using angles t[i,j]";
  
fullRotMat2::usage = "fullRotMat2[x,m,n] - nxn Orthogonal Matrix for (m,n) problem using angles ArcSin[Sqrt[x[i,j]]]";

sumRows::usage = "sumRows[mm,m,n] - For matrix mm, sum the first m columns of each row using the absolute value of the entries";

makeEqn::usage = "makeEqn[ms,n] - Make a list of equations using matrix ms";

trySolve::usage = "trySolve[mm,x,xmin,xmax,m,n] - Try to solve the (m,n) problem";

sumRowsNoAbs::usage = "sumRowsNoAbs[mm,m,n] - For matrix mm, sum the first m columns of each row";

trySolveNoAbs::usage = "trySolveNoAbs[mm,x,xmin,xmax,m,n] - Try to solve the (m,n) problem with not absolute values";

makeCon::usage = "makeCon[x,m,n,xmin,xmax] - Generate a list of constraints on the variable for the (m,n) problem";

Begin["`Private`"]

dof[m_, n_] :=
  Simplify[(n (n - 1) / 2) - ((n - m) (n - m - 1) / 2)];

rot[t_, i_, j_, n_] :=
  Module[{mm, tt},
    (
      mm = IdentityMatrix[n];
      ToExpression[StringJoin[ToString[tt], "=", ToString[t], "[", ToString[
        i], ",", ToString[j], "]"]];
      mm[[i, i]] = Cos[tt];
      mm[[i, j]] = -Sin[tt];
      mm[[j, i]] = Sin[tt];
      mm[[j, j]] = Cos[tt];
      mm
    )
  ];

rot2[x_, i_, j_, n_] :=
  Module[{mm, xx},
    (
      mm = IdentityMatrix[n];
      ToExpression[StringJoin[ToString[xx], "=", ToString[x], "[", ToString[
        i], ",", ToString[j], "]"]];
      mm[[i, i]] = Sqrt[1 - xx];
      mm[[i, j]] = -Sqrt[xx];
      mm[[j, i]] = Sqrt[xx];
      mm[[j, j]] = Sqrt[1 - xx];
      mm
    )
  ];

fullRotMat[t_, m_, n_] :=
  Module[{mm},
    (
      mm = IdentityMatrix[n];
      For[i = 1, i <= n, i++,
        For[j = i + 1, j <= n, j++,
          If[i <= m || j <= m,
            mm = mm . rot[t, i, j, n]
            ,
            Nothing
          ]
        ]
      ];
      mm
    )
  ]

fullRotMat2[x_, m_, n_] :=
  Module[{mm},
    (
      mm = IdentityMatrix[n];
      For[i = 1, i <= n, i++,
        For[j = i + 1, j <= n, j++,
          If[i <= m || j <= m,
            mm = mm . rot2[x, i, j, n]
            ,
            Nothing
          ]
        ]
      ];
      mm
    )
  ]

fullRotMat2[t, 1, 3] // MatrixForm

sumRows[mm_, m_, n_] :=
  Table[Sum[Map[Abs, mm][[j, i]], {i, 1, m}], {j, 1, n}]

makeEqn[ms_, n_] :=
  Table[p == ms[[i]], {i, 1, n}]

trySolve[mm_, x_, xmin_, xmax_, m_, n_] :=
  Module[{mmm, mms, mme},
    (
      mmm = mm;
      mms = sumRows[mmm, m, n];
      mme = makeEqn[mms, n];
      Solve[Join[mme, {0 <= p <= 1}, makeCon[x, m, n, xmin, xmax]], Reals
        ]
    )
  ]

sumRowsNoAbs[mm_, m_, n_] :=
  Table[Sum[mm[[j, i]], {i, 1, m}], {j, 1, n}]

trySolveNoAbs[mm_, x_, xmin_, xmax_, m_, n_] :=
  Module[{mmm, mms, mme},
    (
      mmm = mm;
      mms = sumRowsNoAbs[mmm, m, n];
      mme = makeEqn[mms, n];
      Solve[Join[mme, {0 <= p <= 1}, makeCon[x, m, n, xmin, xmax]], Reals
        ]
    )
  ]

makeCon[x_, m_, n_, xmin_, xmax_] :=
  Flatten[Table[Table[{xmin <= x[i, j] <= xmax}, {j, i + 1, n}], {i, 
    1, m}]]

End[];

EndPackage[];
