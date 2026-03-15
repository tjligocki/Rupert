OneLess[eqn_,digits_] := Module[
  {eqn1Len,eqn2Len,
   newVar,newTerm,
   eqnPlus,initPlus,initLast,
   result,
   initNew,sideInfo,tmp},

  eqn1Len = Length[eqn[[1]]];
  eqn2Len = Length[eqn[[2]]];

  newVar = ToExpression["t" <> ToString[eqn1Len+1]];
  eqnPlus = {newVar^2 == 1/2};

  If [eqn1Len > 0,
    newTerm = ToString[newVar];
    Do[
      newTerm = newTerm <> " + Sqrt[1 - " <> ToString[newVar] <> "^2] ";
      Do[
        newTerm = newTerm <> "t" <> ToString[j] <> " ",
        {j,eqn1Len,i+1,-1}
      ];
      newTerm = newTerm <> "Sqrt[1 - t" <> ToString[i] <> "^2]",
      {i,eqn1Len,1,-1}
    ];
    newTerm = newTerm <> " == t1";
    eqnPlus = Append[eqnPlus,ToExpression[newTerm]];

    eqnPlus = Join[eqnPlus,Rest[eqn[[1]]]];
  ];

  initPlus = {{newVar, Sqrt[2]/2}};

  If [eqn2Len > 0,
    Do[
      initPlus = Append[initPlus,
                        {eqn[[2]][[i]][[1]],eqn[[2]][[i+1]][[2]]}
                 ],
      {i,1,eqn2Len-1}
    ];

    initLast = eqn[[2]][[eqn2Len]][[2]];
    initPlus = Append[initPlus,
                      {eqn[[2]][[eqn2Len]][[1]],9*initLast / (1 + 8*initLast)}
               ];
  ];

  result = FindRoot[eqnPlus,initPlus,WorkingPrecision->digits];

  initNew = {};

  Do[
    tmp = {};
    Do[
      AppendTo[tmp,result[[i]][[j]]],
      {j,1,Length[result[[i]]]}
    ];
    AppendTo[initNew,tmp],
    {i,1,Length[result]}
  ];

  sideInfo = N[1 / result[[Length[result]]][[2]],digits];

  List[eqnPlus,initNew,sideInfo]
]
