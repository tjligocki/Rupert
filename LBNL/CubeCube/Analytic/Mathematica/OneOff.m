OneOff[eqn_,space_,digits_] := Module[
  { eqnLocal,
    eqn1Len,eqn2Len,
    newVar,newTerm,
    eqnPlus,initPlus,initLast,
    result,
    initNew,sideInfo,tmp
  }
  ,
  eqn1Len = Length[eqn[[1]]];
  eqn2Len = Length[eqn[[2]]];

  If [eqn1Len == 0,
    newVar = ToExpression["t" <> ToString[eqn1Len+1]];
    eqnPlus = {newVar^2 == 1/2};

    initPlus = {{newVar, Sqrt[2]/2}};

    result = FindRoot[eqnPlus,initPlus,WorkingPrecision->digits];

    initNew = {};

    Do[
      tmp = {};
      Do[
        AppendTo[tmp,result[[i]][[j]]]
      ,
        {j,1,Length[result[[i]]]}
      ];
      AppendTo[initNew,tmp]
    ,
      {i,1,Length[result]}
    ];

    sideInfo = N[1 / result[[Length[result]]][[2]],digits];

    eqnLocal = List[eqnPlus,initNew,sideInfo];
  ,
    eqnLocal = eqn;
  ];

  eqn1Len = Length[eqnLocal[[1]]];
  eqn2Len = Length[eqnLocal[[2]]];

  newVar = ToExpression["t" <> ToString[eqn1Len+1]];
  eqnPlus = {newVar^2 == 1/2};

  newTerm = ToString[newVar];
  Do[
    newTerm = newTerm <> " + Sqrt[1 - " <> ToString[newVar] <> "^2] ";

    Do[
      newTerm = newTerm <> "t" <> ToString[j] <> " "
    ,
      {j,eqn1Len,i+1,-1}
    ];

    newTerm = newTerm <> "Sqrt[1 - t" <> ToString[i] <> "^2]";

    If [i == 1,
      newTerm = newTerm <> " / Sqrt[" <> ToString[space] <> "]";
    ]
  ,
    {i,eqn1Len,1,-1}
  ];
  newTerm = newTerm <> " == t1";
  eqnPlus = Append[eqnPlus,ToExpression[newTerm]];

  eqnPlus = Join[eqnPlus,Rest[eqnLocal[[1]]]];

  initPlus = {{newVar, Sqrt[2]/2}};

  Do[
    initPlus = Append[initPlus,
                      {eqnLocal[[2]][[i]][[1]],eqnLocal[[2]][[i+1]][[2]]}
               ]
  ,
    {i,1,eqn2Len-1}
  ];

  initLast = eqnLocal[[2]][[eqn2Len]][[2]];
  initPlus = Append[initPlus,
                    {eqnLocal[[2]][[eqn2Len]][[1]],9*initLast / (1 + 8*initLast)}
             ];

  result = FindRoot[eqnPlus,initPlus,WorkingPrecision->digits];

  initNew = {};

  Do[
    tmp = {};
    Do[
      AppendTo[tmp,result[[i]][[j]]]
    ,
      {j,1,Length[result[[i]]]}
    ];
    AppendTo[initNew,tmp]
  ,
    {i,1,Length[result]}
  ];

  sideInfo = N[1 / result[[Length[result]]][[2]],digits];

  List[eqnPlus,initNew,sideInfo]
]
