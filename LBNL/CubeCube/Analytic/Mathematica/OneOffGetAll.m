OneOffGetAll[n_,initSpace_,digits_] := Module[
  {result,curOut,allOut,someOut,space},

  result = {};
  curOut = {};
  allOut = {};
  
  Do[
    someOut = OneOffGetSome[n,space,digits];

    result = Join[result,someOut[[1]]];
    curOut = Join[curOut,someOut[[2]]];
    allOut = Join[allOut,someOut[[3]]];
  ,
    {space,initSpace,n/2}
  ];

  List[result,curOut,allOut]
]
