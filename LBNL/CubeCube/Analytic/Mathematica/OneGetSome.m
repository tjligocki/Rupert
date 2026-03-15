OneGetSome[n_,digits_] := Module[
  {result,curOut,allOut},
  
  result = {};
  curOut = {{},{}};
  allOut = {};

  Do[
    Print[i];
    curOut = OneLess[curOut,digits];
    allOut = Append[allOut,curOut];
    result = Append[result,curOut[[3]]]
  ,
    {i,2,n}
  ];

  List[result,curOut,allOut]
]
