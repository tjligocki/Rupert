TwoGetSome[n_,digits_] := Module[
  {result,curOut,allOut},
  
  result = {};
  curOut = {{},{}};
  allOut = {};

  Do[
    Print[i];
    curOut = TwoLess[curOut,digits];
    allOut = Append[allOut,curOut];
    result = Append[result,curOut[[3]]]
  ,
    {i,5,n,2}
  ];

  List[result,curOut,allOut]
]
