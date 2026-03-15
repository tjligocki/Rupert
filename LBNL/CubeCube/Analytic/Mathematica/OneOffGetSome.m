OneOffGetSome[n_,space_,digits_] := Module[
  {result,curOut,allOut},
  
  result = {};
  curOut = {{},{}};
  allOut = {};

  Do[
    Print["("<>ToString[i-space]<>","<>ToString[i]<>")"];
    curOut = OneOff[curOut,space,digits];
    allOut = Append[allOut,curOut];
    result = Append[result,curOut[[3]]]
  ,
    {i,2*space+1,n,space}
  ];

  List[result,curOut,allOut]
]
