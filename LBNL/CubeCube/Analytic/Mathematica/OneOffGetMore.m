OneOffGetMore[n_,space_,prev_,digits_] := Module[
  {result,curOut,allOut,num},
  
  result = prev[[1]];
  num = Length[result] + 2;

  curOut = prev[[2]];
  allOut = prev[[3]];

  Do[
    Print["("<>ToString[i-space]<>","<>ToString[i]<>")"];
    curOut = OneOff[curOut,space,digits];
    allOut = Append[allOut,curOut];
    result = Append[result,curOut[[3]]]
  ,
    {i,num*space+1,n,space}
  ];

  List[result,curOut,allOut]
]
