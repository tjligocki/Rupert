OneGetMore[n_,prev_,digits_] := Module[
  {result,curOut,allOut,num},
  
  result = prev[[1]];
  num = Length[result] + 1;

  curOut = prev[[2]];
  allOut = prev[[3]];

  Do[
    Print[i + num];
    curOut = OneLess[curOut,digits];
    allOut = Append[allOut,curOut];
    result = Append[result,curOut[[3]]]
  ,
    {i,1,n}
  ];

  List[result,curOut,allOut]
]
