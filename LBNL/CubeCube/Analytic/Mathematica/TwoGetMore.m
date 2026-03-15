TwoGetMore[n_,prev_,digits_] := Module[
  {result,curOut,allOut,num},
  
  result = prev[[1]];
  num = Length[result] + 1;

  curOut = prev[[2]];
  allOut = prev[[3]];

  Do[
    Print[i];
    curOut = TwoLess[curOut,digits];
    allOut = Append[allOut,curOut];
    result = Append[result,curOut[[3]]]
  ,
    {i,2*num+1,n,2}
  ];

  List[result,curOut,allOut]
]
