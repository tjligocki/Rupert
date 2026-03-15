TwoInN[n_,digits_] := Module[
  {curOut},
  
  allOut = {};

  Do[
    Print[i];
    curOut = SetPrecision[Sqrt[(i-1)/2 + 1/8],digits];
    allOut = Append[allOut,curOut];
  ,
    {i,3,n,2}
  ];

  allOut
]
