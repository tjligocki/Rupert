OneInN[n_,digits_] := Module[
  {curOut},
  
  allOut = {};

  Do[
    Print[i];
    curOut = SetPrecision[Sqrt[i],digits];
    allOut = Append[allOut,curOut];
  ,
    {i,2,n}
  ];

  allOut
]
