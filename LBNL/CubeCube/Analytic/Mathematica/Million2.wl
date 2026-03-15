(* ::Package:: *)

TwoGetSome[n_,digits_] := Module[
  {result,curOut,allOut},
  
  result = {};
  curOut = {{},{}};
  allOut = {};

  Do[
    Print[i];
    myOut = Timing[TwoLess[curOut,digits]];
	Print[myOut[[1]]];
	curOut=myOut[[2]];
    allOut = Append[allOut,curOut];
    result = Append[result,curOut[[3]]]
  ,
    {i,5,n,2}
  ];

  List[result,curOut,allOut]
]


t = Timing[TwoGetSome[20,50000]];


Export["e2_0050000.dat",t[[2]][[1]]^2]


t = Timing[TwoGetSome[20,100000]];


Export["e2_0100000.dat",t[[2]][[1]]^2]


t = Timing[TwoGetSome[20,200000]];


Export["e2_0200000.dat",t[[2]][[1]]^2]


t = Timing[TwoGetSome[20,500000]];


Export["e2_0500000.dat",t[[2]][[1]]^2]


t = Timing[TwoGetSome[20,1000000]];


Export["e2_1000000.dat",t[[2]][[1]]^2]


t = Timing[TwoGetSome[20,2000000]];


Export["e2_2000000.dat",t[[2]][[1]]^2]
