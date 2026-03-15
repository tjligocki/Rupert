OutputAll[valList_,file_,digits_] := Module[
  {outList,cur,pos},

  outList = {};
  Do[
    cur = ToString[valList[[i]]];
    pos = StringPosition[cur,"."][[1]][[1]];
    cur = StringTake[cur,pos+digits];
    outList = Append[outList,cur];
  ,
    {i,1,Length[valList]}
  ];

  Export[file,outList,"List"]
]
