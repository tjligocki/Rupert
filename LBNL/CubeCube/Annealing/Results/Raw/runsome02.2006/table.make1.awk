BEGIN {
  findit = 0;
}
NF == 2 {
  m = $1;
  n = $2;
  findit = 1;
}
$1 == "Best" && findit == 1 {
  val = $3;
  val = substr(val,2,length(val)-5);
  printf("%2d %2d   %s\n",m,n,val);
  findit = 0;
}
