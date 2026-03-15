BEGIN {
  findit = 0;
}
NF == 2 {
  m = $1;
  n = $2;
  findit = 1;
  next;
}
$1 == "Error:" && findit == 1 {
  findit = 2;
  next;
}
$1 == "Best" && findit == 2 {
  val = $3;
  val = substr(val,2,length(val)-5);
  printf("%2d %2d   %s\n",m,n,val);
  findit = 0;
  next;
}
