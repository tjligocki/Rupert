BEGIN {
  pmax = 0;
  qmax = 0;
}
NF == 3 {
  val[$1 "," $2] = $3;
  if ($1 > pmax) pmax = $1;
  if ($2 > qmax) qmax = $2;
}
END {
  for (p = 1; p <= pmax; p++) {
    for (r = 0; r < p; r++) {
      empty = 1;
      for (q = p+r; q <= qmax; q += p) {
        f = val[p "," q];
        if (f > 0) {
          printf("%20.17lf %20.17lf %20.17lf %3d %3d %3d\n",q/(1.0*p)-1.0,f*f-1.0,f,p,q,r);
          empty = 0;
        }
      }
      if (empty == 0) {
        printf("\n");
      }
    }
  }
}
