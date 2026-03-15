BEGIN {
  m = 0;
  n = 0;
  nz = 0;
}
NR == 1 {
  m = $1;
  n = $2;
  next;
}
m > 0 && NR <= n+1 {
  for (i = 1; i <= m; i++) {
    if ($i == 0) {
      nz++;
    }
  }
}
END {
  printf("%2d %2d -> %4d - %2d %3d - %4d\n",m,n,(m*(2*n-m-1))/2,n-1,nz,n-1+nz);
}
