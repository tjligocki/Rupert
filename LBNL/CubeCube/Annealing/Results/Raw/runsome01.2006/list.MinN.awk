{
  m = $1;
}
END {
  for (n = 2*m; n <= 49; n++) {
    if (n > 2) {
      p = (m-1.0)/(n-2.0);
      p = p*p*p*p*p*p*p*p;
      p = p + 1.0;
    } else {
      p = 1.0;
    }
    printf("%3d %3d %8d %10d\n",m,n,exp(log(1.0*m)*p)*n*(n-1)/2,exp(log(1.0*m)*p)*n*(n-1)/2*n);
  }
}
