END {
  for (n = 1; n <= 49; n++) {
    for (m = 1; m < n; m++) {
      good = 1;
      for (k = 2; k <= m; k++) {
        if (n % k == 0 && m % k == 0) {
          good = 0;
          break;
        }
      }
      if (good == 1) {
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
  }
}
