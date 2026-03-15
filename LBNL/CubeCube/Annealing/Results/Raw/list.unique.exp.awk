END {
  for (m = 1; m <= 25; m++) {
    for (n = m+1; n <= 25; n++) {
      good = 1;
      for (k = 2; k <= m; k++) {
        if (n % k == 0 && m % k == 0) {
          good = 0;
          break;
        }
      }
      if (good == 1) {
        mul = exp(m*log(300.0))*n*(n-1)/2.0/30.0/5.0;
        printf("%3d %3d %63.0lf %63.0lf\n",m,n,mul,mul*n);
      }
    }
  }
}
