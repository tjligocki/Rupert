BEGIN {
  cur = 0;
  num = 0;
}
NR == 1 {
  for (i = 1; i <= NF; i++) {
    x0[i] = $i;
  }

  num = NF;
}
NR > 1 {
  for (i = 1; i <= num; i++) {
    s[cur,i] = $i - x0[i];
  }

  cur++;
}
END {
  for (i = 0; i < cur; i++) {
    for (j = 0; j < cur; j++) {
      dot = 0.0
      for (k = 1; k <= num; k++) {
        dot += s[i,k]*s[j,k];
      }
      printf("%2d %2d: %20.17lf\n",i,j,dot);
    }
  }
}
