BEGIN {
  cur = 0;
}
NR == 1 {
  for (i = 1; i <= NF; i++) {
    x[cur,i] = $i;
  }
  num = NF;

  cur++;

  looking = 2;
  factor = 1;
}
NR == looking {
  for (i = 1; i <= num; i++) {
    x[cur,i] = $i;
  }

  cur++;

  looking += factor;
  factor *= 2;
}
END {
  for (i = 0; i < cur; i++) {
    for (j = 1; j <= num; j++) {
      printf("%20.17lf ",x[i,j]);
    }
    printf("\n");
  }
}
