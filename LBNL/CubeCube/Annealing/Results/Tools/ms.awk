BEGIN {
  cur = 0;
  num = 0;
}
{
  for (i = 1; i <= NF; i++) {
    x[cur,i] = $i;
  }

  cur++;
  num = NF;
}
END {
  for (k = 1; k <= num; k++) {
    r1 = x[row1-1,k];
    s1 = 0;

    if (r1 < 0) {
      r1 = -r1;
      s1 = 1;
    }

    r2 = x[row2-1,k];
    s2 = 0;

    if (r2 < 0) {
      r2 = -r2;
      s2 = 1;
    }

    r = (r1+r2)/2.0;

    if (s1 == 0) {
      x[row1-1,k] = r;
    } else {
      x[row1-1,k] = -r;
    }

    if (s2 == 0) {
      x[row2-1,k] = r;
    } else {
      x[row2-1,k] = -r;
    }
  }
  for (i = 0; i < cur; i++) {
    for (k = 1; k <= num; k++) {
      printf("%19.16lf ",x[i,k]);
    }
    printf("\n");
  }
}
