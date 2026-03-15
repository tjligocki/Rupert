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
  for (i = 0; i < cur; i++) {
    if (i != row-1) {
      for (k = 1; k <= num; k++) {
        printf("%19.16lf ",x[i,k]);
      }
      printf("\n");
    } else {
      for (j = 0; j < cur; j++) {
        if (j != i) {
          dot = 0.0
          for (k = 1; k <= num; k++) {
            dot += x[i,k]*x[j,k];
          }
          for (k = 1; k <= num; k++) {
            x[i,k] = x[i,k] - dot*x[j,k];
          }
        }
      }
      dot = 0.0
      for (k = 1; k <= num; k++) {
        dot += x[i,k]*x[i,k];
      }
      dot = sqrt(dot);
      for (k = 1; k <= num; k++) {
        x[i,k] = x[i,k]/dot;
      }
      for (k = 1; k <= num; k++) {
        printf("%19.16lf ",x[i,k]);
      }
      printf("\n");
    }
  }
}
