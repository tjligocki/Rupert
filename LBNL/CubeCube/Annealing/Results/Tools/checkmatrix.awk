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
    for (j = 0; j < cur; j++) {
      dot = 0.0
      for (k = 1; k <= num; k++) {
        dot += x[i,k]*x[j,k];
      }
      printf("%19.16lf ",dot);
    }
    printf("\n");
  }
}
