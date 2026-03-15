{
  max = 0.0;
  for (c = 1; c <= m; c++) {
    mat[NR,c] = $c
    if ($c > 0) {
      max += $c;
    } else {
      max += -$c;
    }
  }
}
END {
  for (d = 1; d < m; d++) {
    i[d] = -1;
  }
  i[m] = 1;

  while (1) {
    for (r = 1; r <= n; r++) {
      sum = 0.0;
      for (d = 1; d <= m; d++) {
        sum += i[d]*mat[r,d]
      }
      printf("%10.7f ",sum/max);
    }
    printf("\n");

    i[1] += 2; 

    cur = 1;
    while (i[cur] > 1) {
      i[cur] = -1;
      cur++;

      if (cur > m) {
        break;
      }

      i[cur] += 2;
    }

    if (cur > m) {
      break;
    }
  }
}
