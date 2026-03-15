{
  n = split(row,r);
  comp = NF;
  for (i = 1; i <= comp; i++) {
    for (j = 1; j <= n; j++) {
      if (NR == r[j]) {
        printf("%10.7f ",-$i);
        break;
      }
    }

    if (j > n) {
      printf("%10.7f ",$i);
    }
  }
  printf("\n");
}
