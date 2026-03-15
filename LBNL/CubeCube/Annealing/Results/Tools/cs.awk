{
  n = split(column,c);
  comp = NF;
  for (i = 1; i <= comp; i++) {
    for (j = 1; j <= n; j++) {
      if (i == c[j]) {
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
