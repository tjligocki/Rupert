{
  n = split(column,c);
  comp = NF;
  for (i = 1; i <= comp; i++) {
    for (j = 1; j <= n; j++) {
      if (i == c[j]) {
        if (substr($i,1,1) == "-") {
          printf(" %s ",substr($i,2));
        } else {
          printf("-%s ",$i);
        }
        break;
      }
    }

    if (j > n) {
      if (substr($i,1,1) == "-") {
        printf("%s ",$i);
      } else {
        printf(" %s ",$i);
      }
    }
  }
  printf("\n");
}
