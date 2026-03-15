{
  for (i = 1; i <= NF; i++) {
    j = i;
    if (i == c1) {
      j = c2;
    }
    if (i == c2) {
      j = c1;
    }
    sign = substr($j,1,1);
    if (sign == "-") {
      if ($j == 0.0) {
        printf(" %s ",substr($j,2));
      } else {
        printf("%s ",$j);
      }
    } else {
      printf(" %s ",$j);
    }
  }
  printf("\n");
}
