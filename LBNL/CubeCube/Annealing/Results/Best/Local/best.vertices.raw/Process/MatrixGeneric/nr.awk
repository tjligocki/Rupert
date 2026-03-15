{
  for (i = 1; i <= NF; i++) {
    sign = substr($i,1,1);
    if (sign == "-") {
      printf(" %s ",substr($i,2));
    } else {
      if ($i == 0.0) {
        printf(" %s ",$i);
      } else {
        printf("-%s ",$i);
      }
    }
  }
  printf("\n");
}
