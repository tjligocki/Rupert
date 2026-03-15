{
  comp = NF;
  for (i = 1; i <= comp; i++) {
    if (i == column1) {
      if (substr($column2,1,1) == "-") {
        printf("%s ",$column2);
      } else {
        printf(" %s ",$column2);
      }
    } else
    if (i == column2) {
      if (substr($column1,1,1) == "-") {
        printf("%s ",$column1);
      } else {
        printf(" %s ",$column1);
      }
    } else {
      if (substr($i,1,1) == "-") {
        printf("%s ",$i);
      } else {
        printf(" %s ",$i);
      }
    }
  }
  printf("\n");
}
