{
  comp = NF;
  for (i = 1; i <= comp; i++) {
    if (i == column1) {
      printf("%10.7f ",$column2);
    } else
    if (i == column2) {
      printf("%10.7f ",$column1);
    } else {
      printf("%10.7f ",$i);
    }
  }
  printf("\n");
}
