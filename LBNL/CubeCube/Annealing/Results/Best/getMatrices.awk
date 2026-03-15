NF == n+2 {
  for (i = 1; i <= n; i++) {
    if (substr($(i+1),1,1) == "-") {
      printf("%s ",$(i+1));
    } else {
      printf(" %s ",$(i+1));
    }
  }
  printf("\n");
}
