BEGIN {
  size = 0;
  start = 0;
  stop = 0;
  limit = 1;
}
NF != 0 && start == 0 {
  start = 1;
}
start == 1 && NF == 0 {
  stop = 1;
}
start == 1 && stop == 0 {
  size = NF - 1;

  for (i = 2; i <= size; i++) {
    f[limit,i] = $i;
  }

  limit++;
}
END {
  for (i = 1; i < limit; i++) {
    for (j = 2; j <= size; j++) {
      printf("%s ", f[i,j]);
    }
    printf("\n");
  }
}
