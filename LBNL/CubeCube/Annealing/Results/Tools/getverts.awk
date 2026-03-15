BEGIN {
  max = 0;
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
  size = NF-1;

  for (i = 2; i <= size; i++) {
    f[limit,i] = $i/2.0;
  }

  if ($i > max) max = $i;
  limit++;
}
END {
  for (i = 1; i < limit; i++) {
    for (j = 2; j <= size; j++) {
      printf("%20.17f ", 2*f[i,j]/max);
    }
    printf("\n");
  }
  for (i = limit-1; i >= 1; i--) {
    for (j = 2; j <= size; j++) {
      printf("%20.17f ",-2*f[i,j]/max);
    }
    printf("\n");
  }
}
