{
  m = 2*NR+1;
  n = 2*NR+3;
  printf("%3d %3d    %s\n",m,n,$1);
  for (i = 2; i*n <= 100; i++) {
    printf("%3d %3d    %s\n",i*m,i*n,$1);
  }
}
