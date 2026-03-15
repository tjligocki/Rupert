{
  if (NR == 1) {
    m = step+1;
    n = 2*step+1;
  }

  printf("%3d %3d    %s\n",m,n,$1);

  for (i = 2; i*n <= 100; i++) {
    printf("%3d %3d    %s\n",i*m,i*n,$1);
  }

  m += step;
  n += step;

  if (n > max) {
    step++;
    
    m = step+1;
    n = 2*step+1;
  }
}
