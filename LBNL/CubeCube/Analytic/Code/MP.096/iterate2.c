#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "floatType.h"
#include "matrix.h"

#define DEBUG   1
#define DEBUG2  1

main(int argc, char **argv)
{
  int n;
  REAL** mat;
  char format[80];

  n = 7;

  if (argc > 1) {
    n = atoi(argv[1]);
  }

  getPrimeMatrix(&mat,n);

  printMatrix(mat,n);
  fprintf(stderr,"\n");

  sprintf(format,"Side: %%19.16%s\n\n",PRINTF);
  fprintf(stderr,format,(REAL)1.0/getMaxSum2(mat,n,n-1));

  exit(0);
}
