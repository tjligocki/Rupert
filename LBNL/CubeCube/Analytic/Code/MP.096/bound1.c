#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "floatType.h"
#include "matrix.h"

main(int argc, char **argv)
{
  int m,n;
  int num;
  int div;
  int m1,m2,mm;
  REAL** mat;
  char format[80];

  m = 4;
  n = 7;

  if (argc > 2) {
    m = atoi(argv[1]);
    n = atoi(argv[2]);
  }

  if (2*m < n) {
    fprintf(stderr,"m, %d, must be >= n/2, %.1f\n",m,n/2.0);
    exit(1);
  }

  num = n - m;

  for (div = 1; div <= num; div++) {
    if ((div - num + m) % num == 0) {
      break;
    }
  }

  m1 = (div - num + m) / num + 1;
  m2 = m1 + 1;

  if (div < num) {
    mm = m2;
  } else {
    mm = m1;
  }

  getPrimeMatrix(&mat,mm);
  sprintf(format,"%%3d %%3d     %%19.17%s\n",PRINTF);
  fprintf(stderr,format,m,n,(REAL)1.0/getMaxSum2(mat,mm,mm-1));

  exit(0);
}
