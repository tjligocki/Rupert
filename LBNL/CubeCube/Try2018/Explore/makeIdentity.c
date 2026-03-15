#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"

int main(int argc, char **argv)
{
  int m,n;
  int dim;
  int i;
  double *x;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  dim = (m*(2*n-(m+1)))/2;

  x = pointMake(dim);

  for (i = 0; i < dim; i++) {
    x[i] = 0.0;
  }

  fwrite(x,sizeof(*x),dim,stdout);

  pointFree(x);

  exit(0);
}
