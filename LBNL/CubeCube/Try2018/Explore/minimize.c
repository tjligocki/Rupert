#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"

inline void minimize(double *minx, double range, int m, int n)
{
  int dim;
  double *y;
  double dx;
  int i;
  double cur;

  dx = range / 100.0;

  dim = (m*(2*n-(m+1)))/2;

  y = pointMake(dim);

  cur = rupert(minx,m,n);

  fprintf(stderr,"%18.16lf (%12.6le) - ",cur,dx);
  pointPrint(stderr,"%8.5lf",minx,dim,1);

  while (dx > 1e-11) {
    double lo,hi;

    for (i = 0; i < dim; i++) {
      y[i] = minx[i];
    }

    for (i = 0; i < dim; i++) {
      y[i] = minx[i] - dx/2;
      lo = rupert(y,m,n);

      if (lo < cur) {
        cur = lo;
        minx[i] = y[i];

        break;
      }

      y[i] = minx[i] + dx/2;
      hi = rupert(y,m,n);

      if (hi <= cur) {
        cur = hi;
        minx[i] = y[i];

        break;
      }

      y[i] = minx[i];
    }

    if (i == dim) {
      dx = dx/2;

      fprintf(stderr,"%18.16lf (%12.6le) - ",cur,dx);
      pointPrint(stderr,"%8.5lf",minx,dim,1);
    } else {
      dx = dx*1.0001;

      if (dx > 10*M_PI) {
        break;
      }
    }
  }

  fprintf(stderr,"\n");

  pointFree(y);
}

int main(int argc, char **argv)
{
  int m,n;
  int dim;
  int i;
  double *minx;
  double *minMatrix;
  double minimum;
  double range;
  FILE *output;
  char name[1000];

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  dim = (m*(2*n-(m+1)))/2;

  range = M_PI / 2;

  minx = pointMake(dim);

  i = fread(minx,sizeof(*minx),dim,stdin);

  if (i < dim) {
    for (i = 0; i < dim; i++) {
      minx[i] = 1.0/(i+1);
    }
  }

  minimize(minx,range,m,n);

  pointPrint(stderr,"%8.5lf",minx,dim,0);
  fprintf(stderr," - %20.15lf\n",1.0/rupert(minx,m,n));
  fprintf(stderr,"\n");

  fwrite(minx,sizeof(*minx),dim,stdout);

  pointFree(minx);

  exit(0);
}
