#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "tools.h"

inline printLines(double *x, double range, int np, int nr, int nt, int m, int n)
{
  int dim;
  double dx;
  double *vx;
  double *y;
  double cur,min,imin;
  int i;
  int it;
  double nvx;

  dim = (m*(2*n-(m+1)))/2;

  dx = range / (np-1);

  vx = pointMake(dim);
  y  = pointMake(dim);

  cur = rupert(x,m,n);

  min = cur;
  imin = 0;
  
  for (i = 0; i < dim; i++) {
    y[i] = x[i];
  }

  fprintf(stderr,"%10d %23.16le\n",0,cur);

  do {
    nvx = 0.0;

    for (i = 0; i < dim; i++) {
      double cv;

      cv = 2*drand48() - 1;

      vx[i] = cv;
      nvx += cv*cv;
    }
  } while (nvx > 1.0);

  nvx = sqrt(nvx);

  for (i = 0; i < dim; i++) {
    vx[i] = vx[i] / nvx;
  }

  for (it = 0; it <= nt; it++) {
    for (i = 0; i < dim; i++) {
      y[i] = y[i] + dx*vx[i];
    }

    cur = rupert(y,m,n);

    if (cur < min) {
      min = cur;
      imin = it+1;

      fprintf(stderr,"%10d %23.16le\n",it,min);
    } else {
      if (it % nr == 0) {
        fprintf(stderr,"%10d %23.16le %23.16le\n",it,min,cur);
      }
    }
  }

  for (i = 0; i < dim; i++) {
    x[i] = x[i] + imin*dx*vx[i];
  }

  pointFree(vx);
  pointFree(y);
}

int main(int argc, char **argv)
{
  int m,n;
  int np;
  int dim;
  int i;
  double *minx;
  FILE *input;
  char name[1000];
  double *minMatrix;
  double minimum;
  double range;
  int nr,nt;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  np = atoi(argv[3]);
  nr = atoi(argv[4]);
  nt = atoi(argv[5]);

  srand48(time(0));

  dim = (m*(2*n-(m+1)))/2;

  range = 2*M_PI;

  minx = pointMake(dim);

  i = fread(minx,sizeof(*minx),dim,stdin);

  printLines(minx,range,np,nr,nt,m,n);

  fwrite(minx,sizeof(*minx),dim,stdout);

  pointFree(minx);

  exit(0);
}
