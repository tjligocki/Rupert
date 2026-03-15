#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "tools.h"

inline void printLines1(double *x, double range, int nr, int np, int m, int n)
{
  int dim;
  double dx;
  double *vx;
  double *y;
  double cur;
  int ir;
  int i;
  int ip;

  dim = (m*(2*n-(m+1)))/2;

  dx = range / (np-1);

  vx = pointMake(dim);
  y  = pointMake(dim);

  cur = rupert(x,m,n);

  for (ip = 0; ip < np; ip++) {
    fprintf(stderr,"%.16le %.16le\n",ip/(double)(np-1),cur);
  }
  fprintf(stderr,"\n");

  for (ir = 0; ir < nr; ir++) {
    double nvx;

    nvx = 0.0;

    for (i = 0; i < dim; i++) {
      double cv;

      cv = 2*drand48() - 1;

      vx[i] = cv;
      nvx += cv*cv;
    }

    nvx = sqrt(nvx);

    for (i = 0; i < dim; i++) {
      vx[i] = vx[i] / nvx;
      y[i] = x[i];
    }
    
    for (ip = 0; ip < np; ip++) {
      for (i = 0; i < dim; i++) {
        y[i] = x[i] + (-range/2 + ip*dx)*vx[i];
      }

      cur = rupert(y,m,n);

      fprintf(stderr,"%.16le %.16le\n",ip/(double)(np-1),cur);
    }
    fprintf(stderr,"\n");
  }

  pointFree(vx);
  pointFree(y);
}

int main(int argc, char **argv)
{
  int m,n;
  int np;
  int nr;
  int dim;
  int i;
  double *minx;
  FILE *input;
  char name[1000];
  double *minMatrix;
  double minimum;
  double range;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  np = atoi(argv[3]);
  nr = atoi(argv[4]);

  srand48(time(0));

  dim = (m*(2*n-(m+1)))/2;

  range = 2*M_PI;

  minx = pointMake(dim);

  i = fread(minx,sizeof(*minx),dim,stdin);

  printLines1(minx,range/4,nr,np,m,n);

  pointFree(minx);

  exit(0);
}
