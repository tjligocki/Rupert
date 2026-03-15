#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"

inline void printLines1(double *x, double range, int np, int m, int n)
{
  int dim;
  double dx;
  double *y;
  double cur;
  int i,j;
  int ip;

  dim = (m*(2*n-(m+1)))/2;

  dx = range / (np-1);

  y = pointMake(dim);

  for (i = 0; i < dim; i++) {
    y[i] = x[i];
  }
  
  cur = rupert(x,m,n);

  fprintf(stderr,"%.16le %.16le\n",0.0,cur);
  fprintf(stderr,"%.16le %.16le\n",1.0,cur);
  fprintf(stderr,"\n");

  for (i = 0; i < dim; i++) {
    for (j = i+1; j < dim; j++) {
      for (ip = 0; ip < np; ip++) {
        y[i] = x[i] - range/2 + ip*dx;

        y[j] = x[j] - range/2 + ip*dx;

        cur = rupert(y,m,n);

        fprintf(stderr,"%.16le %.16le\n",ip/(double)(np-1),cur);
      }
      fprintf(stderr,"\n");

      for (ip = 0; ip < np; ip++) {
        y[i] = x[i] - range/2 + ip*dx;

        y[j] = x[j] + range/2 - ip*dx;

        cur = rupert(y,m,n);

        fprintf(stderr,"%.16le %.16le\n",ip/(double)(np-1),cur);
      }
      fprintf(stderr,"\n");

      y[j] = x[j];
    }

    y[i] = x[i];
  }

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

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  np = atoi(argv[3]);

  dim = (m*(2*n-(m+1)))/2;

  range = 2*M_PI;

  minx = pointMake(dim);

  i = fread(minx,sizeof(*minx),dim,stdin);

  printLines1(minx,range,np,m,n);

  pointFree(minx);

  exit(0);
}
