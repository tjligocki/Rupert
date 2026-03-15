#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"

inline double *getGrads(double *x, double dx, int m, int n)
{
  int dim;
  double *y;
  int i;
  double *grads;

  dim = (m*(2*n-(m+1)))/2;

  y = pointMake(dim);

  grads = pointMake(n*dim);

  for (i = 0; i < dim; i++) {
    y[i] = x[i];
  }

  for (i = 0; i < dim; i++) {
    double lo,hi;
    double *rupMat;
    double *sumsLo;
    double *sumsHi;
    int j;

    y[i] = x[i] - dx/2;
    rupMat = rupertMatrix(y,m,n);

    sumsLo = matrixSums(rupMat,m,n);

    matrixFree(rupMat);

    y[i] = x[i] + dx/2;
    rupMat = rupertMatrix(y,m,n);

    sumsHi = matrixSums(rupMat,m,n);

    for (j = 0; j < n; j++) {
      int index;

      index = j*dim + i;

      grads[index] = (sumsHi[j] - sumsLo[j]) / dx;
    }

    matrixFree(rupMat);

    pointFree(sumsLo);
    pointFree(sumsHi);

    y[i] = x[i];
  }

  pointFree(y);

  return(grads);
}

int main(int argc, char **argv)
{
  int m,n;
  int dim;
  int i;
  double *minx;
  FILE *input;
  char name[1000];
  double *minMatrix;
  double minimum;
  double range;
  double *grads;
  double *sumGrads;
  double *x;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  dim = (m*(2*n-(m+1)))/2;

  range = M_PI / 2;

  minx = pointMake(dim);

  i = fread(minx,sizeof(*minx),dim,stdin);

  pointPrint(stderr,"%8.5lf",minx,dim,0);
  fprintf(stderr," - %20.15lf\n",1.0/rupert(minx,m,n));
  fprintf(stderr,"\n");

  grads = getGrads(minx,range/1000,m,n);

  sumGrads = pointMake(dim);
  for (i = 0; i < dim; i++) {
    sumGrads[i] = 0.0;
  }

  for (i = 0; i < n; i++) {
    int j;

    fprintf(stderr,"%2d - ",i);
    pointPrint(stderr,"%8.5lf",grads+i*dim,dim,1);

    for (j = 0; j < dim; j++) {
      sumGrads[j] += grads[i*dim+j];
    }
  }
  fprintf(stderr,"\n");

  pointPrint(stderr,"%8.5lf",sumGrads,dim,1);
  fprintf(stderr,"\n");

  x = pointMake(dim);
  fprintf(stderr,"0: %20.15lf\n",rupert(minx,m,n));
/*
  for (i = -1000; i <= 1000; i++) {
    int j;
    double t;

    t = i/1000.0;

    for (j = 0; j < dim; j++) {
      x[j] = minx[j] + t * sumGrads[j];
    }

    fprintf(stderr,"%20.15lf %20.15lf\n",t,rupert(x,m,n));
  }
  fprintf(stderr,"\n");
*/
  i = 2;/* for (i = n-1; i < n; i++) */ {
    int j;

    for (j = -1000; j <= 1000; j++)
    {
      int k;
      double t;
      double *sums;

      t = j/1000.0;

      for (k = 0; k < dim; k++) {
        x[k] = minx[k] + t * grads[i*dim+k];
      }

      sums = rupertSums(x,m,n);

      fprintf(stderr,"%20.15lf ",t);
      for (k = 0; k < n; k++) {
        fprintf(stderr,"%20.15lf ",sums[k]);
      }
      fprintf(stderr,"\n");

      pointFree(sums);
    }
    fprintf(stderr,"\n");
  }

  pointFree(minx);
  pointFree(grads);
  pointFree(sumGrads);
  pointFree(x);

  exit(0);
}
