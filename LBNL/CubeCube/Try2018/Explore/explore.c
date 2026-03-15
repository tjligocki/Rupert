#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"

inline void printDirs1(double *x, double dx, int m, int n)
{
  int dim;
  double *y;
  int i;
  double cur;

  dim = (m*(2*n-(m+1)))/2;

  y = pointMake(dim);

  cur = rupert(x,m,n);

  fprintf(stderr,"%18.16lf (%12.6le) - ",cur,dx);
  pointPrint(stderr,"%8.5lf",x,dim,1);
  fprintf(stderr,"\n");

  for (i = 0; i < dim; i++) {
    y[i] = x[i];
  }

  for (i = 0; i < dim; i++) {
    double lo,hi;

    fprintf(stderr,"  %2d: ",i);

    y[i] = x[i] - dx/2;
    lo = rupert(y,m,n);

    fprintf(stderr,"lo %12.6le",lo-cur);
    if (lo-cur <= 0.0) {
      fprintf(stderr," ---");
    }
    fprintf(stderr,"\n");

    y[i] = x[i] + dx/2;
    hi = rupert(y,m,n);

    fprintf(stderr,"      ");
    fprintf(stderr,"hi %12.6le",hi-cur);
    if (hi-cur <= 0.0) {
      fprintf(stderr," ---");
    }
    fprintf(stderr,"\n");

    y[i] = x[i];
  }

  fprintf(stderr,"\n");

  pointFree(y);
}

inline void printDirs2(double *x, double dx, int m, int n)
{
  int dim;
  double *y;
  int i,j;
  double cur;

  dim = (m*(2*n-(m+1)))/2;

  y = pointMake(dim);

  cur = rupert(x,m,n);

  fprintf(stderr,"%18.16lf (%12.6le) - ",cur,dx);
  pointPrint(stderr,"%8.5lf",x,dim,1);
  fprintf(stderr,"\n");

  for (i = 0; i < dim; i++) {
    y[i] = x[i];
  }

  for (i = 0; i < dim; i++) {
    for (j = i+1; j < dim; j++) {
      double val;

      y[i] = x[i] - dx/2;

      y[j] = x[j] - dx/2;
      val = rupert(y,m,n);

      fprintf(stderr,"  %2d %2d: ",i,j);
      fprintf(stderr,"lo lo %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

      y[j] = x[j] + dx/2;
      val = rupert(y,m,n);

      fprintf(stderr,"         ");
      fprintf(stderr,"lo hi %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

      y[i] = x[i] + dx/2;

      y[j] = x[j] - dx/2;
      val = rupert(y,m,n);

      fprintf(stderr,"         ");
      fprintf(stderr,"hi lo %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

      y[j] = x[j] + dx/2;
      val = rupert(y,m,n);

      fprintf(stderr,"         ");
      fprintf(stderr,"hi hi %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

      fprintf(stderr,"\n");

      y[j] = x[j];
    }

    y[i] = x[i];
  }

  fprintf(stderr,"\n");

  pointFree(y);
}

inline void printDirs3(double *x, double dx, int m, int n)
{
  int dim;
  double *y;
  int i,j,k;
  double cur;

  dim = (m*(2*n-(m+1)))/2;

  y = pointMake(dim);

  cur = rupert(x,m,n);

  fprintf(stderr,"%18.16lf (%12.6le) - ",cur,dx);
  pointPrint(stderr,"%8.5lf",x,dim,1);
  fprintf(stderr,"\n");

  for (i = 0; i < dim; i++) {
    y[i] = x[i];
  }

  for (i = 0; i < dim; i++) {
    for (j = i+1; j < dim; j++) {
      for (k = j+1; k < dim; k++) {
        double val;

        y[i] = x[i] - dx/2;

        y[j] = x[j] - dx/2;

        y[k] = x[k] - dx/2;
        val = rupert(y,m,n);

        fprintf(stderr,"  %2d %2d %2d: ",i,j,k);
        fprintf(stderr,"lo lo lo %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

        y[k] = x[k] + dx/2;
        val = rupert(y,m,n);

        fprintf(stderr,"            ");
        fprintf(stderr,"lo lo hi %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

        y[j] = x[j] + dx/2;

        y[k] = x[k] - dx/2;
        val = rupert(y,m,n);

        fprintf(stderr,"            ");
        fprintf(stderr,"lo hi lo %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

        y[k] = x[k] + dx/2;
        val = rupert(y,m,n);

        fprintf(stderr,"            ");
        fprintf(stderr,"lo hi hi %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

        y[i] = x[i] + dx/2;

        y[j] = x[j] - dx/2;

        y[k] = x[k] - dx/2;
        val = rupert(y,m,n);

        fprintf(stderr,"            ");
        fprintf(stderr,"hi lo lo %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

        y[k] = x[k] + dx/2;
        val = rupert(y,m,n);

        fprintf(stderr,"            ");
        fprintf(stderr,"hi lo hi %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

        y[j] = x[j] + dx/2;

        y[k] = x[k] - dx/2;
        val = rupert(y,m,n);

        fprintf(stderr,"            ");
        fprintf(stderr,"hi hi lo %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

        y[k] = x[k] + dx/2;
        val = rupert(y,m,n);

        fprintf(stderr,"            ");
        fprintf(stderr,"hi hi hi %12.6le",val-cur);
        if (val-cur <= 0.0) {
          fprintf(stderr," ---");
        }
        fprintf(stderr,"\n");

        fprintf(stderr,"\n");

        y[k] = x[k];
      }

      y[j] = x[j];
    }

    y[i] = x[i];
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
  FILE *input;
  char name[1000];
  double *minMatrix;
  double *minSums;
  double minimum;
  double range;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  dim = (m*(2*n-(m+1)))/2;

  range = M_PI / 2;

  minx = pointMake(dim);

  i = fread(minx,sizeof(*minx),dim,stdin);

  pointPrint(stderr,"%8.5lf",minx,dim,0);
  fprintf(stderr," - %20.15lf\n",1.0/rupert(minx,m,n));
  fprintf(stderr,"\n");

  minMatrix = rupertMatrix(minx,m,n);
  matrixPrint(stderr,"%8.5lf",minMatrix,n);
  fprintf(stderr,"\n");
  matrixFree(minMatrix);

  minSums = rupertSums(minx,m,n);
  pointPrint(stderr,"%8.5lf",minSums,n,1);
  fprintf(stderr,"\n");
  pointFree(minSums);

  printDirs1(minx,range/1000,m,n);

  printDirs2(minx,range/1000,m,n);

  printDirs3(minx,range/1000,m,n);

  pointFree(minx);

  exit(0);
}
