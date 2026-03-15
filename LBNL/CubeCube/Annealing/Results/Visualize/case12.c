#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void allocMatrix(double*** pm, int dim);
void identMatrix(double** pm, int dim);

void printMatrix(double** m, int dim);

void copyMatrix(double** mDst, double** mSrc, int dim);
void multiplyMatrix(double** mDst, double** mSrc1, double** mSrc2, int dim);

double normRupertMatrix(double** m, int indim, int dim);

main(int argc, char** argv)
{
  int i,n;
  double** ident;
  double** trans;
  double** new;
  double tstart,tend,dt;
  double t;

  tstart = 0.0;
  tend = 2*M_PI;

  n = 30;

  if (argc > 1) {
    n = atoi(argv[1]);
  }

  dt = (tend - tstart) / (n - 1);

  allocMatrix(&ident,2);
  identMatrix(ident,2);

  allocMatrix(&trans,2);

  allocMatrix(&new,2);

  for (i = 0; i < n; i++) {
    double value;

    t = i*dt + tstart;

    trans[0][0] =  cos(t);
    trans[0][1] =  sin(t);

    trans[1][0] = -sin(t);
    trans[1][1] =  cos(t);

    multiplyMatrix(new,ident,trans,2);

    value = normRupertMatrix(new,1,2);

    printf("%20.13le %20.13le\n",t,1.0/value);
  }
}

void allocMatrix(double*** pm, int dim)
{
  double** m;
  int i,j;

  m = (double**)malloc(dim*sizeof(double *));

  for (i = 0; i < dim; i++) {
    m[i] = (double *)malloc(dim*sizeof(double));
  }

  *pm = m;
}

void identMatrix(double** m, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if (i == j) {
        m[i][j] = 1.0;
      } else {
        m[i][j] = 0.0;
      }
    }
  }
}

void printMatrix(double** m, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      fprintf(stderr,"%19.16f ",m[i][j]);
    }
    fprintf(stderr,"\n");
  }
}

void copyMatrix(double** mDst, double** mSrc, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mDst[i][j] = mSrc[i][j];
    }
  }
}

void multiplyMatrix(double** mDst, double** mSrc1, double** mSrc2, int dim)
{
  int i,j,k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      double sum;

      sum = 0.0;
      for (k = 0; k < dim; k++) {
        sum += mSrc1[i][k]*mSrc2[k][j];
      }

      mDst[i][j] = sum;
    }
  }
}

double normRupertMatrix(double** m, int indim, int dim)
{
  int i,j;
  double max,cur;

  max = 0.0;

  for (i = 0; i < dim; i++) {

    cur = 0.0;
    for (j = 0; j < indim; j++) {
      cur += fabs(m[i][j]);
    }

    if (cur > max) {
      max = cur;
    }
  }

  return max;
}
