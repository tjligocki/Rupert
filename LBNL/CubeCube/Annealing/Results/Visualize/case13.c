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
  int i,j,k,n;
  double** ident;
  double** trans;
  double** new1;
  double** new2;
  double** new3;
  double tstart,tend,dt;
  double t1,t2,t3;

  tstart = 0.0;
  tend = 2*M_PI;

  n = 10;

  if (argc > 1) {
    n = atoi(argv[1]);
  }

  dt = (tend - tstart) / (n - 1);

  allocMatrix(&ident,3);
  identMatrix(ident,3);

  allocMatrix(&trans,3);

  allocMatrix(&new1,3);
  allocMatrix(&new2,3);
  allocMatrix(&new3,3);

  for (i = 0; i < n; i++) {
    t1 = i*dt + tstart;

    identMatrix(trans,3);

    trans[0][0] =  cos(t1);
    trans[0][1] =  sin(t1);

    trans[1][0] = -sin(t1);
    trans[1][1] =  cos(t1);

    multiplyMatrix(new1,ident,trans,3);

    for (j = 0; j < n; j++) {
      t2 = j*dt + tstart;

      identMatrix(trans,3);

      trans[0][0] =  cos(t2);
      trans[0][2] =  sin(t2);

      trans[2][0] = -sin(t2);
      trans[2][2] =  cos(t2);

      multiplyMatrix(new2,new1,trans,3);

      for (k = 0; k < n; k++) {
        double value;

        t3 = k*dt + tstart;

        identMatrix(trans,3);

        trans[1][1] =  cos(t3);
        trans[1][2] =  sin(t3);

        trans[2][1] = -sin(t3);
        trans[2][2] =  cos(t3);

        multiplyMatrix(new3,new2,trans,3);

        value = normRupertMatrix(new3,1,3);

        printf("%20.13le %20.13le %20.13le %20.13le\n",t1,t2,t3,1.0/value);
      }
   }
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
