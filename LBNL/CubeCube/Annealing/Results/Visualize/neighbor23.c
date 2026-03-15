#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void allocMatrix(double*** pm, int dim);
void identMatrix(double** pm, int dim);

void printMatrix(double** m, int dim);

void copyMatrix(double** mDst, double** mSrc, int dim);
void multiplyMatrix(double** mDst, double** mSrc1, double** mSrc2, int dim);

void orthoNormMatrix(double** m, int dim);
double normRupertMatrix(double** m, int indim, int dim);

#define RANGE   3

main(int argc, char** argv)
{
  int outdim;
  int i,j,k,n;
  double** best;
  double** trans;
  double** new1;
  double** new2;
  double** new3;
  double tstart,tend,dt;
  double t1,t2,t3;
  double range;
  double value;
  int found;
  int i0,j0,k0;
  double v0;
  double** save;

  outdim = 3;

  tstart = 0.0;
  tend = 2*M_PI;

  n = 10000;

  if (argc > 1) {
    n = atoi(argv[1]);
  }

  dt = (tend - tstart) / (n - 1);

  allocMatrix(&best,outdim);
  identMatrix(best,outdim);

  allocMatrix(&trans,outdim);

  allocMatrix(&new1,outdim);
  allocMatrix(&new2,outdim);
  allocMatrix(&new3,outdim);

  allocMatrix(&save,outdim);

  v0 = 1.0;

  range = RANGE;

  do {
    found = 0;

    for (i = -range; i <= range; i++) {
      t1 = i*dt + tstart;

      identMatrix(trans,outdim);

      trans[0][0] =  cos(t1);
      trans[0][1] =  sin(t1);

      trans[1][0] = -sin(t1);
      trans[1][1] =  cos(t1);

      multiplyMatrix(new1,best,trans,outdim);

      for (j = -range; j <= range; j++) {
        t2 = j*dt + tstart;

        identMatrix(trans,outdim);

        trans[0][0] =  cos(t2);
        trans[0][2] =  sin(t2);

        trans[2][0] = -sin(t2);
        trans[2][2] =  cos(t2);

        multiplyMatrix(new2,new1,trans,outdim);

        for (k = -range; k <= range; k++) {
          if (i != 0 || j != 0 || k != 0) {

            t3 = k*dt + tstart;

            identMatrix(trans,outdim);

            trans[1][1] =  cos(t3);
            trans[1][2] =  sin(t3);

            trans[2][1] = -sin(t3);
            trans[2][2] =  cos(t3);

            multiplyMatrix(new3,new2,trans,outdim);

            value = normRupertMatrix(new3,2,outdim);

            // fprintf(stderr,"  %2d %2d %2d %20.13le\n",i,j,k,1.0/value);
            if (v0 > value) {
              found = 1;

              i0 = i;
              j0 = j;
              k0 = k;

              v0 = value;
              copyMatrix(save,new3,outdim);
            }
          }
        }
      }
    }

    // fprintf(stderr,"\n");
    fprintf(stderr,"%4d %4d %4d %30.23le (%20.13le %5.1lf)\n",
            i0,j0,k0,1.0/v0,dt,range);
    // printMatrix(save,outdim);
    // fprintf(stderr,"\n");

    if (found == 0) {
      dt *= 0.65;
      if (range < 20.0) {
        range *= 1.6;
      }
    } else {
      dt *= 2.0;
      if (range/2.0 >= 1.0) {
        range /= 2.0;
      }
    }

    orthoNormMatrix(save,outdim);

    copyMatrix(best,save,outdim);
  } while (dt > 1e-20);
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

void orthoNormMatrix(double** m, int dim)
{
  int i,j,k;

  for (i = dim-1; i >= 0; i--) {
    double norm;

    norm = 0.0;
    for (j = 0; j < dim; j++) {
      norm += m[i][j]*m[i][j];
    }
    norm = sqrt(1.0/norm);

    for (j = 0; j < dim; j++) {
      m[i][j] *= norm;
    }

    for (j = i-1; j >= 0; j--) {
      double dot;

      dot = 0.0;
      for (k = 0; k < dim; k++) {
        dot += m[i][k]*m[j][k];
      }

      for (k = 0; k < dim; k++) {
        m[j][k] -= dot*m[i][k];
      }
    }
  }
}
