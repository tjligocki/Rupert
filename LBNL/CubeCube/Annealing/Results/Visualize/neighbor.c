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

#define RANGE   1

main(int argc, char** argv)
{
  int indim,outdim;
  int size;
  int i,n;
  int* is;
  int axis1,axis2;
  int iter;
  double** best;
  double** trans;
  double*** new;
  double tstart,tend,dt;
  double t;
  double range;
  double value;
  int found;
  int* is0;
  double v0,vs;
  double** save;

  indim = 2;
  outdim = 3;

  tstart = 0.0;
  tend = 2*M_PI;

  n = 10;

  if (argc > 1) {
    indim = atoi(argv[1]);
  }

  if (argc > 2) {
    outdim = atoi(argv[2]);
  }

  if (argc > 3) {
    n = atoi(argv[3]);
  }

  size = (outdim * (outdim-1)) / 2;

  dt = (tend - tstart) / (n - 1);

  allocMatrix(&best,outdim);
  identMatrix(best,outdim);

  allocMatrix(&trans,outdim);

  new = (double ***)malloc((size+1)*sizeof(double **));

  for (i = 0; i <= size; i++) {
    allocMatrix(&new[i],outdim);
  }

  allocMatrix(&save,outdim);

  is  = (int *)malloc(size*sizeof(int));
  is0 = (int *)malloc(size*sizeof(int));

  v0 = 1.0;

  range = RANGE;

  for (i = 0; i < size; i++) {
    is0[i] = 0;
  }

  iter = 0;

  do {
    found = 0;

    for (i = 0; i < size; i++) {
      is[i] = -range;
    }

    vs = v0;

    while (1) {
      int atzero;

      atzero = 1;
      for (i = 0; i < size; i++) {
        if (is[i] != 0) {
          atzero = 0;
          break;
        }
      }

      copyMatrix(new[0],best,outdim);

      i = 0;
      for (axis1 = 0; axis1 < outdim-1; axis1++) {
        for (axis2 = axis1+1; axis2 < outdim; axis2++) {
          t = is[i]*dt + tstart;

          identMatrix(trans,outdim);

          trans[axis1][axis1] =  cos(t);
          trans[axis1][axis2] =  sin(t);

          trans[axis2][axis1] = -sin(t);
          trans[axis2][axis2] =  cos(t);

          multiplyMatrix(new[i+1],new[i],trans,outdim);

          i++;
        }
      }

      if (atzero != 1) {
        value = normRupertMatrix(new[size],indim,outdim);

        /* fprintf(stderr,"  %20.13le\n",1.0/value); */

        if (value <= v0) {
          if (value == vs) {
            found = 2;
          } else {
            found = 1;
          }

          for (i = 0; i < size; i++) {
            is0[i] = is[i];
          }

          v0 = value;
          copyMatrix(save,new[size],outdim);
        }
      }

      for (i = size-1; i >= 0; i--) {
        is[i]++;

        if (is[i] <= range) {
          break;
        } else {
          is[i] = -range;
        }
      }

      if (i == -1) {
        break;
      }
    }

    iter++;

    // fprintf(stderr,"\n");
    // printMatrix(best,outdim);
    // fprintf(stderr,"\n");

    for (i = 0; i < size; i++) {
      fprintf(stderr,"%3d ",is0[i]);
    }
    fprintf(stderr,"%30.23le (%20.13le %5.1lf %d) %d\r",1.0/v0,dt,range,found,iter);

    if (found == 0) {
      dt *= 0.8180339887498948;
      if (pow(2*(range*1.63)+1,size) < 30000.0) {
        range *= 1.63;
      }
    } else {
      if (found == 1) {
        dt *= 1.32;
        if (range/1.62 >= 1.0) {
          range /= 1.62;
        }
      } else {
        dt *= 0.8180339887498948;
        if (pow(2*(range*1.63)+1,size) < 30000.0) {
          range *= 1.63;
        }
      }

      orthoNormMatrix(save,outdim);
      copyMatrix(best,save,outdim);
    }
  } while (dt > 1e-20);

  fprintf(stderr,"\n");
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
