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

#define EPSILON 1e-16

main(int argc, char** argv)
{
  int indim,outdim;
  int size;
  int i,j,k,n;
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
  int dir1,dir2;

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

    vs = v0;

    /* Do coordinate directions */
    for (dir1 = 0; dir1 < size; dir1++) {
      for (j = 0; j < size; j++) {
        is[j] = 0;
      }

      for (j = -range; j <= range; j++) {
        if (j != 0) {

          is[dir1] = j;

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

          value = normRupertMatrix(new[size],indim,outdim);

          if (value <= v0) {
            if ((vs - value) < EPSILON) {
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
      }
    }

    /* Do coordinate direction pairs*/
    for (dir1 = 0; dir1 < size; dir1++) {
      for (dir2 = dir1+1; dir2 < size; dir2++) {
        for (j = 0; j < size; j++) {
          is[j] = 0;
        }

        for (j = -range; j <= range; j++) {
          for (k = -range; k <= range; k++) {
            if (j != 0 && k != 0) {
              is[dir1] = j;
              is[dir2] = k;

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

              value = normRupertMatrix(new[size],indim,outdim);

              if (value <= v0) {
                if (vs - value < EPSILON) {
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
          }
        }
      }
    }

    iter++;

    /* fprintf(stderr,"\n"); */
    /* printMatrix(best,outdim); */
    /* fprintf(stderr,"\n"); */

    /*
    for (i = 0; i < size; i++) {
      fprintf(stderr,"%3d ",is0[i]);
    }
    */
    fprintf(stderr,"%30.23le (%20.13le %5.1lf %d) %d\r",1.0/v0,dt,range,found,iter);

    if (found == 0) {
      double newrange;

      dt *= 0.8180339887498948;
      newrange = 1.63*range;
      if (4.0*newrange*newrange*size*(size-1)/2.0 + 2*newrange*outdim < 30000.0) {
        range = newrange;
      }
    } else {
      if (found == 1) {
        dt *= 1.32;
        if (range/1.62 >= 1.0) {
          range /= 1.62;
        }
      } else {
        double newrange;

        dt *= 0.8180339887498948;
        newrange = 1.63*range;
        if (4.0*newrange*newrange*size*(size-1)/2.0 + 2*newrange*outdim < 30000.0) {
          range = newrange;
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
