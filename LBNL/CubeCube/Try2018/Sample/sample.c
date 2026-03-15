#include <stdio.h>
#include <stdlib.h>
#include <math.h>

inline double *matrixIdentity(int n)
{
  double *matrix,*scan;
  int i,j;

  matrix = (double *)malloc(n*n*sizeof(*matrix));

  scan = matrix;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i == j) {
        *scan = 1.0;
      } else {
        *scan = 0.0;
      }

      scan++;
    }
  }

  return matrix;
}

inline void matrixPrint(char *format, double *matrix, int n)
{
  double *scan;
  int i,j;

  scan = matrix;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      printf(format,*scan);
      printf(" ");
      scan++;
    }
    printf("\n");
  }
}

inline void matrixFree(double *matrix)
{
  free(matrix);
}

inline double *rupertMatrix(double *x, int m, int n)
{
  double *gRot;
  int index;
  int i,j,k;
  double maxRow;

  gRot = matrixIdentity(n);

  index = 0;

  for (i = 0; i < m; i++) {
    for (j = i+1; j < n; j++) {
      double cx,sx;
      double ti,tj;

      cx = cos(x[index]);
      sx = sin(x[index]);

      for (k = 0; k < n; k++) {
        int it,jt;

        it = i*n + k;
        jt = j*n + k;

        ti = gRot[it];
        tj = gRot[jt];

        gRot[it] = cx*ti - sx*tj;
        gRot[jt] = sx*ti + cx*tj;
      }

      index++;
    }
  }

  return gRot;
}

inline double rupert(double *x, int m, int n)
{
  double *gRot;
  int index;
  int i,j;
  double maxRow;

  gRot = rupertMatrix(x,m,n);

  maxRow = 0.0;

  for (i = 0; i < n; i++) {
    double curRow = 0.0;

    index = i*n;

    for (j = 0; j < m; j++) {
      curRow += fabs(gRot[index]);
      index++;
    }

    if (curRow > maxRow) {
      maxRow = curRow;
    }
  }

  matrixFree(gRot);

  return maxRow;
}

inline int *indexMake(int dim)
{
  int *index,*scan;
  int i;

  index = (int *)malloc(dim*sizeof(*index));

  scan = index;
  for (i = 0; i < dim; i++)
  {
    *scan = 0;
    scan++;
  }

  return index;
}

inline int indexValid(int *index, int n)
{
  return (index[0] < n);
}

inline void indexInc(int *index, int dim, int n)
{
  int i;

  for (i = dim-1; i >= 0; i--)
  {
    index[i]++;

    if (index[i] == n) {
      if (i > 0) {
        index[i] = 0;
      }
    } else {
      break;
    }
  }
}

inline void indexPrint(char *format, int *index, int dim)
{
  int i;

  for (i = 0; i < dim; i++) {
    printf(format,index[i]);
    printf(" ");
  }
  printf("\n");
}

void indexFree(int *index)
{
  free(index);
}

inline void pointPrint(char *format, double *point, int dim)
{
  int i;

  for (i = 0; i < dim; i++) {
    printf(format,point[i]);
    printf(" ");
  }
  printf("\n");
}

inline double *minimize(double *minp, double range, int np, int m, int n)
{
  int dim;
  int *index;
  double *x;
  double dp;
  double minimum;
  double *minx;

  dp = range / (np-1);

  dim = (m*(2*n-(m+1)))/2;

  index = indexMake(dim);

  x    = (double *)malloc(dim*sizeof(*x));
  minx = (double *)malloc(dim*sizeof(*minx));

  minimum = 1.0;

  while (indexValid(index,np)) {
    int i;
    double cur;

    for (i = 0; i < dim; i++) {
      x[i] = dp * index[i] + minp[i];
    }

    cur = rupert(x,m,n);

    if (cur < minimum) {
      int j;
      minimum = cur;

      for (j = 0; j < dim; j++) {
        minx[j] = x[j];
      }
    }

    indexInc(index,dim,np);
  }

  indexFree(index);
  free(x);
  
  return(minx);
}

int main(int argc, char **argv)
{
  int m,n;
  int dim;
  int np;
  double *minp;
  int i;
  double *minx;
  double *minMatrix;
  double minimum;
  double range;
  int zooms;
  int z;
  int flates;
  int f;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  dim = (m*(2*n-(m+1)))/2;

  np = atoi(argv[3]);

  minp = (double *)malloc(dim*sizeof(*minp));

  zooms  = atoi(argv[4]);
  flates = atoi(argv[5]);

  range = M_PI/2;

  for (i = 0; i < dim; i++) {
    minp[i] = 0.0;
  }

  minx = minimize(minp,range,np,m,n);

  printf("%11.5le -> %20.15lf\n",range,1.0/rupert(minx,m,n));

  for (f = 0; f < flates; f++) {
    double subrange;
    subrange = range;

    for (z = 1; z < zooms; z++) {
      subrange = 2*subrange / (np-1);

      for (i = 0; i < dim; i++) {
        minp[i] = minx[i] - subrange/2;
      }

      free(minx);

      minx = minimize(minp,subrange,np,m,n);

      printf("%11.5le -> %20.15lf\n",subrange,1.0/rupert(minx,m,n));
    }

    printf("\n");

    range = range / 10;
  }

  free(minx);
  free(minp);

  exit(0);
}
