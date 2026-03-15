#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>

#include "floatType.h"
#if 0
#include "matrix.h"
#endif

extern void allocMatrix(REAL*** pm, int dim);
extern void freeMatrix(REAL** m, int dim);

extern void zeroMatrix(REAL** m, int dim);
extern void identityMatrix(REAL** m, int dim);
extern void initMatrix(REAL** m, int col, int dim);

extern void printMatrix(REAL** m, int dim);
extern void printMatrixSum(REAL** m, int col, int dim);

extern void copyMatrix(REAL** mDst, REAL** mSrc, int dim);
extern void diffMatrix(REAL** mDst, REAL** mSrc, int dim);
extern REAL normMatrix(REAL** m, int dim);

extern void orthoNormMatrix(REAL** m, int row, int dim);
extern void balanceMatrix(REAL** m, int col, int dim);
extern void projectMatrix(REAL** m, REAL** mOld, int dim, REAL factor);

extern void orthoNormMatrix2(REAL** m, int* rows, int* cols, int dim);
extern void balanceMatrix2(REAL** m, int* rows, int* cols, int dim, int num,
                          REAL factor);

extern REAL getMaxSum(REAL** m, int col, int dim);
extern REAL getMaxSum2(REAL** pm, int dim, int num);

extern void insertMatrix(REAL** m1, int dim1, int num, REAL** m2, int dim2,
                         int p1, int p2, int p3);

extern void dotProductMatrix(REAL** dp, REAL** m, int dim);

extern REAL** dotprod;

main(int argc, char **argv)
{
  int m,n;
  int i;
  REAL** mat;
  REAL** saveMat;
  char format[80];
  int row;
  REAL oldSize,newSize;
  int totalCount;

  n = 7;

  if (argc > 1) {
    n = atoi(argv[1]);
  }

  m = n-1;

  allocMatrix(&mat,n);
  initMatrix(mat,m,n);

  allocMatrix(&saveMat,n);
  zeroMatrix(saveMat,n);
  
  oldSize = (REAL)0.0;
  totalCount = 0;

#if 0
  srand48(time(NULL));
#endif

#if 1
      fprintf(stderr,"Start:\n");
      printMatrixSum(mat,m,n);
      fprintf(stderr,"\n");
#endif

  while (1) {
    copyMatrix(saveMat,mat,n);
    
    balanceMatrix(mat,m,n);

    diffMatrix(saveMat,mat,n);
    if (normMatrix(saveMat,n) < TOLER) {
      break;
    }

#if 0
      fprintf(stderr,"Balance:\n");
      printMatrixSum(mat,m,n);
      fprintf(stderr,"\n");
#endif

#if 0
    row = drand48() * n;
#else
    row = 0;
#endif
    orthoNormMatrix(mat,row,n);

#if 0
      fprintf(stderr,"Ortho:\n");
      printMatrixSum(mat,m,n);
      fprintf(stderr,"\n");
#endif

    oldSize = newSize;
    newSize = (REAL)1.0/getMaxSum(mat,m,n) - (REAL)1.0;

    if (totalCount % (1000000000 / n) == 0) {
      printMatrixSum(mat,m,n);
      fprintf(stderr,"\n");

      sprintf(format,"%%3d %%3d     %%19.17%s\n",PRINTE);
      fprintf(stderr,format,m,n,newSize);
      fprintf(stderr,"\n");
    } else if (totalCount % (10000000 / n) == 0) {
      sprintf(format,"%%3d %%3d     %%19.17%s\n",PRINTE);
      fprintf(stderr,format,m,n,newSize);
      fprintf(stderr,"\n");
    }

    totalCount++;
  }

  printMatrixSum(mat,m,n);
  fprintf(stderr,"\n");

  fprintf(stderr,"Iterations: %d\n",totalCount);
  fprintf(stderr,"\n");
  sprintf(format,"%%3d %%3d     %%19.17%s\n",PRINTE);
  fprintf(stderr,format,m,n,newSize);
  fprintf(stderr,"\n");

  freeMatrix(mat,n);

  exit(0);
}

#define DEBUG   1
#define DEBUG2  1

REAL** dotprod;

void allocMatrix(REAL*** pm, int dim)
{
  REAL** m;
  int i,j;

  m = (REAL**)malloc(dim*sizeof(REAL *));

  for (i = 0; i < dim; i++) {
    m[i] = (REAL *)malloc(dim*sizeof(REAL));
  }

  *pm = m;
}

void freeMatrix(REAL** m, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    free(m[i]);
  }

  free(m);
}

void zeroMatrix(REAL** m, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      m[i][j] = (REAL)0.0;
    }
  }
}

void identityMatrix(REAL** m, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      m[i][j] = (REAL)0.0;
    }
  }

  for (i = 0; i < dim; i++) {
    m[i][i] = (REAL)1.0;
  }
}

void initMatrix(REAL** m, int col, int dim)
{
  int i,j;

  m[0][0] = SQRT((REAL)0.5);
  for (j = 1; j < dim; j++) {
    m[0][j] = SQRT((REAL)(0.5/(dim-1)));
  }

  m[1][0] = SQRT((REAL)0.5);
  for (j = 1; j < dim; j++) {
    m[1][j] = -SQRT((REAL)(0.5/(dim-1)));
  }

  for (i = 2; i < dim; i++) {
    for (j = 0; j < i-1; j++) {
      m[i][j] = (REAL)0.0;
    }
    m[i][i-1]   =  SQRT((REAL)(      0.9/(dim-i)));
    for (j = i; j < dim; j++) {
      m[i][j] = -SQRT((REAL)((1.0 - 0.9/(dim-i))/(dim-i)));
    }
  }
}

void printMatrix(REAL** m, int dim)
{
  int i,j;
  char format[80];

  sprintf(format,"%%19.16%s ",PRINTF);

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      fprintf(stderr,format,m[i][j]);
    }
    fprintf(stderr,"\n");
  }
}

void printMatrixSum(REAL** m, int col, int dim)
{
  int i,j;
  REAL sum;
  char format[80];

  sprintf(format,"%%19.16%s ",PRINTF);

  for (i = 0; i < dim; i++) {
    sum = (REAL)0.0;
    for (j = 0; j < dim; j++) {
      fprintf(stderr,format,m[i][j]);

      if (j < col) {
        sum += fabs(m[i][j]);
      }
    }

    fprintf(stderr," - ");
    fprintf(stderr,format,sum);
    fprintf(stderr,"\n");
  }
}

void copyMatrix(REAL** mDst, REAL** mSrc, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mDst[i][j] = mSrc[i][j];
    }
  }
}

void diffMatrix(REAL** mDst, REAL** mSrc, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mDst[i][j] -= mSrc[i][j];
    }
  }
}

REAL normMatrix(REAL** m, int dim)
{
  int i,j;
  REAL max;

  max = (REAL)0.0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if (fabs(m[i][j]) > max) {
        max = fabs(m[i][j]);
      }
    }
  }

  return max;
}

void orthoNormMatrix(REAL** m, int row, int dim)
{
  int i,j,k;
  int ri,rj;

  for (i = dim-1; i >= 1; i--) {
    REAL norm;

    ri = (i + row) % dim;

    if (i < 2) {
      norm = (REAL)0.0;
      for (j = 1; j < dim; j++) {
        norm += m[ri][j]*m[ri][j];
      }
      norm = SQRT((REAL)0.5/norm);

      for (j = 1; j < dim; j++) {
        m[ri][j] *= norm;
      }
    } else {
      norm = (REAL)0.0;
      for (j = 0; j < dim; j++) {
        norm += m[ri][j]*m[ri][j];
      }
      norm = SQRT((REAL)1.0/norm);

      for (j = 0; j < dim; j++) {
        m[ri][j] *= norm;
      }
    }

    for (j = i-1; j >= 0; j--) {
      REAL dot;

      rj = (j + row) % dim;

      dot = (REAL)0.0;
      for (k = 0; k < dim; k++) {
        dot += m[ri][k]*m[rj][k];
      }

      for (k = 0; k < dim; k++) {
        m[rj][k] -= dot*m[ri][k];
      }
    }
  }

  m[0][0] = m[1][0];
  for (j = 1; j < dim; j++) {
    m[0][j] = -m[1][j];
  }
}

void balanceMatrix(REAL** m, int col, int dim)
{
  int i,j;
  REAL sumSum,aver;
  int done;

  sumSum = (REAL)0.0;
  for (i = 0; i < dim; i++) {
    REAL sum;

    sum = (REAL)0.0;
    for (j = 0; j < col; j++) {
      sum += fabs(m[i][j]);
    }

    sumSum += sum;
  }

  aver = sumSum / dim;

  for (i = 0; i < dim; i++) {
    REAL sum,amt;

    for (j = 0; j < i-1; j++) {
      m[i][j] = (REAL)0.0;
    }

    sum = (REAL)0.0;
    for (j = 1; j < col; j++) {
      sum += fabs(m[i][j]);
    }

    while (fabs(sum - aver) > TOLER) {
      if (i < 2) {
        amt = (aver - sum) / (col-1);
        
        for (j = 1; j < col; j++) {
          if (m[i][j] < 0) {
            m[i][j] -= amt;
          } else {
            m[i][j] += amt;
          }
        }
      } else {
        amt = (aver - sum) / (col-i+1);

        for (j = i-1; j < col; j++) {
          if (m[i][j] < 0) {
            m[i][j] -= amt;
          } else {
            m[i][j] += amt;
          }
        }
      }

      sum = (REAL)0.0;
      for (j = 0; j < col; j++) {
        sum += fabs(m[i][j]);
      }
    }
  }
}

void projectMatrix(REAL** m, REAL** mOld, int dim, REAL factor)
{
  int i,j;

  for (i = 0; i < dim-1; i++) {
    for (j = 0; j < dim; j++) {
      m[i][j] = m[i][j] + factor*(m[i][j] - mOld[i][j]);
    }
  }
}

void orthoNormMatrix2(REAL** m, int* rows, int* cols, int dim)
{
  int l;
  for (l = 0; l < 5; l++) {
    int i,j,k;

    for (i = dim-1; i >= 0; i--) {
      REAL sum0,sum1,sum2,sum3;

      if (rows[i] == -1) {
        sum0 = (REAL)0.0;
        sum1 = (REAL)0.0;
        sum2 = (REAL)0.0;
        sum3 = (REAL)0.0;
        for (j = 0; j < dim; j++) {
          if (cols[j] == 0) {
            sum0 += m[i][j]*m[i][j];
          } else if (cols[j] == 1) {
            sum1 += m[i][j]*m[i][j];
          } else {
            sum2 += m[i][j]*m[i][j];
            sum3 += m[i][j];
          }
        }

        if (sum0 + sum1 < (REAL)1.0) {
          REAL delta;
          REAL scale;

          delta = (1 - (sum0+sum1+sum2)) / (2*sum3);

          sum2 = (REAL)0.0;
          for (j = 0; j < dim; j++) {
            if (cols[j] == 2) {
              m[i][j] += delta;
              sum2 += m[i][j]*m[i][j];
            }
          }

          scale = SQRT(((REAL)1.0-(sum0+sum1))/sum2);

          for (j = 0; j < dim; j++) {
            if (cols[j] == 2) {
              m[i][j] *= scale;
            }
          }
        } else {
          REAL scale;

          scale = SQRT(((REAL)1.0-sum0)/(sum1+sum2));

          for (j = 0; j < dim; j++) {
            if (cols[j] == 1 || cols[j] == 2) {
              m[i][j] *= scale;
            }
          }
        }
      }
#if DEBUG2
  fprintf(stderr,"i (normal): %d\n",i);
  printMatrix(m,dim);
  fprintf(stderr,"\n");
  dotProductMatrix(dotprod,m,dim);
  printMatrix(dotprod,dim);
  fprintf(stderr,"\n");
#endif

      if (rows[i] == -1) {
        for (j = i-1; j >= 0; j--) {
          REAL dot,sum;

          dot = (REAL)0.0;
          sum = (REAL)0.0;
          for (k = 0; k < dim; k++) {
            dot += m[i][k]*m[j][k];
            if (cols[k] == 2) {
              sum += m[i][k]*m[i][k];
            }
          }

          for (k = 0; k < dim; k++) {
            if (cols[k] == 2) {
              m[j][k] -= dot*m[i][k]/sum;
            }
          }
#if DEBUG2
  fprintf(stderr,"i (ortho): %d\n",i);
  printMatrix(m,dim);
  fprintf(stderr,"\n");
  dotProductMatrix(dotprod,m,dim);
  printMatrix(dotprod,dim);
  fprintf(stderr,"\n");
#endif
        }
      }
    }

    for (i = 0; i < dim; i++) {
      if (rows[i] >= 0) {
        for (j = rows[i]+1; j < dim; j++) {
          m[i][j] = -m[i+1][j];
        }
      }
    }
#if DEBUG 
  fprintf(stderr,"Orthonormal:\n");
  printMatrix(m,dim);
  fprintf(stderr,"\n");
  dotProductMatrix(dotprod,m,dim);
  printMatrix(dotprod,dim);
  fprintf(stderr,"\n");
#endif
  }
}

void balanceMatrix2(REAL** m, int* rows, int* cols, int dim, int num, REAL factor)
{
  int i,j;
  REAL sumSum,aver;
  int done;

  sumSum = (REAL)0.0;
  for (i = 0; i < dim; i++) {
    REAL sum;

    sum = (REAL)0.0;
    for (j = 0; j < num; j++) {
      sum += fabs(m[i][j]);
    }

    sumSum += sum;
  }

  aver = sumSum / dim;

  for (i = 0; i < dim; i++) {
    REAL sum;

    if (rows[i] == -1) {
      sum = (REAL)0.0;
      for (j = 0; j < num; j++) {
        sum += fabs(m[i][j]);
      }

      if (m[i][cols[i]] == (REAL)0.0) {
        if (i % 2 == 0) {
          m[i][cols[i]] -= factor*(aver-sum);
        } else {
          m[i][cols[i]] += factor*(aver-sum);
        }
      } else
      if (m[i][cols[i]] < (REAL)0.0) {
        m[i][cols[i]] -= factor*(aver-sum);
      } else {
        m[i][cols[i]] += factor*(aver-sum);
      }
    }
  }

  for (i = 0; i < dim; i++) {
    if (rows[i] >= 0) {
      for (j = rows[i]+1; j < dim; j++) {
        m[i][j] = -m[i+1][j];
      }
    }
  }
}

REAL getMaxSum(REAL** m, int col, int dim)
{
  int i,j;
  REAL sum,maxSum;

  maxSum = (REAL)0.0;
  for (i = 0; i < dim; i++) {
    REAL sum;

    sum = (REAL)0.0;
    for (j = 0; j < col; j++) {
      sum += fabs(m[i][j]);
    }

    if (sum > maxSum) {
      maxSum = sum;
    }
  }

  return maxSum;
}

REAL getMaxSum2(REAL** m, int dim, int num)
{
  int i,j;
  REAL sum,maxSum;

  maxSum = (REAL)0.0;
  for (i = 0; i < dim; i++) {
    REAL sum;

    sum = (REAL)0.0;
    for (j = 0; j < num; j++) {
      sum += fabs(m[i][j]);
    }

    if (sum > maxSum) {
      maxSum = sum;
    }
  }

  return maxSum;
}

void insertMatrix(REAL** m1, int dim1, int num, REAL** m2, int dim2,
                  int p1, int p2, int p3)
{
  int i,j;

  for (i = 0; i < dim2; i++) {
    for (j = 0; j < dim2-1; j++) {
      m1[i+p1][j+p2] = m2[i][j];
    }
  }

  for (i = 0; i < dim2; i++) {
    m1[i+p1][p3] = m2[i][dim2-1];
  }
}

void dotProductMatrix(REAL** dp, REAL** m, int dim)
{
  int i,j,k;

  for (i = 0; i < dim; i++) {
    REAL dot;

    for (j = 0; j < dim; j++) {
      dot = (REAL)0.0;
      for (k = 0; k < dim; k++) {
        dot += m[i][k]*m[j][k];
      }

      dp[i][j] = dot;
    }
  }
}
