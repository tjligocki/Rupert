#include "matrix.h"

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

void initMatrix(REAL** m, int dim)
{
  int i,j;

  for (i = 1; i < dim; i++) {
    for (j = 0; j < i-1; j++) {
      m[i][j] = (REAL)0.0;
    }

    m[i][i-1] = (REAL)1.0;

    for (j = i; j < dim; j++) {
      m[i][j] = (REAL)-1.0;
    }
  }

  m[1][0] = SQRT((REAL)2.0)/(REAL)2.0;

  m[0][0] = m[1][0];
  for (j = 1; j < dim; j++) {
    m[0][j] = -m[1][j];
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

void printFullMatrix(REAL** m, int dim)
{
  int i,j;
  char format[80];

  sprintf(format,"%%19.16%s ",PRINTF);

  fprintf(stderr,format,m[0][0]);
  for (j = 1; j < dim; j++) {
    fprintf(stderr,format,-m[0][j]);
  }
  fprintf(stderr,"\n");

  for (i = 0; i < dim-1; i++) {
    for (j = 0; j < dim; j++) {
      fprintf(stderr,format,m[i][j]);
    }
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

void getPrimeMatrix(REAL*** pm, int dim)
{
  REAL** m;
  REAL** mLast;
  REAL** mOld;
  int count;
  REAL factor;
  int copyOnce;

  allocMatrix(&m,dim);
  initMatrix(m,dim);

  allocMatrix(&mLast,dim);
  copyMatrix(mLast,m,dim);

  allocMatrix(&mOld,dim);
  copyOnce = 0;

  count = 0;
  while (1) {
    if ((count % 10) == 5 && copyOnce == 1 && count > 100000) {
      projectPrimeMatrix(m,mOld,dim,(REAL)1.00);
    }

    orthoNormPrimeMatrix(m,dim);

    if ((count % 100000) == 1) {
      diffMatrix(mLast,m,dim);

      if (normMatrix(mLast,dim) < TOLER) {
        break;
      }
    }

    if ((count % 100000) == 0) {
      copyMatrix(mLast,m,dim);
    }

    if ((count % 10) == 5) {
      copyOnce = 1;
      copyMatrix(mOld,m,dim);
    }

    factor = (REAL)1.0 + (REAL)0.5*exp(-count/(REAL)1000000000.0);
    balancePrimeMatrix(m,dim,factor);

    count++;
  }

  *pm = m;

  freeMatrix(mLast,dim);
  freeMatrix(mOld,dim);
}

void orthoNormPrimeMatrix(REAL** m, int dim)
{
  int i,j,k;

  for (i = dim-1; i >= 1; i--) {
    if (i == 1) {
      REAL norm;

      norm = (REAL)0.0;
      for (j = 1; j < dim; j++) {
        norm += m[i][j]*m[i][j];
      }
      norm = SQRT((REAL)0.5/norm);

      for (j = 1; j < dim; j++) {
        m[i][j] *= norm;
      }
    } else {
      REAL norm;

      norm = (REAL)0.0;
      for (j = 0; j < dim; j++) {
        norm += m[i][j]*m[i][j];
      }
      norm = SQRT((REAL)1.0/norm);

      for (j = 0; j < dim; j++) {
        m[i][j] *= norm;
      }
    }

    for (j = i-1; j >= 1; j--) {
      REAL dot;

      dot = (REAL)0.0;
      for (k = i-1; k < dim; k++) {
        dot += m[i][k]*m[j][k];
      }

      for (k = i-1; k < dim; k++) {
        m[j][k] -= dot*m[i][k];
      }
    }
  }

  m[0][0] = m[1][0];
  for (j = 1; j < dim; j++) {
    m[0][j] = -m[1][j];
  }
}

void balancePrimeMatrix(REAL** m, int dim, REAL factor)
{
  int i,j;
  REAL sumSum,aver;
  int done;

  sumSum = (REAL)0.0;
  for (i = 0; i < dim; i++) {
    REAL sum;

    sum = (REAL)0.0;
    for (j = 0; j < dim-1; j++) {
      sum += fabs(m[i][j]);
    }

    sumSum += sum;
  }

  aver = sumSum / dim;

  for (i = 1; i < dim; i++) {
    REAL sum;

    sum = (REAL)0.0;
    for (j = 0; j < dim-1; j++) {
      sum += fabs(m[i][j]);
    }

    if (i == 1) {
      if (m[i][i] < 0) {
        m[i][i] -= (aver-sum);
      } else {
        m[i][i] += (aver-sum);
      }
    } else {
      if (m[i][i-1] < 0) {
        m[i][i-1] -= factor*(aver-sum);
      } else {
        m[i][i-1] += factor*(aver-sum);
      }
    }
  }

  m[0][0] = m[1][0];
  for (j = 1; j < dim; j++) {
    m[0][j] = -m[1][j];
  }
}

void projectPrimeMatrix(REAL** m, REAL** mOld, int dim, REAL factor)
{
  int i,j;

  for (i = 1; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      m[i][j] = m[i][j] + factor*(m[i][j] - mOld[i][j]);
    }
  }

  m[0][0] = m[1][0];
  for (j = 1; j < dim; j++) {
    m[0][j] = -m[1][j];
  }
}

void orthoNormMatrix(REAL** m, int dim)
{
  int i,j,k;

  for (i = dim-2; i >= 0; i--) {
    if (i == 0) {
      REAL norm;

      norm = (REAL)0.0;
      for (j = 1; j < dim; j++) {
        norm += m[i][j]*m[i][j];
      }
      norm = SQRT((REAL)0.5/norm);

      for (j = 1; j < dim; j++) {
        m[i][j] *= norm;
      }
    } else {
      REAL norm;

      norm = (REAL)0.0;
      for (j = 0; j < dim; j++) {
        norm += m[i][j]*m[i][j];
      }
      norm = SQRT((REAL)1.0/norm);

      for (j = 0; j < dim; j++) {
        m[i][j] *= norm;
      }
    }

    for (j = i-1; j >= 0; j--) {
      REAL dot;

      dot = (REAL)0.0;
      for (k = i; k < dim; k++) {
        dot += m[i][k]*m[j][k];
      }

      for (k = i; k < dim; k++) {
        m[j][k] -= dot*m[i][k];
      }
    }
  }
}

void balanceMatrix(REAL** m, int dim, REAL factor)
{
  int i,j;
  REAL sumSum,aver;
  int done;

  sumSum = (REAL)0.0;
  for (i = 0; i < dim-1; i++) {
    REAL sum;

    sum = fabs(m[i][i]);
    for (j = i+1; j < dim-1; j++) {
      sum += fabs(m[i][j]);
    }

    if (i == 0) {
      sumSum += 2*sum;
    } else {
      sumSum += sum;
    }
  }

  aver = sumSum / dim;

  for (i = 0; i < dim-1; i++) {
    REAL sum;

    sum = fabs(m[i][i]);
    for (j = i+1; j < dim-1; j++) {
      sum += fabs(m[i][j]);
    }

    if (i == 0) {
      if (m[i][i+1] < 0) {
        m[i][i+1] -= (aver-sum);
      } else {
        m[i][i+1] += (aver-sum);
      }
    } else {
      if (m[i][i] < 0) {
        m[i][i] -= factor*(aver-sum);
      } else {
        m[i][i] += factor*(aver-sum);
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

REAL getMaxSum(REAL** m, int dim)
{
  int i,j;
  REAL sum,maxSum;

  maxSum = (REAL)0.0;
  for (i = 0; i < dim-1; i++) {
    REAL sum;

    sum = (REAL)0.0;
    for (j = i; j < dim-1; j++) {
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
