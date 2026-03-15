#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "floatType.h"

void allocMatrix(REAL*** pm, int dim);

void initMatrix(REAL** pm, int dim);

void printMatrix(REAL** m, int dim);
void printFullMatrix(REAL** m, int dim);

void copyMatrix(REAL** mDst, REAL** mSrc, int dim);
void diffMatrix(REAL** mDst, REAL** mSrc, int dim);
REAL normMatrix(REAL** m, int dim);

void orthoNormMatrix(REAL** m, int dim);
void balanceMatrix(REAL** m, int dim, REAL factor);
void projectMatrix(REAL** m, REAL** mOld, int dim, REAL factor);

REAL getMaxSum(REAL** m, int dim);

main(int argc, char **argv)
{
  int dim;
  REAL** m;
  REAL** mLast;
  REAL** mOld;
  int count;
  REAL factor;
  int copyOnce;
  char format[80];

  dim = 3;

  if (argc > 1) {
    dim = atoi(argv[1]);
  }

  fprintf(stderr,"Dim: %d\n",dim);
  fprintf(stderr,"\n");

  allocMatrix(&m,dim);
  initMatrix(m,dim);

  allocMatrix(&mLast,dim);
  copyMatrix(mLast,m,dim);

  allocMatrix(&mOld,dim);
  copyOnce = 0;

  count = 0;
  while (1) {
    if ((count % 10) == 5 && copyOnce == 1 && count > 100000) {
      projectMatrix(m,mOld,dim,(REAL)1.00);
    }

    orthoNormMatrix(m,dim);

    if ((count % 1000000) == 0) {
      fprintf(stderr,"%d:\n",count);
      printMatrix(m,dim);
      fprintf(stderr,"\n");
    }

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
    balanceMatrix(m,dim,factor);

    count++;
  }

  fprintf(stderr,"%d iterations\n\n",count);

  printFullMatrix(m,dim);

  sprintf(format,"\nSide: %%19.16%s\n",PRINTF);
  fprintf(stderr,format,(REAL)1.0/getMaxSum(m,dim));
}

void allocMatrix(REAL*** pm, int dim)
{
  REAL** m;
  int i,j;

  m = (REAL**)malloc((dim-1)*sizeof(REAL *));

  for (i = 0; i < dim-1; i++) {
    m[i] = (REAL *)malloc(dim*sizeof(REAL));
  }

  *pm = m;
}

void initMatrix(REAL** m, int dim)
{
  int i,j;

  for (i = 0; i < dim-1; i++) {
    for (j = 0; j < i; j++) {
      m[i][j] = (REAL)0.0;
    }

    m[i][i] = (REAL)1.0;

    for (j = i+1; j < dim; j++) {
      m[i][j] = (REAL)-1.0;
    }
  }

  m[0][0] = SQRT((REAL)2.0)/(REAL)2.0;
}

void printMatrix(REAL** m, int dim)
{
  int i,j;
  char format[80];

  sprintf(format,"%%19.16%s ",PRINTF);

  for (i = 0; i < dim-1; i++) {
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

  for (i = 0; i < dim-1; i++) {
    for (j = 0; j < dim; j++) {
      mDst[i][j] = mSrc[i][j];
    }
  }
}

void diffMatrix(REAL** mDst, REAL** mSrc, int dim)
{
  int i,j;

  for (i = 0; i < dim-1; i++) {
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

  for (i = 0; i < dim-1; i++) {
    for (j = 0; j < dim; j++) {
      if (fabs(m[i][j]) > max) {
        max = fabs(m[i][j]);
      }
    }
  }

  return max;
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
