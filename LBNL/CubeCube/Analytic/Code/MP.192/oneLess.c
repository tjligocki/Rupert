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
extern void normMatrix(REAL* norm, REAL** m, int dim);

extern void orthoNormMatrix(REAL** m, int row, int dim);
extern void balanceMatrix(REAL** m, int col, int dim);
extern void projectMatrix(REAL** m, REAL** mOld, int dim, REAL factor);

#if 0
extern void orthoNormMatrix2(REAL** m, int* rows, int* cols, int dim);
extern void balanceMatrix2(REAL** m, int* rows, int* cols, int dim, int num,
                          REAL factor);
#endif

extern void getMaxSum(REAL* maxSum, REAL** m, int col, int dim);
#if 0
extern void getMaxSum2(REAL* maxSum, REAL** pm, int dim, int num);
#endif

#if 0
extern void insertMatrix(REAL** m1, int dim1, int num, REAL** m2, int dim2,
                         int p1, int p2, int p3);
#endif

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
  REAL norm,sum;
  int totalCount;

  n = 7;

  if (argc > 1) {
    n = atoi(argv[1]);
  }

  m = n-1;

  mpf_init2(oldSize,PRES);
  mpf_init2(newSize,PRES);

  mpf_init2(norm,PRES);
  mpf_init2(sum,PRES);

  allocMatrix(&mat,n);
  initMatrix(mat,m,n);

  allocMatrix(&saveMat,n);
  zeroMatrix(saveMat,n);
  
  mpf_set_d(oldSize,0.0);
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
    
    balanceMatrix(mat,m,n);

    copyMatrix(saveMat,mat,n);
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

    diffMatrix(saveMat,mat,n);

    normMatrix(&norm,saveMat,n);
#if 0
      fprintf(stderr,"Norm     ");
      mpf_out_str(stderr,10,0,norm);
      fprintf(stderr,"\n");
      fprintf(stderr,"\n");
#endif
    if (mpf_cmp_d(norm,TOLER) < 0) {
      break;
    }

#if 0
      fprintf(stderr,"Ortho:\n");
      printMatrixSum(mat,m,n);
      fprintf(stderr,"\n");
#endif

    mpf_set(oldSize,newSize);

    getMaxSum(&sum,mat,m,n);
    mpf_ui_div(newSize,1,sum);
    mpf_sub_ui(newSize,newSize,1);

#if 0
      fprintf(stderr,"%3d %3d     ",m,n);
      mpf_out_str(stderr,10,0,newSize);
      fprintf(stderr,"\n");
      fprintf(stderr,"\n");
#endif
    if (totalCount % (1000000 / n) == 0) {
      printMatrixSum(mat,m,n);
      fprintf(stderr,"\n");

      fprintf(stderr,"%3d %3d     ",m,n);
      mpf_out_str(stderr,10,0,newSize);
      fprintf(stderr,"\n");
      fprintf(stderr,"\n");
    }

    totalCount++;
  }

  printMatrixSum(mat,m,n);
  fprintf(stderr,"\n");

  fprintf(stderr,"%3d %3d     ",m,n);
  mpf_out_str(stderr,10,0,newSize);
  fprintf(stderr,"\n");
  fprintf(stderr,"\n");

  mpf_clear(oldSize);
  mpf_clear(newSize);

  mpf_clear(norm);
  mpf_clear(sum);

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

    for (j = 0; j < dim; j++) {
      mpf_init2(m[i][j],PRES);
    }
  }

  *pm = m;
}

void freeMatrix(REAL** m, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mpf_clear(m[i][j]);
    }

    free(m[i]);
  }

  free(m);
}

void zeroMatrix(REAL** m, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mpf_set_d(m[i][j],0.0);
    }
  }
}

void identityMatrix(REAL** m, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mpf_set_d(m[i][j],0.0);
    }
  }

  for (i = 0; i < dim; i++) {
    mpf_set_d(m[i][i],1.0);
  }
}

void initMatrix(REAL** m, int col, int dim)
{
  int i,j;

  for (j = 0; j < dim; j++) {
    mpf_set_d(m[0][j],1.0);
  }

  for (i = 1; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if (j < i-1) {
        mpf_set_d(m[i][j],0.0);
      } else if (j == i-1) {
        mpf_set_d(m[i][j],1.0);
      } else {
        mpf_set_d(m[i][j],-1.0);
      }
    }
  }
}

void printMatrix(REAL** m, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mpf_out_str(stderr,10,0,m[i][j]);
      fprintf(stderr," ");
    }
    fprintf(stderr,"\n");
  }
}

void printMatrixSum(REAL** m, int col, int dim)
{
  int i,j;
  REAL sum,value;

  mpf_init2(sum,PRES);
  mpf_init2(value,PRES);

  for (i = 0; i < dim; i++) {
    mpf_set_d(sum,0.0);
    for (j = 0; j < dim; j++) {
      mpf_out_str(stderr,10,0,m[i][j]);
      fprintf(stderr," ");

      if (j < col) {
        mpf_abs(value,m[i][j]);
        mpf_add(sum,sum,value);
      }
    }

    fprintf(stderr," - ");
    mpf_out_str(stderr,10,0,sum);
    fprintf(stderr,"\n");
  }

  mpf_clear(sum);
  mpf_clear(value);
}

void copyMatrix(REAL** mDst, REAL** mSrc, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mpf_set(mDst[i][j],mSrc[i][j]);
    }
  }
}

void diffMatrix(REAL** mDst, REAL** mSrc, int dim)
{
  int i,j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mpf_sub(mDst[i][j],mDst[i][j],mSrc[i][j]);
    }
  }
}

void normMatrix(REAL* norm, REAL** m, int dim)
{
  int i,j;
  REAL value;

  mpf_init2(value,PRES);

  mpf_set_d(*norm,0.0);

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mpf_abs(value,m[i][j]);
      if (mpf_cmp(value,*norm) > 0) {
        mpf_set(*norm,value);
      }
    }
  }

  mpf_clear(value);
}

void orthoNormMatrix(REAL** m, int row, int dim)
{
  int i,j,k;
  int ri,rj;
  REAL dot;
  REAL prod;

  mpf_init2(dot,PRES);
  mpf_init2(prod,PRES);

  for (i = dim-1; i >= 0; i--) {
    ri = (i + row) % dim;

    mpf_set_d(dot,0.0);
    for (j = 0; j < dim; j++) {
      mpf_mul(prod,m[ri][j],m[ri][j]);
      mpf_add(dot,dot,prod);
    }
    mpf_sqrt(dot,dot);
    mpf_ui_div(dot,1,dot);

    for (j = 0; j < dim; j++) {
      mpf_mul(m[ri][j],m[ri][j],dot);
    }

    for (j = i-1; j >= 0; j--) {
      rj = (j + row) % dim;

      mpf_set_d(dot,0.0);
      for (k = 0; k < dim; k++) {
        mpf_mul(prod,m[ri][k],m[rj][k]);
        mpf_add(dot,dot,prod);
      }

      for (k = 0; k < dim; k++) {
        mpf_mul(prod,dot,m[ri][k]);
        mpf_sub(m[rj][k],m[rj][k],prod);
      }
    }
  }

  mpf_clear(dot);
  mpf_clear(prod);
}

void balanceMatrix(REAL** m, int col, int dim)
{
  int i,j;
  REAL sumSum,aver;
  REAL sum,amt;
  REAL value,diff;
  int done;

  mpf_init2(sumSum,PRES);
  mpf_init2(aver,PRES);

  mpf_init2(sum,PRES);
  mpf_init2(amt,PRES);

  mpf_init2(value,PRES);
  mpf_init2(diff,PRES);

  mpf_set_d(sumSum,0.0);
  for (i = 0; i < dim; i++) {
    mpf_set_d(sum,0.0);
    for (j = 0; j < col; j++) {
      mpf_abs(value,m[i][j]);
      mpf_add(sum,sum,value);
    }

    mpf_add(sumSum,sumSum,sum);
  }

  mpf_div_ui(aver,sumSum,dim);

  for (i = 0; i < dim; i++) {
    for (j = 0; j < i-1; j++) {
      mpf_set_d(m[i][j],0.0);
    }

    mpf_set_d(sum,0.0);
    for (j = 0; j < col; j++) {
      mpf_abs(value,m[i][j]);
      mpf_add(sum,sum,value);
    }

    mpf_sub(diff,aver,sum);
    mpf_abs(value,diff);

    while (mpf_cmp_d(value,TOLER) > 0) {
      if (i == 0) {
        mpf_div_ui(amt,diff,col);
        
        for (j = 0; j < col; j++) {
          if (mpf_cmp_d(m[i][j],0.0) < 0) {
            mpf_sub(m[i][j],m[i][j],amt);
          } else {
            mpf_add(m[i][j],m[i][j],amt);
          }
        }
      } else {
        mpf_div_ui(amt,diff,col-i+1);

        for (j = i-1; j < col; j++) {
          if (mpf_cmp_d(m[i][j],0.0) < 0) {
            mpf_sub(m[i][j],m[i][j],amt);
          } else {
            mpf_add(m[i][j],m[i][j],amt);
          }
        }
      }

      mpf_set_d(sum,0.0);
      for (j = 0; j < col; j++) {
        mpf_abs(value,m[i][j]);
        mpf_add(sum,sum,value);
      }

      mpf_sub(diff,aver,sum);
      mpf_abs(value,diff);
    }
  }

  mpf_clear(sumSum);
  mpf_clear(aver);

  mpf_clear(sum);
  mpf_clear(amt);

  mpf_clear(value);
  mpf_clear(diff);
}

void projectMatrix(REAL** m, REAL** mOld, int dim, REAL factor)
{
  int i,j;
  REAL inc;

  mpf_init2(inc,PRES);

  for (i = 0; i < dim-1; i++) {
    for (j = 0; j < dim; j++) {
      mpf_sub(inc,m[i][j],mOld[i][j]);
      mpf_mul(inc,inc,factor);

      mpf_add(m[i][j],m[i][j],inc);
    }
  }

  mpf_clear(inc);
}

#if 0
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
#endif

void getMaxSum(REAL* maxSum, REAL** m, int col, int dim)
{
  int i,j;
  REAL sum;
  REAL value;

  mpf_init2(sum,PRES);
  mpf_init2(value,PRES);

  mpf_set_d(*maxSum,0.0);
  for (i = 0; i < dim; i++) {

    mpf_set_d(sum,0.0);
    for (j = 0; j < col; j++) {
      mpf_abs(value,m[i][j]);
      mpf_add(sum,sum,value);
    }

    if (mpf_cmp(sum,*maxSum) > 0) {
      mpf_set(*maxSum,sum);
    }
  }

  mpf_clear(sum);
  mpf_clear(value);
}

#if 0
void getMaxSum2(REAL** m, int dim, int num)
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
#endif

#if 0
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
#endif

void dotProductMatrix(REAL** dp, REAL** m, int dim)
{
  int i,j,k;
  REAL dot;
  REAL prod;

  mpf_init2(dot,PRES);
  mpf_init2(prod,PRES);

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      mpf_set_d(dot,0.0);
      for (k = 0; k < dim; k++) {
        mpf_mul(prod,m[i][k],m[j][k]);
        mpf_add(dot,dot,prod);
      }

      mpf_set(dp[i][j],dot);
    }
  }

  mpf_clear(dot);
  mpf_clear(prod);
}
