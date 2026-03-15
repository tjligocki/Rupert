#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "floatType.h"
#include "matrix.h"

#define DEBUG   1
#define DEBUG2  1

main(int argc, char **argv)
{
  int m,n;
  int num;
  int div;
  int m1,m2;
  REAL** mat1;
  REAL** mat2;
  REAL** full;
  int* balanceCols;
  int* orthoNormCols;
  int* rows;
  int i,j;
  int p1,p2,p3;
  char format[80];

  m = 4;
  n = 7;

  if (argc > 2) {
    m = atoi(argv[1]);
    n = atoi(argv[2]);
  }

  if (2*m < n) {
    fprintf(stderr,"m, %d, must be >= n/2, %.1f\n",m,n/2.0);
    exit(1);
  }

  num = n - m;

  for (div = 1; div <= num; div++) {
    if ((div - num + m) % num == 0) {
      break;
    }
  }

  m1 = (div - num + m) / num + 1;
  m2 = m1 + 1;

  fprintf(stderr,"Dim: (%d,%d)\n",m,n);
  fprintf(stderr,"  Num: %d\n",num);
  fprintf(stderr,"  Div: %d\n",div);
  fprintf(stderr,"  M: (%d,%d)\n",m1,m2);
  fprintf(stderr,"\n");

  getPrimeMatrix(&mat1,m1);

  if (num > 1) {
    getPrimeMatrix(&mat2,m2);
  }

  allocMatrix(&full,n);
  zeroMatrix(full,n);

  allocMatrix(&dotprod,n);

  balanceCols = (int *)malloc(n*sizeof(*balanceCols));
  orthoNormCols = (int *)malloc(n*sizeof(*orthoNormCols));

  rows = (int *)malloc(n*sizeof(*rows));

  for (i = 0; i < n; i++) {
    rows[i] = -1;
  }

  p1 = 0;
  p2 = 0;
  p3 = m;
  for (i = 0; i < num; i++) {
    if (i < div) {
      insertMatrix(full,n,m,mat1,m1,p1,p2,p3);

      if (m1 > 2) {
        balanceCols[p1]   = p2+1;
        balanceCols[p1+1] = p2+1;
        for (j = 2; j < m1; j++) {
          balanceCols[p1+j] = p2+j-1;
        }
      } else {
        balanceCols[p1]   = div+1;
        balanceCols[p1+1] = div+1;
      }

      orthoNormCols[p2] = 0;
      for (j = 1; j < m1-1; j++) {
        orthoNormCols[p2+j] = 1;
      }

      rows[p1] = p2;

      p1 += m1;
      p2 += m1-1;
    } else {
      insertMatrix(full,n,m,mat2,m2,p1,p2,p3);

      balanceCols[p1]   = p2+1;
      balanceCols[p1+1] = p2+1;
      for (j = 2; j < m2; j++) {
        balanceCols[p1+j] = p2+j-1;
      }

      orthoNormCols[p2] = 0;
      for (j = 1; j < m2-1; j++) {
        orthoNormCols[p2+j] = 1;
      }

      rows[p1] = p2;

      p1 += m2;
      p2 += m2-1;
    }

    p3++;
  }

  for (i = 0; i < m-1; i++) {
    orthoNormCols[i] = 0;
  }

  orthoNormCols[m-1] = 1;

  for (i = m; i < n; i++) {
    orthoNormCols[i] = 2;
  }

  printMatrix(mat1,m1);
  fprintf(stderr,"\n");

  sprintf(format,"Side: %%19.16%s\n\n",PRINTF);
  fprintf(stderr,format,(REAL)1.0/getMaxSum2(mat1,m1,m1-1));

  if (num > 1) {
    printMatrix(mat2,m2);
    fprintf(stderr,"\n");

    sprintf(format,"Side: %%19.16%s\n\n",PRINTF);
    fprintf(stderr,format,(REAL)1.0/getMaxSum2(mat2,m2,m2-1));
  }

  printMatrix(full,n);
  fprintf(stderr,"\n");

  sprintf(format,"Side: %%19.16%s\n\n",PRINTF);
  fprintf(stderr,format,(REAL)1.0/getMaxSum2(full,n,m));

  for (i = 0; i < 50000; i++) {
#if DEBUG
    fprintf(stderr,"Rows: ");
    for (j = 0; j < n; j++) {
      fprintf(stderr,"%1d ",rows[j]);
    }
    fprintf(stderr,"\n\n");
#endif

#if DEBUG
    fprintf(stderr,"Balance column:\n");
    for (j = 0; j < n; j++) {
      fprintf(stderr,"    %2d: %2d\n",j,balanceCols[j]);
    }
    fprintf(stderr,"\n");
#endif

    balanceMatrix2(full,rows,balanceCols,n,m,(REAL)1.0);

#if DEBUG
    printMatrix(full,n);
    fprintf(stderr,"\n");

    sprintf(format,"Side: %%19.16%s\n\n",PRINTF);
    fprintf(stderr,format,(REAL)1.0/getMaxSum2(full,n,m));
#else
    sprintf(format,"Side: %%19.16%s ... ",PRINTF);
    fprintf(stderr,format,(REAL)1.0/getMaxSum2(full,n,m));
#endif

#if DEBUG
    fprintf(stderr,"OrthoNorm columns: ");
    for (j = 0; j < n; j++) {
      fprintf(stderr,"%1d ",orthoNormCols[j]);
    }
    fprintf(stderr,"\n\n");
#endif

    orthoNormMatrix2(full,rows,orthoNormCols,n);

#if DEBUG
    printMatrix(full,n);
    fprintf(stderr,"\n");

    sprintf(format,"Side: %%19.16%s\n\n",PRINTF);
    fprintf(stderr,format,(REAL)1.0/getMaxSum2(full,n,m));
#else
    sprintf(format,"Side: %%19.16%s\r",PRINTF);
    fprintf(stderr,format,(REAL)1.0/getMaxSum2(full,n,m));
#endif

#if DEBUG
    dotProductMatrix(dotprod,full,n);
    printMatrix(dotprod,n);
    fprintf(stderr,"\n");
#endif
  }

  fprintf(stderr,"\n\n");

  printMatrix(full,n);
  fprintf(stderr,"\n");

  sprintf(format,"Side: %%19.16%s\n\n",PRINTF);
  fprintf(stderr,format,(REAL)1.0/getMaxSum2(full,n,m));

  dotProductMatrix(dotprod,full,n);
  printMatrix(dotprod,n);
  fprintf(stderr,"\n");

  exit(0);
}
