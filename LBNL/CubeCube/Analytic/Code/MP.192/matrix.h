#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "floatType.h"

extern void allocMatrix(REAL*** pm, int dim);
extern void freeMatrix(REAL** m, int dim);

extern void zeroMatrix(REAL** m, int dim);
extern void identityMatrix(REAL** m, int dim);
extern void initMatrix(REAL** m, int dim);

extern void printMatrix(REAL** m, int dim);
extern void printFullMatrix(REAL** m, int dim);

extern void copyMatrix(REAL** mDst, REAL** mSrc, int dim);
extern void diffMatrix(REAL** mDst, REAL** mSrc, int dim);
extern REAL normMatrix(REAL** m, int dim);

extern void getPrimeMatrix(REAL*** pm, int dim);
extern void orthoNormPrimeMatrix(REAL** m, int dim);
extern void balancePrimeMatrix(REAL** m, int dim, REAL factor);
extern void projectPrimeMatrix(REAL** m, REAL** mOld, int dim, REAL factor);

extern void orthoNormMatrix(REAL** m, int dim);
extern void balanceMatrix(REAL** m, int dim, REAL factor);
extern void projectMatrix(REAL** m, REAL** mOld, int dim, REAL factor);

extern void orthoNormMatrix2(REAL** m, int* rows, int* cols, int dim);
extern void balanceMatrix2(REAL** m, int* rows, int* cols, int dim, int num,
                          REAL factor);

extern REAL getMaxSum(REAL** m, int dim);
extern REAL getMaxSum2(REAL** pm, int dim, int num);

extern void insertMatrix(REAL** m1, int dim1, int num, REAL** m2, int dim2,
                         int p1, int p2, int p3);

extern void dotProductMatrix(REAL** dp, REAL** m, int dim);

extern REAL** dotprod;
