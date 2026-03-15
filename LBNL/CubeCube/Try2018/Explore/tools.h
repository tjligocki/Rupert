#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern inline double *matrixIdentity(int n);
extern inline void    matrixCopy(double *matrixOut, double *matrixIn, int n);
extern inline void    matrixRotLeft(double *matrix, double angle, int i1, int i2, int n);
extern inline void    matrixRotRight(double *matrix, double angle, int i1, int i2, int n);
extern inline double *matrixSums(double *matrix, int m, int n);
extern inline void    matrixPrint(FILE *fd, char *format, double *matrix, int n);
extern inline void    matrixFree(double *matrix);

extern inline double *rupertMatrix(double *x, int m, int n);
extern inline double *rupertRep(double *matrix, int m, int n);
extern inline double *rupertSums(double *x, int m, int n);
extern inline double  rupert(double *x, int m, int n);

extern inline int *indexMake(int dim);
extern inline int  indexValid(int *index, int n);
extern inline void indexInc(int *index, int dim, int n);
extern inline void indexPrint(FILE *fd, char *format, int *index, int dim, int final);
extern inline void indexFree(int *index);

extern inline double *pointMake(int dim);
extern inline void    pointPrint(FILE *fd, char *format, double *point, int dim, int final);
extern inline void    pointFree(double *point);

extern inline double brentZero(double x1, double x2, double (*f)(double), double tol);
extern inline double brentMin(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
