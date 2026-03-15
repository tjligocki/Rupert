#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"

int main(int argc, char **argv)
{
  int m,n;
  int dim;
  int i;
  int r1,r2;
  double a;
  double *x;
  double *matrix;
  double *sums;
  double done = 0;

  r1 = atoi(argv[1]);
  r2 = atoi(argv[2]);
  a  = atof(argv[3]);

  m = atoi(argv[4]);
  n = atoi(argv[5]);

  if (r1 < 0 || r2 < 0 || r1 >= n || r2 >= n) {
    fprintf(stderr,"The rotation indices, %d %d, must both be non-negative and less than n = %d\n",r1,r2,n);
    exit(1);
  }

  dim = (m*(2*n-(m+1)))/2;

  x = pointMake(dim);

  i = fread(x,sizeof(*x),dim,stdin);

  pointPrint(stderr,"%8.5lf",x,dim,1);
  fprintf(stderr,"\n");

  matrix = rupertMatrix(x,m,n);

  matrixPrint(stderr,"%8.5lf",matrix,n);
  fprintf(stderr,"\n");

  sums = matrixSums(matrix,m,n);
  pointPrint(stderr,"%8.5lf",sums,n,1);
  fprintf(stderr,"\n");
  pointFree(sums);

  fprintf(stderr,"\n");

  matrixRotLeft(matrix,a,r1,r2,n);

  x = rupertRep(matrix,m,n);

  pointPrint(stderr,"%8.5lf",x,dim,1);
  fprintf(stderr,"\n");

  matrixPrint(stderr,"%8.5lf",matrix,n);
  fprintf(stderr,"\n");

  sums = matrixSums(matrix,m,n);
  pointPrint(stderr,"%8.5lf",sums,n,1);
  fprintf(stderr,"\n");
  pointFree(sums);

  fwrite(x,sizeof(*x),dim,stdout);

  pointFree(x);
  matrixFree(matrix);

  exit(0);
}
