#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"

double *gmatrixRef;
double *gmatrixTmp;
double *gsums;
int gimin,gimax;
double gsmin,gsmax;
int gm,gn;

double diffSum(double a)
{
  double gap;

  matrixCopy(gmatrixTmp,gmatrixRef,gn);

  matrixRotLeft(gmatrixTmp,a,gimin,gimax,gn);

  gsums = matrixSums(gmatrixTmp,gm,gn);

  gsmin = gsums[gimin];
  gsmax = gsums[gimax];

  gap = gsmax - gsmin;

  pointFree(gsums);

  return(gap);
}

double printSum(double a)
{
  double gap;

  matrixCopy(gmatrixTmp,gmatrixRef,gn);

  matrixRotLeft(gmatrixTmp,a,gimin,gimax,gn);

  gsums = matrixSums(gmatrixTmp,gm,gn);

  fprintf(stderr,"   %8.5lf %8.5lf %8.5lf\n",a,gsums[gimax],gsums[gimin]);

  pointFree(gsums);

  return(gap);
}

int matchRowSums(double *matrix, int m, int n)
{
  double *sums;
  double smin,smax;
  int imin,imax;
  double amin;
  int i;
  int retval;

  sums = matrixSums(matrix,m,n);

  smin = sums[0];
  smax = sums[0];

  imin = 0;
  imax = 0;

  for (i = 1; i < n; i++) {
    double sum;

    sum = sums[i];

    if (sum < smin) {
      smin = sum;
      imin = i;
    }

    if (sum > smax) {
      smax = sum;
      imax = i;
    }
  }

  free(sums);

  fprintf(stderr,"   %2d - %8.5lf\n",imin,smin);
  fprintf(stderr,"   %2d - %8.5lf\n",imax,smax);
  fprintf(stderr,"      - %12.5le\n",smax-smin);
  fprintf(stderr,"\n");

  retval = 1;

  if (imin != imax) {
    double a;

    gmatrixRef = matrix;
    gmatrixTmp = matrixIdentity(n);

    gsums = NULL;

    gm = m;
    gn = n;

    gimax = imax;

    imin = -1;

    for (i = 0; i < n; i++) {
      if (i != gimax) {
        gimin = i;

        fprintf(stderr,"   %2d %2d\n",gimax,gimin);
        fprintf(stderr,"   %8.5lf - %8.5lf\n",0.0,diffSum(0.0));
        fprintf(stderr,"   %8.5lf - %8.5lf\n",M_PI/2.0,diffSum(M_PI/2.0));
        fprintf(stderr,"\n");

        a = brentZero(0.0,M_PI/2.0,diffSum,1e-16);

        diffSum(a);

        if (gsmin == 1.0 && gsmax == 1.0) {
          a = M_PI/4.0;
          diffSum(a);
        }

        if (imin == -1 || gsmin < smin) {
          imin = i;
          smin = gsmin;
          amin = a;
        }

        printSum(a);
        fprintf(stderr,"\n");
      }
    }

    matrixFree(gmatrixTmp);
    matrixRotLeft(matrix,amin,imin,imax,n);
  }

  fprintf(stderr,"-- %2d %2d - %8.5lf\n",imax,imin,amin);
  fprintf(stderr,"\n");

  if (imin != imax && amin != 0.0 && amin != M_PI/2) {
    retval = 0;
  }

  return(retval);
}

int main(int argc, char **argv)
{
  int m,n;
  int dim;
  int i;
  double *x;
  double *matrix;
  double *sums;
  double done = 0;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  dim = (m*(2*n-(m+1)))/2;

  x = pointMake(dim);

  for (i = 0; i < dim; i++) {
    x[i] = 0.0;
  }

  matrix = rupertMatrix(x,m,n);

  matrixPrint(stderr,"%8.5lf",matrix,n);
  fprintf(stderr,"\n");

  sums = matrixSums(matrix,m,n);
  pointPrint(stderr,"%8.5lf",sums,n,1);
  fprintf(stderr,"\n");
  pointFree(sums);

  i = 0;
  while (done == 0) {
    fprintf(stderr,"i = %d\n",i);
    done = matchRowSums(matrix,m,n);

    matrixPrint(stderr,"%8.5lf",matrix,n);
    fprintf(stderr,"\n");

    sums = matrixSums(matrix,m,n);
    pointPrint(stderr,"%8.5lf",sums,n,1);
    fprintf(stderr,"\n");
    pointFree(sums);

    i++;
  }

  x = rupertRep(matrix,m,n);

  pointPrint(stderr,"%8.5lf",x,dim,1);
  fwrite(x,sizeof(*x),dim,stdout);

  pointFree(x);
  matrixFree(matrix);

  exit(0);
}
