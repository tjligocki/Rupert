#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"

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

inline void matrixCopy(double *matrixOut, double *matrixIn, int n)
{
  double *scanIn,*scanOut;
  int i,j;

  scanIn  = matrixIn;
  scanOut = matrixOut;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      *scanOut = *scanIn;
      scanIn++;
      scanOut++;
    }
  }
}

inline void matrixRotLeft(double *matrix, double angle, int i1, int i2, int n)
{
  double cx,sx;
  int k;

  cx = cos(angle);
  sx = sin(angle);

  for (k = 0; k < n; k++) {
    int ii1,ii2;
    double ti1,ti2;

    ii1 = i1*n + k;
    ii2 = i2*n + k;

    ti1 = matrix[ii1];
    ti2 = matrix[ii2];

    matrix[ii1] = cx*ti1 - sx*ti2;
    matrix[ii2] = sx*ti1 + cx*ti2;
  }
}

inline void matrixRotRight(double *matrix, double angle, int i1, int i2, int n)
{
  double cx,sx;
  int k;

  cx = cos(angle);
  sx = sin(angle);

  for (k = 0; k < n; k++) {
    int ii1,ii2;
    double ti1,ti2;

    ii1 = i1 + k*n;
    ii2 = i2 + k*n;

    ti1 = matrix[ii1];
    ti2 = matrix[ii2];

    matrix[ii1] =  cx*ti1 + sx*ti2;
    matrix[ii2] = -sx*ti1 + cx*ti2;
  }
}

inline double *matrixSums(double *matrix, int m, int n)
{
  double *sums;
  int i;

  sums = pointMake(n);

  for (i = 0; i < n; i++) {
    double *scan;
    double curRow = 0.0;
    int j;

    scan = matrix + i*n;

    for (j = 0; j < m; j++) {
      curRow += fabs(*scan);
      scan++;
    }

    sums[i] = curRow;
  }

  return(sums);
}

inline void matrixPrint(FILE *fd, char *format, double *matrix, int n)
{
  double *scan;
  int i,j;

  scan = matrix;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      fprintf(fd,format,*scan);
      fprintf(fd," ");
      scan++;
    }
    fprintf(fd,"\n");
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
      matrixRotRight(gRot,x[index],i,j,n);
      index++;
    }
  }

  return gRot;
}

inline double *rupertRep(double *matrix, int m, int n)
{
  int dim;
  double *x;
  int i,j;
  int index;
  double *matrixTmp;

  matrixTmp = matrixIdentity(n);
  matrixCopy(matrixTmp,matrix,n);

  dim = (m*(2*n-(m+1)))/2;

  x = pointMake(dim);

  for (i = 0; i < dim; i++) {
    x[i] = 0.0;
  }

  index = 0;

  for (i = 0; i < m; i++) {
    for (j = i+1; j < n; j++) {
      double e1,e2;
      double a;

      e1 = matrixTmp[i*n + i];
      e2 = matrixTmp[j*n + i];
      
      a = atan2(e2,e1);

      matrixRotLeft(matrixTmp,-a,i,j,n);

      x[index] = a;
      index++;
    }
  }

  matrixFree(matrixTmp);

  return(x);
}

inline double *rupertSums(double *x, int m, int n)
{
  double *gRot;
  double *sums;
  int index;
  int i,j;

  gRot = rupertMatrix(x,m,n);

  sums = matrixSums(gRot,m,n);

  matrixFree(gRot);

  return(sums);
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

inline void indexPrint(FILE *fd, char *format, int *index, int dim, int final)
{
  int i;

  for (i = 0; i < dim; i++) {
    fprintf(fd,format,index[i]);
    fprintf(fd," ");
  }

  if (final) {
    fprintf(fd,"\n");
  }
}

void indexFree(int *index)
{
  free(index);
}

inline double *pointMake(int dim)
{
  double *point,*scan;
  int i;

  point = (double *)malloc(dim*sizeof(*point));

  scan = point;
  for (i = 0; i < dim; i++)
  {
    *scan = 0.0;
    scan++;
  }

  return point;
}

inline void pointPrint(FILE *fd, char *format, double *point, int dim, int final)
{
  int i;

  for (i = 0; i < dim; i++) {
    fprintf(fd,format,point[i]);
    fprintf(fd," ");
  }

  if (final) {
    fprintf(fd,"\n");
  }
}

inline void pointFree(double *point)
{
  free(point);
}

inline double brentZero(double x1, double x2, double (*f)(double), double tol) {
  int ITMAX = 1000;
  double EPS = 1e-16;
  double xzero;
  double a,b;
  double fa,fb;

  a = x1;
  b = x2;

  fa = f(a);
  fb = f(b);

  if ((fa > 0 && fb > 0) || (fa < 0 && fb < 0)) {
    fprintf(stderr,"The root must be bracketed to find a root\n");
  } else {
    int i;
    double c,d,e;
    double fc;
    double tol1;
    double xm;
    double p,q,r,s;

    c  = b;
    fc = fb;

    for (i = 0; i < ITMAX; i++) {
      if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
        c  = a;
        fc = fa;

        d = b-a;
        e = d;
      }

      if (fabs(fc) < fabs(fb)) {
        a = b;
        b = c;
        c = a;

        fa = fb;
        fb = fc;
        fc = fa;
      }

      tol1 = 2.0*EPS*fabs(b) + 0.5*tol;

      xm = 0.5*(c-b);

      if ((fabs(xm) <= tol1) || (fb == 0)) {
        break;
      }

      if (fabs(e) <= tol1 && fabs(fa) > fabs(fb)) {
        s = fb/fa;

        if (a == c) {
          p = 2.0*xm*s;
          q = 1.0 - s;
        } else {
          q = fa/fc;
          r = fb/fc;

          p = s*(2.0*xm*q*(q-r) - (b-a)*(r-1.0));
          q = (q-1.0)*(r-1.0)*(s-1.0);
        }

        if (p > 0) {
          q = -q;
        }

        p = fabs(p);

        if (2.0*p < fmin(3.0*xm*q - fabs(tol1*q),fabs(e*q))) {
          e = d;
          d = p/q;
        } else {
          d = xm;
          e = d;
        }
      } else {
        d = xm;
        e = d;
      }

      a  = b;
      fa = fb;

      if (fabs(d) > tol1) {
        b = b + d;
      } else {
        b = b + copysign(tol1,xm);
      }

      fb = f(b);
    }

    if (i == ITMAX) {
      fprintf(stderr,"Brent's algorithm didn't converge\n");
    }
  }

  return(b);
}

inline double brentMin(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin) {
  int ITMAX = 1000;
  double CGOLD = 1.5 - 0.5*sqrt(5.0);
  double ZEPS = 1e-10;

  double a,b;
  double u,v,w,x;
  double d,e;

  double fx,fu,fv,fw;

  int i;

  a = fmin(ax,cx);
  b = fmax(ax,cx);

  v = bx;
  w = v;
  x = v;
  e = 0.0;

  fx = f(x);
  fv = fx;
  fw = fx;

  for (i = 0; i < ITMAX; i++) {
    double xm;
    double tol1,tol2;

    xm = 0.5 * (a+b);

    tol1 = tol*fabs(x) + ZEPS;
    tol2 = 2.0 * tol1;

    if (fabs(x-xm) <= (tol2 * 0.5*(b-a))) {
      break;
    }

    if (fabs(e) > tol1) {
      double r,p,q;
      double etemp;

      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q - (x-w)*r;

      q = 2.0 * (q-r);

      if (q > 0.0) {
        p = -p;
      }

      q = fabs(q);

      etemp = e;
      e = d;

      if ((fabs(p) >= fabs(0.5*q*etemp)) ||
          (p <= q*(a-x)) ||
          (p >= q*(b-x))
         ) {
        if (x >= xm) {
          e = a-x;
        } else {
          e = b-x;
        }
        d = CGOLD*e;
      } else {
        d = p/q;
        u = x+d;
        if (((u-a) < tol2) || ((b-u) < tol2)) {
          d = copysign(tol1,xm-x);
        }
      }
    } else {
      if (x >= xm) {
        e = a-x;
      } else {
        e = b-x;
      }
      d = CGOLD*e;
    }

    if (fabs(d) >= tol1) {
      u = x+d;
    } else {
      u = x+copysign(tol1,d);
    }

    fu = f(u);

    if (fu <= fx) {
      if (u >= x) {
        a = x;
      } else {
        b = x;
      }

      v = w;
      fv = fw;

      w = x;
      fw = fx;

      x = u;
      fx = fu;
    } else {
      if (u < x) {
        a = u;
      } else {
        b = u;
      }

      if ((fu <= fw) || (w == x)) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if ((fu >= fv) || (v == x) || (v == w)) {
        v = u;
        fv = fu;
      }
    }
  }

  if (i == ITMAX) {
    fprintf(stderr,"Brent's algorithm didn't converge\n");
  }

  *xmin = x;
  return(fx);
}
