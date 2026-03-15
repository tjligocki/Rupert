#include "f.h"
#include "f1d.h"
#include "brent.h"

static double g1d_b;
static double g1d_m;

void g1d_init(double x)
{
  g1d_b = f1d(x);
  g1d_m = df1d(x);
}

double g1d(double x)
{
  return(f1d(x) - (g1d_m*x + g1d_b));
}

void findInter(double xm[D], double vm[D], double xm_new[D])
{
  int j;
  int n1d;
  double x1d;
  double x1d_0;
  double dx1d;
  double df1d_0;
  double a,b,c;
  double ga,gb;

  f1d_init(xm,vm);

  x1d_0 = -0.001;

  n1d  = 10000;
  dx1d = 0.002 / n1d;

  for (j = 0; j <= n1d; j++) {
    x1d = x1d_0 + j*dx1d;
    fprintf(stdout,"%19.12le %19.12le %19.12le\n",x1d,f1d(x1d),df1d(x1d));
  }
  fprintf(stdout,"\n");

  g1d_init(0.0);

  df1d_0 = df1d(0.0);

  if (df1d_0 > 0.0) {
    b = -0.1;
  } else {
    b =  0.1;
  }

  gb = g1d(b);

  a = b;
  ga = gb;

  while (ga*gb > 0.0) {
    a /= 1.01;
    ga = g1d(a);
  }

  b = 1.01*a;
  gb = g1d(b);

  fprintf(stderr,"a: %19.12le %19.12le\n",a,ga);
  fprintf(stderr,"b: %19.12le %19.12le\n",b,gb);

  c = brentRootFinder(a,b,g1d);

  fprintf(stdout,"%19.12le %19.12le %19.12le\n",c,f1d(c),g1d(c));
  fprintf(stdout,"\n");

  fprintf(stderr,"c: %19.12le %19.12le %19.12le\n",c,f1d(c),g1d(c));
  fprintf(stderr,"\n");

  f1d_to_3d(c,xm_new);
}

int main()
{
  int i,n;
  double x[D];
  double sdx,dx[D];
  double mx;
  double fm;
  int im;
  double xm[D];
  double p;
  bool found;
  bool first,firstm;

  mx = 100000.0;
  n = 200000001;

  sdx = mx/(n-1);
  p = 2.0;

  for (i = 0; i < D; i++) {
    x[i] = 0.0;
    dx[i] = sdx * p;
    p = sqrt(p);
#if 0
    fprintf(stderr,"%19.12le ",dx[i]);
#endif
  }
#if 0
  fprintf(stderr,"\n");
  fprintf(stderr,"\n");
#endif

  fm = 100.0;
  found = true;
  first = true;

  for (i = 0; i < n; i++) {
    int j;
    double fi;

    fi = f(x);

    if (i % 10000000 == 0) {
      fprintf(stderr,"+++ %d %19.12le\n",i,fi);
    }

    if (fi < fm) {
      fm = fi;
      im = i;
      for (j = 0; j < D; j++) {
        xm[j] = x[j];
      }
      firstm = first;
      found = true;
    }
    else
    {
      if (found) {
        int k;
        double dxm[D];
        double xm_new[D];
        double xm_newer[D];

        fprintf(stderr,"%12d           %19.12le\n",i,fm);

        df(dxm,xm);

        fprintf(stderr,"xm:  ");
        for (k = 0; k < D; k++) {
          fprintf(stderr,"%19.12le ",xm[k]);
        }
        fprintf(stderr,"\n");

        fprintf(stderr,"dxm: ");
        for (k = 0; k < D; k++) {
          fprintf(stderr,"%19.12le ",dxm[k]);
        }
        fprintf(stderr,"\n");

        findInter(xm,dxm,xm_new);
#if 1
        df(dxm,xm_new);

        fprintf(stderr,"xm:  ");
        for (k = 0; k < D; k++) {
          fprintf(stderr,"%19.12le ",xm_new[k]);
        }
        fprintf(stderr,"\n");

        fprintf(stderr,"dxm: ");
        for (k = 0; k < D; k++) {
          fprintf(stderr,"%19.12le ",dxm[k]);
        }
        fprintf(stderr,"\n");

        findInter(xm_new,dxm,xm_newer);
#endif

        found = false;
      }
    }

    for (j = 0; j < D; j++) {
      double f;
      f = (x[j] + dx[j])/(2*M_PI);

      x[j] = 2*M_PI * (f - floor(f));
    }

    first = false;
  }
}
