#include "f.h"

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
  double known;

#if 1
  known = 1.0/sqrt(9.0/8.0);
#else
  known = 1.0/sqrt(3.0);
#endif

  mx = 1000000.0;
  n = 1000000000;

  srand48(12345678);

  for (i = 0; i < D; i++) {
    x[i] = 2*M_PI*drand48();
  }

  fm = 100.0;

  for (i = 0; i < n; i++) {
    int j;
    double fi;

    fi = f(x);

    if (fi < fm) {
      fm = fi;
      im = i;
      for (j = 0; j < D; j++) {
        xm[j] = x[j];
      }
      fprintf(stderr,"%12d ",im);
      for (j = 0; j < D; j++) {
        fprintf(stderr,"%.3f ",xm[j]);
      }
      fprintf(stderr,"%.12f (%.12f)\n",fm,known);
    }

    for (j = 0; j < D; j++) {
      x[j] = 2*M_PI*drand48();
    }
  }
}
