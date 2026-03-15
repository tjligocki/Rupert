#include "f.h"

inline void im(double m[N][N])
{
  int i,j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (i != j) {
        m[i][j] = 0.0;
      } else {
        m[i][j] = 1.0;
      }
    }
  }
}

inline void zm(double m[N][N])
{
  int i,j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      m[i][j] = 0.0;
    }
  }
}

inline void cm(double m1[N][N],
               double m2[N][N])
{
  int i,j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      m1[i][j] = m2[i][j];
    }
  }
}

inline void mm(double m1[N][N],
               double m2[N][N],
               double m3[N][N])
{
  int i,j,k;
  double mt[N][N];

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      double sum = 0.0;
      for (k = 0; k < N; k++) {
        sum += m2[i][k]*m3[k][j];
      }

      mt[i][j] = sum;
    }
  }

  cm(m1,mt);
}

void pm(double m[N][N])
{
  int i,j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      fprintf(stderr,"%19.12le ",m[i][j]);
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
}

bool f_debug = false;

double f(double x[D])
{
  int d1,d2,k;
  double m[N][N];
  double mr[N][N];
  double max;

  im(m);

  k = 0;
  for (d1 = 0; d1 < M; d1++) {
    for (d2 = d1+1; d2 < N; d2++) {
      double cx,sx;

      cx = cos(x[k]);
      sx = sin(x[k]);

      im(mr);

      mr[d1][d1] = cx; mr[d1][d2] = -sx;
      mr[d2][d1] = sx; mr[d2][d2] =  cx;

      mm(m,mr,m);

      k++;
    }
  }

  max = 0.0;

  for (d1 = 0; d1 < N; d1++) {
    double cur = 0.0;

    for (d2 = 0; d2 < M; d2++) {
      cur += fabs(m[d1][d2]);
    }

    if (cur > max) {
      max = cur;
    }
  }

  return max;
}

void df(double d[D], double x[D])
{
  int d1,d2,i,k;
  double m[N][N];
  double dm[D][N][N];
  double mr[N][N];
  double dmr[N][N];
  double max;

  im(m);
  for (i = 0; i < D; i++) {
    im(dm[i]);
  }

  k = 0;
  for (d1 = 0; d1 < M; d1++) {
    for (d2 = d1+1; d2 < N; d2++) {
      double cx,sx;

      cx = cos(x[k]);
      sx = sin(x[k]);

      im(mr);

      mr[d1][d1] = cx; mr[d1][d2] = -sx;
      mr[d2][d1] = sx; mr[d2][d2] =  cx;

      zm(dmr);

      dmr[d1][d1] = -sx; dmr[d1][d2] = -cx;
      dmr[d2][d1] =  cx; dmr[d2][d2] = -sx;

      mm(m,mr,m);

      for (i = 0; i < D; i++) {
        if (i == k) {
          mm(dm[i],dmr,dm[i]);
        } else {
          mm(dm[i],mr,dm[i]);
        }
      }

      k++;
    }
  }

#if 0
  for (d1 = 0; d1 < N; d1++) {
    for (d2 = 0; d2 < N; d2++) {
      fprintf(stderr,"%10.3e ",m[d1][d2]);
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  for (i = 0; i < D; i++) {
    for (d1 = 0; d1 < N; d1++) {
      for (d2 = 0; d2 < N; d2++) {
        fprintf(stderr,"%10.3e ",dm[i][d1][d2]);
      }
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
  }
#endif

  max = 0.0;

  for (d1 = 0; d1 < N; d1++) {
    double cur = 0.0;
    double dcur[D];

    for (i = 0; i < D; i++) {
      dcur[i] = 0.0;
    }

    for (d2 = 0; d2 < M; d2++) {
      if (m[d1][d2] > 0) {
        cur += m[d1][d2];
        for (i = 0; i < D; i++) {
          dcur[i] += dm[i][d1][d2];
        }
      } else {
        cur += -m[d1][d2];
        for (i = 0; i < D; i++) {
          dcur[i] += -dm[i][d1][d2];
        }
      }
    }

    if (cur > max) {
      max = cur;
      for (i = 0; i < D; i++) {
        d[i] = dcur[i];
      }
    }
  }

  if (f_debug) {
    fprintf(stderr,"df:\n");
    for (i = 0; i < D; i++) {
      fprintf(stderr,"%19.12le ",d[i]);
    }
    fprintf(stderr,"\n\n");
  }
}
