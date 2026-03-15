#include "f1d.h"

static double f1d_x0[D];
static double f1d_v0[D];

void f1d_init(double x0[D],
              double v0[D])
{
  int i;
  double lv0 = 0.0;

  for (i = 0; i < D; i++) {
    lv0 += v0[i]*v0[i];
  }
  lv0 = sqrt(lv0);

  for (i = 0; i < D; i++) {
    f1d_x0[i] = x0[i];
    f1d_v0[i] = v0[i]/lv0;
  }
}

void f1d_to_3d(double x, double x3d[D])
{
  int i;

  for (i = 0; i < D; i++) {
    x3d[i] = f1d_x0[i] + x*f1d_v0[i];
  }
}

double f1d(double x)
{
  double x3d[D];

  f1d_to_3d(x,x3d);

  return f(x3d);
}

double df1d(double x)
{
  int i;
  double x3d[D];
  double df3d[D];
  double m;

  f1d_to_3d(x,x3d);

  df(df3d,x3d);

  m = 0.0;
  for (i = 0; i < D; i++) {
    m += df3d[i] * f1d_v0[i];
  }

  return m;
}
