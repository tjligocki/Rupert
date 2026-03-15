#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "rupert.h"

int main(int argc, char **argv)
{
  RUPERT *rup;
  double r;

  double a01,a02,a12;
  double a[3];
  ROTATE arot;

  double da01,da02,da12;
  double da2;
  double da[3];
  ROTATE darot;

  double b[3];
  ROTATE brot;

  double delta;

  int i,n;

  double smin,rmin;
  double bmin[3];

  int tries;

  double det;

  int count;

#if 0
  a01 = M_PI/4;
  a02 = 0.0;
  a12 = 4 * atan(sqrt(5-2*sqrt(6)));

  printf("Exploring 'lines' in the Rupert spaces\n");
  printf("\n");
#endif

#if 0
  srand48(9453679);
#endif
  srand48(time(0));

  a01 = M_PI * (2 * drand48() - 1);
  a02 = M_PI * (2 * drand48() - 1);
  a12 = M_PI * (2 * drand48() - 1);

  a[0] = a12;
  a[1] = a02;
  a[2] = a01;
#if 0
  printf("a01: %9.6f (%9.4f)\n",a[2],rad2deg(a[2]));
  printf("a02: %9.6f (%9.4f)\n",a[1],rad2deg(a[1]));
  printf("a12: %9.6f (%9.4f)\n",a[0],rad2deg(a[0]));
  printf("\n");
#endif
  arot.rot = a;
  arot.m   = 2;
  arot.n   = 3;

  brot.rot = b;
  brot.m   = 2;
  brot.n   = 3;

  n = 10000;
  delta = 0.001;

  count = 0;

  while (count < 1000)
  {
    da01 = (2 * drand48() - 1);
    da02 = (2 * drand48() - 1);
    da12 = (2 * drand48() - 1);

    da2 = da01*da01 + da02*da02 + da12*da12;

    tries = 1;
    while (da2 > 1.0 || da2 == 0.0)
    {
      da01 = (2 * drand48() - 1);
      da02 = (2 * drand48() - 1);
      da12 = (2 * drand48() - 1);

      da2 = da01*da01 + da02*da02 + da12*da12;

      tries++;
    }

    da2 = sqrt(da2);
#if 0
    printf("da2: %9.6f (%9.4f), tries %d\n",da2,rad2deg(da2),tries);
    printf("\n");
#endif
    da[0] = da12 / da2;
    da[1] = da02 / da2;
    da[2] = da01 / da2;
#if 0
    printf("da01: %9.6f (%9.4f)\n",da[2],rad2deg(da[2]));
    printf("da02: %9.6f (%9.4f)\n",da[1],rad2deg(da[1]));
    printf("da12: %9.6f (%9.4f)\n",da[0],rad2deg(da[0]));
    printf("\n");
#endif
    for (i = -n/2; i <= n/2; i++)
    {
      b[0] = a[0] + i * delta*da[0];
      b[1] = a[1] + i * delta*da[1];
      b[2] = a[2] + i * delta*da[2];

      rup = rupert_new(&brot);

      r = rupert_value(rup);
#if 0
      printf("%9.6f %9.6f\n",i*delta,r);
#endif

      if (i == 0)
      {
        smin = i * delta;
        rmin = r;

        bmin[0] = b[0];
        bmin[1] = b[1];
        bmin[2] = b[2];
      }

      if (r < rmin)
      {
        smin = i * delta;
        rmin = r;

        bmin[0] = b[0];
        bmin[1] = b[1];
        bmin[2] = b[2];
      }

      rupert_free(rup);
    }

    if (smin != 0.0)
    {
#if 0
      printf("%19.16f %19.16f (%19.16f), %d\n",smin,rmin,1.0/rmin,count);
#endif
      count = 0;
    }
    else
    {
      count++;
    }

    a[0] = bmin[0];
    a[1] = bmin[1];
    a[2] = bmin[2];

    delta *= 0.99;
  }

  a[0] -= floor(a[0]/(2*M_PI)) * (2*M_PI);
  if (a[0] > M_PI) { a[0] -= 2*M_PI; }
  a[1] -= floor(a[1]/(2*M_PI)) * (2*M_PI);
  if (a[1] > M_PI) { a[1] -= 2*M_PI; }
  a[2] -= floor(a[2]/(2*M_PI)) * (2*M_PI);
  if (a[2] > M_PI) { a[2] -= 2*M_PI; }

  printf("%19.16f (%19.16f), %23.16e, ",rmin,1.0/rmin,delta);
  printf("%20.16f, ",rad2deg(a[2]));
  printf("%20.16f, ",rad2deg(a[1]));
  printf("%20.16f\n",rad2deg(a[0]));

  return 0;
}
