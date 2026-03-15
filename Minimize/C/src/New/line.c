#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "rupert.h"

int main(int argc, char **argv)
{
  int m,n,total;

  RUPERT *rup;
  double r;

  ROTATE *arot;

  ROTATE *darot;
  double da2;

  ROTATE *brot;

  double delta;

  int i,j;

  double smin,rmin;
  ROTATE *bminrot;

  int steps;
  int tries;

  m = 2;
  if (argc > 1)
  {
    m = atoi(argv[1]);
  }

  n = 3;
  if (argc > 2)
  {
    n = atoi(argv[2]);
  }

  arot    = rotate_new(m,n);
  darot   = rotate_new(m,n);
  brot    = rotate_new(m,n);
  bminrot = rotate_new(m,n);

  total = arot->total;

  printf("%d, %d -> %d\n",m,n,total);

#if 0
  a01 = M_PI/4;
  a02 = 0.0;
  a12 = 4 * atan(sqrt(5-2*sqrt(6)));

  printf("Exploring 'lines' in the Rupert spaces\n");
  printf("\n");
#endif

#if 1
  srand48(9453679);
#endif
  srand48(time(0));

  for (j = 0; j < total; j++)
  {
    arot->rot[j] = M_PI * (2 * drand48() - 1);
  }
#if 1
  for (j = 0; j < total; j++)
  {
    printf("%9.6lf (%9.4lf), ",arot->rot[j],rad2deg(arot->rot[j]));
  }
  printf("\n");
  printf("\n");
#endif

  steps = 10000;
  delta = 0.001;

  rmin  = 2.0;

  while (delta > 1e-20)
  {
    for (j = 0; j < total; j++)
    {
      darot->rot[j] = (2 * drand48() - 1);
    }

    da2 = 0.0;
    for (j = 0; j < total; j++)
    {
      da2 += darot->rot[j]*darot->rot[j];
    }

    tries = 1;
    while (da2 > 1.0 || da2 == 0.0)
    {
      for (j = 0; j < total; j++)
      {
        darot->rot[j] = (2 * drand48() - 1);
      }

      da2 = 0.0;
      for (j = 0; j < total; j++)
      {
        da2 += darot->rot[j]*darot->rot[j];
      }

      tries++;
    }

    da2 = sqrt(da2);
#if 0
    printf("da2: %9.6lf (%9.4lf), tries %d\n",da2,rad2deg(da2),tries);
    printf("\n");
#endif
    for (j = 0; j < total; j++)
    {
      darot->rot[j] = darot->rot[j] / da2;
    }
#if 0
    for (j = 0; j < total; j++)
    {
      printf("%9.6lf (%9.4lf), ",darot->rot[j],rad2deg(darot->rot[j]));
    }
    printf("\n");
    printf("\n");
#endif
    for (i = -steps/2; i <= steps/2; i++)
    {
      for (j = 0; j < total; j++)
      {
        brot->rot[j] = arot->rot[j] + i * delta * darot->rot[j];
      }

      rup = rupert_new(brot);

      r = rupert_value(rup);
#if 0
      printf("%9.6lf %9.6lf\n",i*delta,r,rmin);
#endif

      if (r < rmin)
      {
#if 0
      printf("%d %19.16lf %19.16lf, %le\n",i,r,rmin,delta);
#endif
        smin = i * delta;
        rmin = r;

        for (j = 0; j < total; j++)
        {
          bminrot->rot[j] = brot->rot[j];
        }
      }

      rupert_free(rup);
    }

    if (smin != 0.0)
    {
#if 1
      printf("%19.16lf %19.16lf (%19.16lf), %le, ",smin,rmin,1.0/rmin,delta);
      for (j = 0; j < total; j++)
      {
        printf("%19.14lf, ",rad2deg(bminrot->rot[j]));
      }
      printf("\n");
#endif
      smin = 0.0;

      delta *= 0.9999;
    }
#if 0
    printf("%le\n",delta);
#endif
    for (j = 0; j < total; j++)
    {
      arot->rot[j] = bminrot->rot[j];
    }

    delta *= 0.9;
#if 0
    printf("%le\n",delta);
#endif
  }

  for (j = 0; j < total; j++)
  {
    arot->rot[j] -= floor(arot->rot[j]/(2*M_PI)) * (2*M_PI);
    if (arot->rot[j] > M_PI) { arot->rot[j] -= 2*M_PI; }
  }

  printf("%19.16lf (%19.16lf), %23.16e, ",rmin,1.0/rmin,delta);
  for (j = 0; j < total; j++)
  {
    printf("%9.4lf, ",rad2deg(arot->rot[j]));
  }
  printf("\n");
  printf("\n");

  rotate_free(arot);
  rotate_free(darot);
  rotate_free(brot);
  rotate_free(bminrot);

  return 0;
}
