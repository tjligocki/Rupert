#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "rupert.h"

int main(int argc, char **argv)
{
  RUPERT *rup;
  double r;

  ROTATE *arot;

  double da2;
  ROTATE *darot;

  ROTATE *brot;

  double delta;

  int i,j;
  int numpts;

  int m,n;
  int total;

  double smin,rmin;
  int imin;
  ROTATE *bminrot;

  int tries;

  double det;

  int count;

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

#if 0
  a01 = M_PI/4;
  a02 = 0.0;
  a12 = 4 * atan(sqrt(5-2*sqrt(6)));

  printf("Exploring 'lines' in the Rupert spaces\n");
  printf("\n");
#endif

#if 0
  srand48(9453679);
#else
  srand48(time(0));
#endif

  arot = rotate_new(m,n);

  total = arot->rot->n;

  printf("%d, %d -> %d\n",m,n,total);

  for (j = 0; j < total; j++)
  {
    arot->rot->vec[j] = M_PI * (2 * drand48() - 1);
  }

  brot    = rotate_new(m,n);
  darot   = rotate_new(m,n);
  bminrot = rotate_new(m,n);

  numpts = 10000;
  delta = 0.001;

  count = 0;

  while (count < 1000 || delta > 1e-20)
  {
    double minlen2;

    for (j = 0; j < total; j++)
    {
      darot->rot->vec[j] = (2 * drand48() - 1);
    }

    da2 = 0.0;
    for (j = 0; j < total; j++)
    {
      da2 += darot->rot->vec[j]*darot->rot->vec[j];
    }
#if 0
    tries = 1;
    minlen2 = -1.0;
    while (da2 > 1.0 || da2 == 0.0)
    {
      for (j = 0; j < total; j++)
      {
        darot->rot->vec[j] = (2 * drand48() - 1);
      }

      da2 = 0.0;
      for (j = 0; j < total; j++)
      {
        da2 += darot->rot->vec[j]*darot->rot->vec[j];
      }

      if (da2 < minlen2 || minlen2 == -1.0)
      {
        minlen2 = da2;
      }

      tries++;

      if (tries % 100000 == 0)
      {
        printf("  %d, %e\n",tries,sqrt(minlen2));
      }
    }
#endif
    da2 = sqrt(da2);

    for (j = 0; j < total; j++)
    {
      darot->rot->vec[j] /= da2;
    }

    for (i = -numpts/2; i <= numpts/2; i++)
    {
      for (j = 0; j < total; j++)
      {
        brot->rot->vec[j] = arot->rot->vec[j] + i * delta*darot->rot->vec[j];
      }
      
      rup = rupert_new(brot);

      r = rupert_value(rup);
#if 0
      printf("  %9.6f %9.6f\n",i*delta,r);
#endif

      if (i == 0)
      {
        smin = i * delta;
        rmin = r;

        rotate_copy(bminrot,brot);
      }

      if (r < rmin)
      {
        smin = i * delta;
        rmin = r;
        imin = i;

        rotate_copy(bminrot,brot);
      }

      rupert_free(rup);
    }

    if (smin != 0.0)
    {
#if 1
      printf("%19.16f (%d, %e) %19.16f (%19.16f), %d\n",smin,imin,delta,rmin,1.0/rmin,count);

      if (delta < 1.0)
      {
        delta *= 1.1;
      }
#endif
      count = 0;
    }
    else
    {
      count++;
    }

    rotate_copy(arot,bminrot);

    delta *= 0.999;
  }

  for (j = 0; j < total; j++)
  {
    arot->rot->vec[j] -= floor(arot->rot->vec[j]/(2*M_PI)) * (2*M_PI);
    if (arot->rot->vec[j] > M_PI) { arot->rot->vec[j] -= 2*M_PI; }
  }

  printf("%19.16f (%19.16f), %23.16e, ",rmin,1.0/rmin,delta);
  for (j = 0; j < total; j++)
  {
    printf("%20.16f, ",rad2deg(arot->rot->vec[j]));
  }
  printf("\n");

  rotate_free(arot);
  rotate_free(darot);
  rotate_free(brot);
  rotate_free(bminrot);

  return 0;
}
