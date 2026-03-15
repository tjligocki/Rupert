#include <cstdio>
#include <cmath>
#include <ctime>

#include "angle.h"
#include "rotate.h"
#include "rupert.h"

int main(int argc, char **argv)
{
  int m = 2;
  if (argc > 1)
  {
    m = atoi(argv[1]);
  }

  int n = 3;
  if (argc > 2)
  {
    n = atoi(argv[2]);
  }

#if 0
  srand48(9453679);
#else
  srand48(time(0));
#endif

  ROTATE arot(m,n);

  int total = arot.m_total;

  printf("%d, %d -> %d\n",m,n,total);

  for (int i = 0; i < total; i++)
  {
    arot.m_rot[i] = M_PI * (2 * drand48() - 1);
  }

  ROTATE brot(m,n);
  ROTATE darot(m,n);

  int numpts = 10000;
  double delta = 0.001;

  double smin;
  double rmin;
  int imin;
  ROTATE bminrot(m,n);

  int count = 0;

  while (count < 1000 || delta > 1e-20)
  {
    double minlen2;

    for (int i = 0; i < total; i++)
    {
      darot.m_rot[i] = (2 * drand48() - 1);
    }

    darot.norm();

    for (int i = -numpts/2; i <= numpts/2; i++)
    {
      for (int j = 0; j < total; j++)
      {
        brot.m_rot[j] = arot.m_rot[j] + i * delta*darot.m_rot[j];
      }
      
      RUPERT rup(brot);

      double r = rup.value();
#if 0
      printf("  %9.6f %9.6f\n",i*delta,r);
#endif

      if (i == 0)
      {
        smin = i * delta;
        rmin = r;

        bminrot.copy(brot);
      }

      if (r < rmin)
      {
        smin = i * delta;
        rmin = r;
        imin = i;

        bminrot.copy(brot);
      }
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

    arot.copy(bminrot);

    delta *= 0.999;
  }

  for (int j = 0; j < total; j++)
  {
    arot.m_rot[j] -= floor(arot.m_rot[j]/(2*M_PI)) * (2*M_PI);
    if (arot.m_rot[j] > M_PI) { arot.m_rot[j] -= 2*M_PI; }
  }

  printf("%19.16f (%19.16f), %23.16e, ",rmin,1.0/rmin,delta);
  for (int j = 0; j < total; j++)
  {
    printf("%20.16f, ",rad2deg(arot.m_rot[j]));
  }
  printf("\n");

  return 0;
}
