#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "tools.h"

static double *gx;
static double *gy;
static double *gvx;

static int gm,gn;
static int gdim;

double rupertLine(double t)
{
  int i;

  for (i = 0; i < gdim; i++) {
    gy[i] = gx[i] + t*gvx[i];
  }

  return rupert(gy,gm,gn);
}

inline checkLines(double *x, double range, int np, int nr, int nt, int m, int n)
{
  double dx;
  double cur;
  double min;
  double tmin;
  double gmin;
  double at,bt,ct;
  int i;
  int it;
  int ip;
  double nvx;

  gm = m;
  gn = n;

  gdim = (m*(2*n-(m+1)))/2;

  gx = x;
  gy = pointMake(gdim);

  gvx = pointMake(gdim);

  gmin = rupert(gx,gm,gn);

  for (it = 0; it < nt; it++) {
    cur = rupert(gx,gm,gn);

    min = cur;
    tmin = 0.0;
    
    fprintf(stderr,"--- %23.16le %23.16le\n",tmin,1.0/min);

    do {
      nvx = 0.0;

      for (i = 0; i < gdim; i++) {
        double cv;

        cv = 2*drand48() - 1;

        gvx[i] = cv;
        nvx += cv*cv;
      }
    /* } while (nvx > 1.0); */
    } while (0.0);

    nvx = sqrt(nvx);

    fprintf(stderr,"\n");
    for (i = 0; i < gdim; i++) {
      gvx[i] = gvx[i] / nvx;
      fprintf(stderr,"* %23.16le\n",gvx[i]);
    }
    fprintf(stderr,"\n");

    at = -0.10;
    bt = -0.01;
    ct =  0.10;

    fprintf(stderr,"=== %23.16le %23.16le\n",at,1.0/rupertLine(at));
    fprintf(stderr,"=== %23.16le %23.16le\n",bt,1.0/rupertLine(bt));
    fprintf(stderr,"=== %23.16le %23.16le\n",ct,1.0/rupertLine(ct));
    fprintf(stderr,"\n");

    min = brentMin(at,bt,ct,rupertLine,1e-16,&tmin);

    if (min < gmin) {
      fprintf(stderr,"+++ %23.16le %23.16le\n",tmin,1.0/min);
      fprintf(stderr,"\n");

      for (i = 0; i < gdim; i++) {
        gx[i] = gx[i] + tmin*gvx[i];
      }

      gmin = min;
    }
  }

  pointFree(gvx);
  pointFree(gy);
}

int main(int argc, char **argv)
{
  int m,n;
  int np;
  int dim;
  int i;
  double *minx;
  FILE *input;
  char name[1000];
  double *minMatrix;
  double minimum;
  double range;
  int nr,nt;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  np = atoi(argv[3]);
  nr = atoi(argv[4]);
  nt = atoi(argv[5]);

  srand48(time(0));

  dim = (m*(2*n-(m+1)))/2;

  range = 2*M_PI;

  minx = pointMake(dim);

  i = fread(minx,sizeof(*minx),dim,stdin);

  checkLines(minx,range/10.0,np,nr,nt,m,n);

  fwrite(minx,sizeof(*minx),dim,stdout);

  pointFree(minx);

  exit(0);
}
