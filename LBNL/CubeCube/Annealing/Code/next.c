#include "define.h"
#include "global.h"

static REAL move(CHAIN *, CHAIN *, int, int, int, REAL, PARAMETERS *);

int next(CHAIN *chain, CHAIN *newchain, int length, REAL sd, PARAMETERS *params, REAL *ener)
{
  int indim,outdim;
  int axis1,axis2;
  int i,j;
  int count;
  REAL cur;

  indim = params->indim;
  outdim = params->outdim;

  if (params->distribution == UNIFORM) {
    sd = sqrt(3.0) * sd;
  }

  count = 0;

  if (outdim == 1) {
    copychain(chain,newchain,length,outdim);

    cur = 1.0;

    count++;
  } else
  if (params->movement == LOCAL) {
    axis1 = drand48() * indim;

    axis2 = drand48() * (outdim - 1);
    if (axis2 >= axis1) {
      axis2++;
    }

    cur = move(chain,newchain,length,axis1,axis2,sd,params);

    count++;
  }

  *ener = -1.0/cur;

  return(count);
}

#define RANDOM  (2*drand48() - 1)

static REAL normal()
{
  static int flag = 0;
  static REAL r1,r2;
  REAL r;

  if (flag == 0) {
    do {
      r1 = RANDOM;
      r2 = RANDOM;

      r = r1*r1 + r2*r2;
    } while (r > 1);
    
    r = sqrt(-2*log(r)/r);

    r1 = r*r1;
    r2 = r*r2;

    flag = 1;

    return(r1);
  } else {
    flag = 0;

    return(r2);
  }
}

static REAL move(CHAIN *chain, CHAIN *newchain, int length, int axis1, int axis2, REAL scale, PARAMETERS *params)
{
  int indim,outdim;
  REAL phi;
  REAL cphi,sphi;
  int i,j;
  REAL *s1,*s2,*s3,*s4;
  REAL sum1,sum2;
  REAL maxsum;
  REAL **tr1,**tr2;

  indim = params->indim;
  outdim = params->outdim;

  if (params->distribution == UNIFORM) {
    phi = scale * RANDOM;
  } else {
    phi = scale * normal();
  }

  cphi = cos(phi);
  sphi = sin(phi);

  tr1 = chain->trans;
  tr2 = newchain->trans;

  memcpy(newchain->insum,chain->insum,outdim*sizeof(REAL));

  for (i = 0; i < outdim; i++) {
    s1 = *(tr1++);
    s2 = *(tr2++);

    memcpy(s2,s1,outdim*sizeof(REAL));
  }

  s1 = chain->trans[axis1];
  s2 = chain->trans[axis2];

  s3 = newchain->trans[axis1];
  s4 = newchain->trans[axis2];

  sum1 = 0.0;
  sum2 = 0.0;

  for (i = 0; i < outdim; i++) {
    REAL m1,m2,m3,m4;

    m1 = *(s1++);
    m2 = *(s2++);

    *(s3++) = m3 = cphi * m1 - sphi * m2;
    *(s4++) = m4 = sphi * m1 + cphi * m2;

    if (i < indim) {
      sum1 += fabs(m3);
      sum2 += fabs(m4);
    }
  }

  newchain->insum[axis1] = sum1;
  newchain->insum[axis2] = sum2;

  maxsum = 0.0;

  for (i = 0; i < outdim; i++) {
    REAL cursum;

    cursum = newchain->insum[i];

    if (cursum > maxsum) maxsum = cursum;
  }

  newchain->maxsum = maxsum;

  return(maxsum);
}
