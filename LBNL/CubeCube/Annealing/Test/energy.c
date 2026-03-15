#include "define.h"
#include "global.h"

REAL energy(CHAIN *newchain, int length, PARAMETERS *params)
{
  int indim,outdim;
  int i,j;
  REAL *scan;
  REAL maxsum,cursum;

  indim = params->indim;
  outdim = params->outdim;

  maxsum = 0.0;

  for (i = 0; i < outdim; i++) {
    scan = newchain->trans[i];

    cursum = 0.0;
    for (j = 0; j < indim; j++) {
      cursum += fabs(*(scan++));
    }

    newchain->insum[i] = cursum;
    
    if (cursum > maxsum) {
      maxsum = cursum;
    }
  }

  newchain->maxsum = maxsum;

  return(-1.0/maxsum);
}
