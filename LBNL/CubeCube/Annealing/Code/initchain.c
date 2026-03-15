#include "define.h"
#include "global.h"

void initchain(CHAIN *chain, int *plength, PARAMETERS *params)
{
  int indim,outdim;
  int length;
  int i;

  indim = params->indim;
  outdim = params->outdim;

  length = 1;

  newchain(chain,length,outdim);

  for (i = 0; i < indim; i++) {
    chain->insum[i] = 1.0;
  }

  chain->maxsum = 1.0;

  *plength = length;
}
