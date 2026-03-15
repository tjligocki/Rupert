#include "define.h"
#include "global.h"

static void copydouble(REAL *energy, REAL *newenergy, int length)
{
  if (length > 0) {
    while (length--) {
      *(newenergy++) = *(energy++);
    }
  }
}

void copychain(CHAIN *chain, CHAIN *newchain, int length, int dim)
{
  int i;

  if (length > 0) {
    for (i = 0; i < dim; i++) {
      copydouble(chain->trans[i],newchain->trans[i],dim);
    }

    copydouble(chain->insum,newchain->insum,dim);

    newchain->maxsum = chain->maxsum;
  }
}
