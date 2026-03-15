#include "define.h"
#include "global.h"

void newchain(CHAIN *chain, int length, int dim)
{
  REAL **cur_t,*cur_r;
  REAL *cur_is;
  int i,j;

  cur_t = (REAL **)malloc((unsigned)(dim*sizeof(*cur_t)));

  if (cur_t == NULL) {
    fprintf(stderr,"Not enough memory for chain transform\n");
    exit(1);
  }

  cur_is = (REAL *)malloc((unsigned)(dim*sizeof(*cur_is)));

  for (i = 0; i < dim; i++) {
    cur_r = (REAL *)malloc((unsigned)(dim*sizeof(*cur_r)));

    if (cur_r == NULL) {
      fprintf(stderr,"Not enough memory for chain transform\n");
      exit(1);
    }

    for (j = 0; j < dim; j++) {
      if (j != i) {
        cur_r[j] = 0.0;
      } else {
        cur_r[j] = 1.0;
      }
    }

    cur_t[i]  = cur_r;
    cur_is[i] = 0.0;
  }

  chain->trans = cur_t;
  chain->insum = cur_is;
  chain->maxsum = 0.0;
}
