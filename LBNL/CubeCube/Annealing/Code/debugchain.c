#include "define.h"
#include "global.h"

void debugchain(FILE *file, CHAIN *chain, int length, int dim)
{
  int i,j;

  fprintf(file,"Matrix:\n");

  for (i = 0; i < dim; i++) {
    fprintf(file,"  ");
    for (j = 0; j < dim; j++) {
      fprintf(file,"%30.27Lf ",chain->trans[i][j]);
    }
    fprintf(file,"- %30.27Lf\n",chain->insum[i]);
  }

  fprintf(file,"Max sum: %30.27Le\n",chain->maxsum);
  fprintf(file,"\n");
}
