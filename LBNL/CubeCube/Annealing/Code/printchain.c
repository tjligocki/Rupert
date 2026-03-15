#include "define.h"
#include "global.h"

void printchain(FILE *file, CHAIN *chain, int length, PARAMETERS* params)
{
  int dim;
  int i,j;

  dim = params->outdim;

  for (i = 0; i < dim; i++) {
    fprintf(file,"%3d ",i);

    for (j = 0; j < dim; j++) {
      fprintf(file,"%30.27Lf ",chain->trans[i][j]);
    }

    fprintf(file,"%30.27Lf ",chain->insum[i]);
    fprintf(file,"\n");
  }
  fprintf(file,"\n");

  fprintf(file,"%30.27Lf ",chain->maxsum);
  fprintf(file,"\n");
}
