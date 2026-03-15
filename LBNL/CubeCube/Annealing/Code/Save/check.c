#include "define.h"
#include "global.h"

static REAL row_dot(REAL *, REAL *, int);

void check(FILE *file, CHAIN *curchain, int length, PARAMETERS* params)
{
  REAL err,max;
  int dim;
  int i,j;

  dim = params->outdim;

  max = 0.0;

  for (i = 0; i < length; i++) {
    REAL length;

    length = sqrt(row_dot(curchain->trans[i],curchain->trans[i],dim));

    err = fabs(1.0 - length);

    if (err > max) max = err;
  }

  fprintf(file,"Error: length: %34.27Le, ",max);

  max = 0.0;

  for (i = 0; i < length; i++) {
    REAL length1,length2;
    REAL angle;

    for (j = 0; j < length-1; j++) {
      length1 = sqrt(row_dot(curchain->trans[i],curchain->trans[i],dim));
      length2 = sqrt(row_dot(curchain->trans[j],curchain->trans[j],dim));

      angle = acos(row_dot(curchain->trans[i],curchain->trans[j],dim)
                  / (length1 * length2));

      err = fabs(angle);

      if (err > max) max = err;
    }
  }

  fprintf(file,"angle: %34.27Le\n",max);
}

static REAL row_dot(REAL *row1, REAL *row2, int dim)
{
  REAL dot;
  int i;

  dot = 0.0;
  for (i = 0; i < dim; i++) {
    dot += *(row1++) * *(row2++); 
  }

  return(dot);
}
