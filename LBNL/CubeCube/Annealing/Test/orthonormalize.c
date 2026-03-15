#include "define.h"
#include "global.h"

void orthonormalize(REAL **trans, int dim)
{
  int i,j,k;
  REAL inner,eps;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < i; j++) {
      inner = 0;

      for (k = 0; k < dim; k++) {
        inner += trans[i][k] * trans[j][k];
      }

      for (k = 0; k < dim; k++) {
        trans[j][k] -= inner*trans[i][k];
      }
    }

    inner = 0;

    for (k = 0; k < dim; k++) {
        inner += trans[i][k] * trans[i][k];
    }

    inner = sqrt(1.0/inner);

    for (k = 0; k < dim; k++) {
      trans[i][k] *= inner;
    }
  }
}
