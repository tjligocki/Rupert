#include "matrix.h"

typedef struct
{
  MATRIX *mat;
  int m,n;
} RUPERT;

RUPERT *rupert_new(ROTATE *rot);

double rupert_value(RUPERT *mat);

void rupert_free(RUPERT *rup);
