#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rupert.h"

RUPERT *rupert_new(ROTATE *rot)
{
  RUPERT *rup;
  MATRIX *mat;

  int m,n;
  int total;

  double *scan;

  int i,j;

  rup = (RUPERT *)malloc(sizeof(*rup));

  m     = rot->m;
  n     = rot->n;
  total = rot->total;

  mat = matrix_new(n,n);
  matrix_identity(mat);

  scan = rot->rot + rot->total;
  for (i = m-1; i >= 0; i--)
  {
    for (j = n-1; j > i; j--)
    {
      scan--;
      matrix_rotate(mat,mat,*scan,i,j);
    }
  }

  rup->mat = mat;
  rup->m   = m;
  rup->n   = n;

  return rup;
}

double rupert_value(RUPERT *rup)
{
  double **scan;
  int m,n;

  int i,j;

  double r;

  scan = rup->mat->mat;
  m = rup->m;
  n = rup->n;

  r = 0.0;

  for (i = 0; i < n; i++)
  {
    double sum = 0.0;

    for (j = 0; j < m; j++)
    {
      sum += fabs(scan[i][j]);
    }

    if (sum > r)
    {
      r = sum;
    }
  }

  return r;
}

void rupert_free(RUPERT *rup)
{
  matrix_free(rup->mat);

  free(rup);
}
