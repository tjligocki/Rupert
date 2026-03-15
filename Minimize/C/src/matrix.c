#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"

double rad2deg(double rad)
{
  return rad / M_PI * 180.0;
}

double deg2rad(double deg)
{
  return deg / 180.0 * M_PI;
}

MATRIX *matrix_new(int n, int m)
{
  int i;
  MATRIX *temp;
  double **scan;

  temp = (MATRIX *)malloc(sizeof(MATRIX));

  temp->mat = (double **)malloc(n*sizeof(double *));
  temp->n   = n;
  temp->m   = m;

  scan = temp->mat;
  for (i = 0; i < m; i++)
  {
    scan[i] = (double *)malloc(m*sizeof(double));
  }

  return temp;
}

void matrix_zero(MATRIX *mat1)
{
  double **scan1;
  int n1,m1;

  int i,j;

  scan1 = mat1->mat;
  n1    = mat1->n;
  m1    = mat1->m;

  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < m1; j++)
    {
      scan1[i][j] = 0.0;
    }
  }
}

void matrix_identity(MATRIX *mat1)
{
  double **scan1;
  int n1,m1;

  int i,j;

  scan1 = mat1->mat;
  n1    = mat1->n;
  m1    = mat1->m;

  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < m1; j++)
    {
      if (i != j)
      {
        scan1[i][j] = 0.0;
      } else {
        scan1[i][j] = 1.0;
      }
    }
  }
}

void matrix_add(MATRIX *mat1, MATRIX *mat2, MATRIX *mat3)
{
  double **scan1;
  int n1,m1;

  double **scan2;
  int n2,m2;

  double **scan3;
  int n3,m3;

  int i,j;

  scan1 = mat1->mat;
  n1    = mat1->n;
  m1    = mat1->m;

  scan2 = mat2->mat;
  n2    = mat2->n;
  m2    = mat2->m;

  scan3 = mat3->mat;
  n3    = mat3->n;
  m3    = mat3->m;

  if (n1 != n2 || n2 != n3 || m1 != m2 || m2 != m3)
  {
    fprintf(stderr,"matrix_add: matrix dimensions (%d x %d, %d x %d, %d x %d) incompatible\n",n1,m1,n2,m2,n3,m3);
    exit(1);
  }

  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < m1; j++)
    {
      scan1[i][j] = scan2[i][j] + scan3[i][j];
    }
  }
}

void matrix_scalar(MATRIX *mat1, MATRIX *mat2, double a)
{
  double **scan1;
  int n1,m1;

  double **scan2;
  int n2,m2;

  int i,j;

  scan1 = mat1->mat;
  n1    = mat1->n;
  m1    = mat1->m;

  scan2 = mat2->mat;
  n2    = mat2->n;
  m2    = mat2->m;

  if (n1 != n2 || m1 != m2)
  {
    fprintf(stderr,"matrix_scalar: matrix dimensions (%d x %d, %d x %d) incompatible\n",n1,m1,n2,m2);
    exit(1);
  }

  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < m1; j++)
    {
      scan1[i][j] = a * scan2[i][j];
    }
  }
}

void matrix_rotate(MATRIX *mat1, MATRIX *mat2, double a, int i1, int i2)
{
  double **scan1;
  int n1,m1;

  double **scan2;
  int n2,m2;

  double e1,e2;

  double sa,ca;

  int i,j;

  scan1 = mat1->mat;
  n1    = mat1->n;
  m1    = mat1->m;

  scan2 = mat2->mat;
  n2    = mat2->n;
  m2    = mat2->m;

  sa = sin(a);
  ca = cos(a);

  if (n1 != n2 || m1 != m2 || n1 != m1)
  {
    fprintf(stderr,"matrix_rotate: matrix dimensions (%d x %d, %d x %d) incompatible\n",n1,m1,n2,m2);
    exit(1);
  }

  for (i = 0; i < n1; i++)
  {
    e1 = ca * scan2[i1][i] - sa * scan2[i2][i];
    e2 = sa * scan2[i1][i] + ca * scan2[i2][i];

    scan1[i1][i] = e1;
    scan1[i2][i] = e2;
  }
}

void matrix_multiply(MATRIX *mat1, MATRIX *mat2, MATRIX *mat3)
{
  double **scan1;
  int n1,m1;

  double **scan2;
  int n2,m2;

  double **scan3;
  int n3,m3;

  int i,j,k;

  scan1 = mat1->mat;
  n1    = mat1->n;
  m1    = mat1->m;

  scan2 = mat2->mat;
  n2    = mat2->n;
  m2    = mat2->m;

  scan3 = mat3->mat;
  n3    = mat3->n;
  m3    = mat3->m;

  if (n1 != n2 || m1 != m3 || m2 != n3)
  {
    fprintf(stderr,"matrix_multiply: matrix dimensions (%d x %d, %d x %d, %d x %d) incompatible\n",n1,m1,n2,m2,n3,m3);
    exit(1);
  }

  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < m1; j++)
    {
      double sum = 0.0;

      for (k = 0; k < m2; k++)
      {
        sum += scan2[i][k] * scan3[k][j];
      }

      scan1[i][j] = sum;
    }
  }
}

double matrix_det(MATRIX *mat1)
{
  double **scan1;
  int n1,m1;

  int i;

  double sum;

  scan1 = mat1->mat;
  n1    = mat1->n;
  m1    = mat1->m;

  if (n1 != m1)
  {
    fprintf(stderr,"matrix_det: matrix dimensions (%d x %d) incompatible\n",n1,m1);
    exit(1);
  }

  sum = 0.0;

  if (n1 == 2)
  {
    sum = scan1[0][0] * scan1[1][1] - scan1[1][0] * scan1[0][1];
  } else {
    double parity;
    MATRIX *reduced;
    double **scan;

    parity = 1.0;
    reduced = matrix_new(n1-1,m1-1);
    scan = reduced->mat;

    for (i = 0; i < n1; i++)
    {
      double sub;
      int j,k;

      for (j = 1; j < n1; j++)
      {
        for (k = 0; k < m1; k++)
        {
          if (k < i)
          {
            scan[j-1][k] = scan1[j][k];
          } else if (k > i) {
            scan[j-1][k-1] = scan1[j][k];
          }
        }
      }

      sub = matrix_det(reduced);

      sum += parity * scan1[0][i] * sub;

      parity *= -1.0;
    }

    matrix_free(reduced);
  }

  return sum;
}

void matrix_print(MATRIX *mat1)
{
  double **scan1;
  int n1,m1;

  int i,j;

  scan1 = mat1->mat;
  n1    = mat1->n;
  m1    = mat1->m;

  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < m1; j++)
    {
      printf("%9.6f ",scan1[i][j]);
    }

    printf("\n");
  }
}

void matrix_free(MATRIX *mat)
{
  double **scan;
  int n,m;

  int i;

  scan = mat->mat;
  n    = mat->n;
  m    = mat->m;

  for (i = 0; i < n; i++)
  {
    free(scan[i]);
  }

  free(scan);
  free(mat);
}

VECTOR *vector_new(int n)
{
  VECTOR *vec;

  vec = (VECTOR *)malloc(sizeof(VECTOR));

  vec->vec = (double *)malloc(n*sizeof(double));
  vec->n   = n;

  return(vec);
}

void vector_copy(VECTOR *out, VECTOR *in)
{
  int i,n;

  n = in->n;

  for (i = 0; i < n; i++)
  {
    out->vec[i] = in->vec[i];
  }
}

void vector_free(VECTOR *vec)
{
  free(vec->vec);
  free(vec);
}

ROTATE *rotate_new(int m, int n)
{
  ROTATE *rot;

  int total;

  rot = (ROTATE *)malloc(sizeof(ROTATE));

  total = (m * (2*n - m - 1)) / 2;

  rot->rot = vector_new(total);
  rot->n   = n;
  rot->m   = m;

  return(rot);
}

void rotate_copy(ROTATE *out, ROTATE *in)
{
  vector_copy(out->rot,in->rot);
  out->m = in->m;
  out->n = in->n;
}

void rotate_free(ROTATE *rot)
{
  vector_free(rot->rot);
  free(rot);
}
