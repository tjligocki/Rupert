#include <stdio.h>
#include <math.h>

#include "matrix.h"

int main(int argc, char **argv)
{
  MATRIX *mat1;
  MATRIX *mat2;
  MATRIX *mat3;
  MATRIX *mat4;
  MATRIX *mat5;

  double a01,a02,a12;

  double det;

  printf("Starting to write C code to work on the Generalized Rupert Problem.\n");
  printf("\n");

  mat1 = matrix_new(3,3);
  mat2 = matrix_new(3,3);
  mat3 = matrix_new(3,3);
  mat4 = matrix_new(3,3);
  mat5 = matrix_new(3,3);

  a01 = M_PI/4;
  a02 = 0.0;
  a12 = 4 * atan(sqrt(5-2*sqrt(6)));

  printf("a01: %9.6f (%9.4f)\n",a01,rad2deg(a01));
  printf("a02: %9.6f (%9.4f)\n",a02,rad2deg(a02));
  printf("a12: %9.6f (%9.4f)\n",a12,rad2deg(a12));
  printf("\n");

  matrix_identity(mat4);
  matrix_print(mat4);
  printf("\n");

  matrix_rotate(mat4,mat4,a12,1,2);
  matrix_print(mat4);
  printf("\n");

  matrix_rotate(mat4,mat4,a02,0,2);
  matrix_print(mat4);
  printf("\n");

  matrix_rotate(mat4,mat4,a01,0,1);
  matrix_print(mat4);
  printf("\n");

  matrix_free(mat1);
  matrix_free(mat2);
  matrix_free(mat3);
  matrix_free(mat4);
  matrix_free(mat5);

  return 0;
}
