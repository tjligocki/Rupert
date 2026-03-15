#include <cstdio>
#include <cmath>

#include "angle.h"
#include "matrix.h"

int main(int argc, char **argv)
{
  MATRIX mat1(3,3);
  MATRIX mat2(3,3);
  MATRIX mat3(3,3);
  MATRIX mat4(3,3);
  MATRIX mat5(3,3);

  double a01,a02,a12;

  double det;

  printf("Starting to write C code to work on the Generalized Rupert Problem.\n");
  printf("\n");

  a01 = M_PI/4;
  a02 = 0.0;
  a12 = 4 * atan(sqrt(5-2*sqrt(6)));

  printf("a01: %9.6f (%9.4f)\n",a01,rad2deg(a01));
  printf("a02: %9.6f (%9.4f)\n",a02,rad2deg(a02));
  printf("a12: %9.6f (%9.4f)\n",a12,rad2deg(a12));
  printf("\n");

  mat4.identity();
  mat4.print();
  printf("\n");

  mat4.rotate(a12,1,2);
  mat4.print();
  printf("\n");

  mat4.rotate(a02,0,2);
  mat4.print();
  printf("\n");

  mat4.rotate(a01,0,1);
  mat4.print();
  printf("\n");

  return 0;
}
