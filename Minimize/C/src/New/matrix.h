typedef struct
{
  double **mat;
  int n,m;
} MATRIX;

typedef struct
{
  double *rot;
  int m,n;
  int total;
} ROTATE;

extern double rad2deg(double deg);
extern double deg2rad(double rad);

extern MATRIX *matrix_new(int n, int m);

extern void matrix_zero(MATRIX *mat);
extern void matrix_identity(MATRIX *mat);
extern void matrix_add(MATRIX *mat1, MATRIX *mat2, MATRIX *mat3);
extern void matrix_scalar(MATRIX *mat1, MATRIX *mat2, double a);
extern void matrix_multiply(MATRIX *mat1, MATRIX *mat2, MATRIX *mat3);
extern void matrix_rotate(MATRIX *mat1, MATRIX *mat2, double a, int i1, int i2);

extern double matrix_det(MATRIX *mat1);

extern void matrix_print(MATRIX *mat);

extern void matrix_free(MATRIX *mat);

extern ROTATE *rotate_new(int m, int n);

extern void rotate_free(ROTATE *rot);
