typedef struct
{
  double **mat;
  int n,m;
} MATRIX;

typedef struct
{
  double *vec;
  int n;
} VECTOR;

typedef struct
{
  VECTOR *rot;
  int n,m;
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

extern VECTOR *vector_new(int n);
extern void vector_copy(VECTOR *out, VECTOR *in);
extern void vector_free(VECTOR *vec);

extern ROTATE *rotate_new(int n, int m);
extern void rotate_copy(ROTATE *out, ROTATE *in);
extern void rotate_free(ROTATE *rot);
