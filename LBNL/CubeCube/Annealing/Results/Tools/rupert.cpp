#include "rupert.h"

Number::Number() {
  mpfr_init2(number,NUMBER_PRECISION);
}

Number::Number(const Number& in) {
  mpfr_init2(number,NUMBER_PRECISION);
  mpfr_set(number,((Number &)in).Get(),GMP_RNDN);
}

Number::~Number() {
  mpfr_clear(number);
}

const Number& Number::operator=(const Number& in) {
  if (this != &in) {
    mpfr_set(number,((Number &)in).Get(),GMP_RNDN);
  }

  return *this;
}


Rupert::Rupert(int m0, int n0) {
  m = m0;
  n = n0;

  matrix.resize(n);
  for (int i = 0; i < n; i++) {
    matrix[i].resize(n);
  }

  sum.resize(n);
}

Rupert::Rupert(const Rupert& in) {
  m = in.m;
  n = in.n;

  matrix.resize(n);
  for (int i = 0; i < n; i++) {
    matrix[i].resize(n);
    for (int j = 0; j < n; j++) {
      matrix[i][j] = in.matrix[i][j];
    }
  }

  sum.resize(n);
  for (int i = 0; i < n; i++) {
    sum[i] = in.sum[i];
  }

  maxSum = in.maxSum;
}

Rupert::~Rupert() {
}

const Rupert& Rupert::operator=(const Rupert& in) {
  if (this != &in) {
    m = in.m;
    n = in.n;

    matrix.resize(n);
    for (int i = 0; i < n; i++) {
      matrix[i].resize(n);
      for (int j = 0; j < n; j++) {
        matrix[i][j] = in.matrix[i][j];
      }
    }

    sum.resize(n);
    for (int i = 0; i < n; i++) {
      sum[i] = in.sum[i];
    }

    maxSum = in.maxSum;
  }

  return *this;
}

void Rupert::Read() {
  mpfr_set_d(maxSum.Get(),0.0,GMP_RNDN);
  for (int i = 0; i < n; i++) {
    Number absSum;

    scanf("%*s");

    for (int j = 0; j < n; j++) {
      mpfr_inp_str((matrix[i][j]).Get(),stdin,10,GMP_RNDN);
    }

    mpfr_inp_str((sum[i]).Get(),stdin,10,GMP_RNDN);

    mpfr_abs(absSum.Get(),(sum[i]).Get(),GMP_RNDN);
    if (mpfr_cmp(absSum.Get(),maxSum.Get()) > 0) {
      mpfr_set(maxSum.Get(),absSum.Get(),GMP_RNDN);
    }
  }
}

void Rupert::Write() {
  printf("\n");

  for (int i = 0; i < n; i++) {
    printf("%3d",i);

    for (int j = 0; j < n; j++) {
      printf(" ");
      mpfr_out_str(stdout,10,0,(matrix[i][j]).Get(),GMP_RNDN);
    }

    printf(" ");
    mpfr_out_str(stdout,10,0,(sum[i]).Get(),GMP_RNDN);
    printf("\n");
  }
  printf("\n");

  mpfr_out_str(stdout,10,0,maxSum.Get(),GMP_RNDN);
  printf("\n");

  printf("\n");
}

void Rupert::DotRows(Number& dot, int i1, int i2) {
  int i;
  Number mul;

  mpfr_set_d(dot.Get(),0.0,GMP_RNDN);
  for (int j = 0; j < n; j++) {
    mpfr_mul(mul.Get(),(matrix[i1][j]).Get(),(matrix[i2][j]).Get(),GMP_RNDN);
    mpfr_add(dot.Get(),dot.Get(),mul.Get(),GMP_RNDN);
  }
}

void Rupert::Compute() {
  mpfr_set_d(maxSum.Get(),0.0,GMP_RNDN);

  for (int i = 0; i < n; i++) {
    mpfr_set_d(sum[i].Get(),0.0,GMP_RNDN);

    for (int j = 0; j < m; j++) {
      Number absEntry;

      mpfr_abs(absEntry.Get(),matrix[i][j].Get(),GMP_RNDN);
      mpfr_add(sum[i].Get(),sum[i].Get(),absEntry.Get(),GMP_RNDN);
    }

    if (mpfr_cmp(sum[i].Get(),maxSum.Get()) > 0) {
      maxSum = sum[i];
    }
  }

  mpfr_ui_div(side.Get(),1,maxSum.Get(),GMP_RNDN);
}

void Rupert::OrthoNormal() {
  for (int i = 0; i < n; i++) {
    Number norm;

    DotRows(norm,i,i);
    mpfr_sqrt(norm.Get(),norm.Get(),GMP_RNDN);

    for (int j = 0; j < n; j++) {
      mpfr_div(matrix[i][j].Get(),matrix[i][j].Get(),norm.Get(),GMP_RNDN);
    }

    for (int i2 = i+1; i2 < n; i2++) {
      Number dot;

      DotRows(dot,i,i2);
      for (int j = 0; j < n; j++) {
        Number proj;

        mpfr_mul(proj.Get(),dot.Get(),matrix[i][j].Get(),GMP_RNDN);
        mpfr_sub(matrix[i2][j].Get(),matrix[i2][j].Get(),proj.Get(),GMP_RNDN);
      }
    }
  }
}

void Rupert::Gradient(Number& grad, int i1, int i2, double dt) {
  Rupert delta(m,n);
  Number angle;
  Number fplus,fminus;

  mpfr_set_d(angle.Get(),-dt,GMP_RNDN);
  Transform(delta,i1,i2,angle);
  delta.Compute();
  fminus = delta.GetSide();

  mpfr_neg(angle.Get(),angle.Get(),GMP_RNDN);
  Transform(delta,i1,i2,angle);
  delta.Compute();
  fplus = delta.GetSide();

  mpfr_sub(grad.Get(),fplus.Get(),fminus.Get(),GMP_RNDN);
  mpfr_mul_ui(angle.Get(),angle.Get(),2,GMP_RNDN);
  mpfr_div(grad.Get(),grad.Get(),angle.Get(),GMP_RNDN);
}

void Rupert::Transform(Rupert& trans, int i1, int i2, const Number& angle) {
  Number ct,st;

  mpfr_cos(ct.Get(),((Number &)angle).Get(),GMP_RNDN);
  mpfr_sin(st.Get(),((Number &)angle).Get(),GMP_RNDN);

  trans = *this;

  for (int j = 0; j < n; j++) {
    Number t1,t2;

    mpfr_mul(t1.Get(),ct.Get(),matrix[i1][j].Get(),GMP_RNDN); 
    mpfr_mul(t2.Get(),st.Get(),matrix[i2][j].Get(),GMP_RNDN); 
    mpfr_sub(trans.matrix[i1][j].Get(),t1.Get(),t2.Get(),GMP_RNDN);

    mpfr_mul(t1.Get(),st.Get(),matrix[i1][j].Get(),GMP_RNDN); 
    mpfr_mul(t2.Get(),ct.Get(),matrix[i2][j].Get(),GMP_RNDN); 
    mpfr_add(trans.matrix[i2][j].Get(),t1.Get(),t2.Get(),GMP_RNDN);
  }
}

vector<vector<Number> >& Rupert::GetMatrix() {
  return matrix;
}

vector<Number>& Rupert::GetSum() {
  return sum;
}

Number& Rupert::GetMaxSum() {
  return maxSum;
}

Number& Rupert::GetSide() {
  return side;
}
