#include "rupert.h"

main(int argc, char** argv) {
  int m,n;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  if (m > n) {
    fprintf(stderr,"m (%d) must be <= n (%d)\n",m,n);
    exit(1);
  }

  Rupert rupertIn(m,n);

  rupertIn.Read();

  Number maxLengthError;
  Number maxOrthoError;

  mpfr_set_d(maxLengthError.Get(),0.0,GMP_RNDN);
  mpfr_set_d(maxOrthoError.Get(),0.0,GMP_RNDN);

  for (int i = 0; i < n; i++) {
    Number norm;
    Number aCur;

    int i2;

    rupertIn.DotRows(norm,i,i);
    mpfr_sqrt(norm.Get(),norm.Get(),GMP_RNDN);
    mpfr_ui_sub(norm.Get(),1,norm.Get(),GMP_RNDN);

    if (mpfr_cmp(norm.Get(),maxLengthError.Get()) > 0) {
      mpfr_set(maxLengthError.Get(),norm.Get(),GMP_RNDN);
    }

    for (int i2 = i+1; i2 < n; i2++) {
      rupertIn.DotRows(norm,i,i2);
      mpfr_abs(norm.Get(),norm.Get(),GMP_RNDN);

      if (mpfr_cmp(norm.Get(),maxOrthoError.Get()) > 0) {
        mpfr_set(maxOrthoError.Get(),norm.Get(),GMP_RNDN);
      }
    }
  }

  mpfr_out_str(stdout,10,0,maxLengthError.Get(),GMP_RNDN);
  printf(" ");

  mpfr_out_str(stdout,10,0,maxOrthoError.Get(),GMP_RNDN);
  printf(" ");

  Rupert rupertOut = rupertIn;
  rupertOut.Compute();

  Number maxSumError;
  mpfr_set_d(maxSumError.Get(),0.0,GMP_RNDN);

  for (int i = 0; i < n; i++) {
    Number diff;

    mpfr_sub(diff.Get(),rupertOut.GetSum()[i].Get(),rupertIn.GetSum()[i].Get(),GMP_RNDN);
    mpfr_abs(diff.Get(),diff.Get(),GMP_RNDN);

    if (mpfr_cmp(diff.Get(),maxSumError.Get()) > 0) {
      mpfr_set(maxSumError.Get(),diff.Get(),GMP_RNDN);
    }
  }

  mpfr_out_str(stdout,10,0,maxSumError.Get(),GMP_RNDN);
  printf("\n");

  printf("\n");

  exit(0);
}
