#include "rupert.h"

main(int argc, char** argv) {
  int m,n;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  if (m > n) {
    fprintf(stderr,"m (%d) must be <= n (%d)\n",m,n);
    exit(1);
  }

  Rupert rupert(m,n);

  rupert.Read();
  rupert.Compute();

  printf("%3d %3d     ",m,n);
  mpfr_out_str(stdout,10,0,rupert.GetSide().Get(),GMP_RNDN);
  printf("\n");

  exit(0);
}
