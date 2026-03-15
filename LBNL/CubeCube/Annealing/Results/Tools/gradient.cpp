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

  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      Number grad;

      double dt = 1e-45;
      rupert.Gradient(grad,i,j,dt);
    
      printf("  %3d %3d %5.0le  ",i,j,dt);
      mpfr_out_str(stdout,10,0,grad.Get(),GMP_RNDN);
      printf("\n");
    }
  }

  exit(0);
}
