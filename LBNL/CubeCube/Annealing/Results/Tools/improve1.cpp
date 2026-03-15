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
  rupert.OrthoNormal();
  rupert.Compute();
  rupert.Write();

  exit(0);
}
