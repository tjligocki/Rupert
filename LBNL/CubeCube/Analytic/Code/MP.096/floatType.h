#if 0

#if 0

#define REAL   double
#define PRINTF "lf"
#define PRINTE "le"
#define TOLER  1e-15
#define SQRT   sqrt
#define FABS   fabs

#else

#define REAL   long double
#define PRINTF "Lf"
#define PRINTE "Le"
#define TOLER  1e-17
#define SQRT   sqrtl
#define FABS   fabsl

#endif

#else

#include <gmp.h>

#define REAL   mpf_t
#define PRINTF "lf"
#define PRINTE "le"
#define TOLER  1e-25
#define SQRT   sqrtl
#define FABS   fabsl
#define PRES   96

#endif
