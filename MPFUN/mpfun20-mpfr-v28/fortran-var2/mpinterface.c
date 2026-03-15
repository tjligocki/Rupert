/*
****************************************************************************

MPFUN2020: A thread-safe arbitrary precision package
MPFUN-MPFR version

Revision date: 9 Jan 2022

AUTHOR:
David H. Bailey
Lawrence Berkeley National Lab (retired) and University of California, Davis
Email: dhbailey@lbl.gov
 
COPYRIGHT AND DISCLAIMER:
All software in this package (c) 2022 David H. Bailey. By downloading or using this software you agree to the copyright, disclaimer and license agreement in the accompanying file DISCLAIMER.txt.

PURPOSE OF THIS ROUTINE:
This is interface between the MPFUN-MPFR Fortran modules and some of the C functions in the MPFR package. The main issue is to convert Fortran calls, which are always by address, to
MPFR calls, which in most cases include at least one variable that is call by value.

These routines are listed alphabetically.

*/

#include <stdio.h>
#include <mpfr.h>

void mpfrabs (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_abs (r1, r2, n1);
}

void mpfracos (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_acos (r1, r2, n1);
}

void mpfracosh (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_acosh (r1, r2, n1);
}

void mpfradd (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_add (r1, r2, r3, n1);
}

void mpfraddd (mpfr_t r1, mpfr_t r2, double *d, mpfr_rnd_t *irnd)
{double d1; int n1;
 d1 = *d; n1 = *irnd;
 mpfr_add_d (r1, r2, d1, n1);
}

void mpfragm (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_agm (r1, r2, r3, n1);
}

void mpfrai (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_ai (r1, r2, n1);
}

void mpfrasin (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_asin (r1, r2, n1);
}

void mpfrasinh (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_asinh (r1, r2, n1);
}

void mpfratan (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_atan (r1, r2, n1);
}

void mpfratanh (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_atanh (r1, r2, n1);
}

void mpfratan2 (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_atan2 (r1, r2, r3, n1);
}

void mpfrbesseljn (mpfr_t r1, int *iexp, mpfr_t r2, mpfr_rnd_t *irnd)
{long int n1; int n2;
 n1 = *iexp; n2 = *irnd;
 mpfr_jn (r1, n1, r2, n2);
}

void mpfrbesselyn (mpfr_t r1, int *iexp, mpfr_t r2, mpfr_rnd_t *irnd)
{long int n1; int n2;
 n1 = *iexp; n2 = *irnd;
 mpfr_yn (r1, n1, r2, n2);
}

int mpfrcmp (mpfr_t r1, mpfr_t r2)
{int i1; 
i1 = mpfr_cmp (r1, r2); return i1;
}

int mpfrcmpd (mpfr_t r1, double *d)
{int i1; double d1; d1 = *d;
i1 = mpfr_cmp_d (r1, d1); return i1;
}

void mpfrconsteuler (mpfr_t r1, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_const_euler (r1, n1);
}

void mpfrconstlog2 (mpfr_t r1, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_const_log2 (r1, n1);
}

void mpfrconstpi (mpfr_t r1, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_const_pi (r1, n1);
}

void mpfrcos (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_cos (r1, r2, n1);
}

void mpfrcosh (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_cosh (r1, r2, n1);
}

void mpfrdigamma (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_digamma (r1, r2, n1);
}

void mpfrdiv (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_div (r1, r2, r3, n1);
}

void mpfrddiv (mpfr_t r1, double *d, mpfr_t r2, mpfr_rnd_t *irnd)
{double d1; int n1;
 d1 = *d; n1 = *irnd;
 mpfr_d_div (r1, d1, r2, n1);
}

void mpfrdsub (mpfr_t r1, double *d, mpfr_t r2, mpfr_rnd_t *irnd)
{double d1; int n1;
 d1 = *d; n1 = *irnd;
 mpfr_d_sub (r1, d1, r2, n1);
}

void mpfrdivd (mpfr_t r1, mpfr_t r2, double *d, mpfr_rnd_t *irnd)
{double d1; int n1;
 d1 = *d; n1 = *irnd;
 mpfr_div_d (r1, r2, d1, n1);
}

void mpfreint (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_eint (r1, r2, n1);
}

void mpfrerf (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_erf (r1, r2, n1);
}

void mpfrerfc (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_erfc (r1, r2, n1);
}

void mpfrexp (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_exp (r1, r2, n1);
}

void mpfrfmod (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_fmod (r1, r2, r3, n1);
}

void mpfrfms (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_t r4, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_fms (r1, r2, r3, r4, n1);
}

void mpfrgamma (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_gamma (r1, r2, n1);
}

void mpfrgammainc (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_gamma_inc (r1, r2, r3, n1);
}

double mpfrgetd (mpfr_t r1, mpfr_rnd_t *irnd)
{double d1; int n1; n1 = *irnd; 
d1 = mpfr_get_d (r1, n1); return d1;
}

double mpfrgetd2exp (long *exp, mpfr_t r1, mpfr_rnd_t *irnd)
{double d1; int n1; n1 = *irnd; 
d1 = mpfr_get_d_2exp (exp, r1, n1); return d1;
}

void mpfrgetstr (mpfr_t r1, char *chr1, int *nc1, mpfr_exp_t *iexp, mpfr_rnd_t *irnd)
{int n1; int n2;
 n1 = *nc1; n2 = *irnd; 
 mpfr_get_str (chr1, iexp, 10, n1, r1, n2);
}

void mpfrhypot (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_hypot (r1, r2, r3, n1);
}

void mpfrj0 (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_j0 (r1, r2, n1);
}

void mpfrj1 (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_j1 (r1, r2, n1);
}

void mpfrli2 (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_li2 (r1, r2, n1);
}

void mpfrlngamma (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_lngamma (r1, r2, n1);
}

void mpfrlog (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_log (r1, r2, n1);
}

void mpfrmax (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_max (r1, r2, r3, n1);
}

void mpfrmin (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_min (r1, r2, r3, n1);
}

void mpfrmul (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_mul (r1, r2, r3, n1);
}

void mpfrmuld (mpfr_t r1, mpfr_t r2, double *d, mpfr_rnd_t *irnd)
{double d1; int n1;
 d1 = *d; n1 = *irnd;
 mpfr_mul_d (r1, r2, d1, n1);
}

void mpfrneg (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_neg (r1, r2, n1);
}
 
void mpfrpow (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_pow (r1, r2, r3, n1);
}

void mpfrpowsi (mpfr_t r1, mpfr_t r2, int *iexp, mpfr_rnd_t *irnd)
{long int n1; int n2;
 n1 = *iexp; n2 = *irnd;
 mpfr_pow_si (r1, r2, n1, n2);
}

void mpfrprecround (mpfr_t r1, int *prec, mpfr_rnd_t *irnd)
{int n1; int n2; n1 = *prec*64; n2 = *irnd;
n1 = mpfr_prec_round (r1, n1, n2);
}

void mpfrroot (mpfr_t r1, mpfr_t r2, int *nrt, mpfr_rnd_t *irnd)
{unsigned long int n1; int n2;
 n1 = *nrt; n2 = *irnd;
 mpfr_rootn_ui (r1, r2, n1, n2);
}

void mpfrround (mpfr_t r1, mpfr_t r2)
{
mpfr_round (r1, r2);
}

void mpfrset (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_set (r1, r2, n1);
}

void mpfrsetd (mpfr_t r1, double *d, mpfr_rnd_t *irnd)
{double d1; int n1;
 d1 = *d; n1 = *irnd;
 mpfr_set_d (r1, d1, n1);
}

void mpfrsetprec (mpfr_t r1, int *prec)
{int n1; n1 = *prec*64;
mpfr_set_prec (r1, n1);
}

void mpfrsetstr (mpfr_t r1, char *chr1, mpfr_rnd_t *irnd)
{int n1;
 n1 = *irnd;
 mpfr_set_str (r1, chr1, 10, n1);
}

int mpfrsgn (mpfr_t r1)
{int i1; 
i1 = mpfr_sgn (r1); return i1;
}

void mpfrsin (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_sin (r1, r2, n1);
}

void mpfrsincos (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_sin_cos (r1, r2, r3, n1);
}

void mpfrsinh (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_sinh (r1, r2, n1);
}

void mpfrsinhcosh (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_sinh_cosh (r1, r2, r3, n1);
}

void mpfrsqrt (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_sqrt (r1, r2, n1);
}

void mpfrsub (mpfr_t r1, mpfr_t r2, mpfr_t r3, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_sub (r1, r2, r3, n1);
}

void mpfrsubd (mpfr_t r1, mpfr_t r2, double *d, mpfr_rnd_t *irnd)
{double d1; int n1;
 d1 = *d; n1 = *irnd;
 mpfr_sub_d (r1, r2, d1, n1);
}

void mpfrtan (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_tan (r1, r2, n1);
}

void mpfrtanh (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_tanh (r1, r2, n1);
}

void mpfrtrunc (mpfr_t r1, mpfr_t r2)
{
mpfr_trunc (r1, r2);
}

void mpfry0 (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_y0 (r1, r2, n1);
}

void mpfry1 (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_y1 (r1, r2, n1);
}

void mpfrzeta (mpfr_t r1, mpfr_t r2, mpfr_rnd_t *irnd)
{int n1; n1 = *irnd;
mpfr_zeta (r1, r2, n1);
}

