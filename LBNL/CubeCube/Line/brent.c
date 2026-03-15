#include "define.h"

/*
 *  The following is an implementation of "Brent's Method"
 *  for one-dimensional root finding. Pseudo-code for this
 *  algorithm can be found on p. 253 of "Numerical Recipes"
 *  ISBN 0-521-30811-9
 */

double brentRootFinder(double ina,
                       double inb,
                       double (*f)(double))
{
  double tol = 1.0e-12;

  /* Max allowed iterations and floating point precision */
  unsigned int MAXITER = 100;
  double       EPS     = 3.0e-15;

  unsigned int i;
  double a;
  double b;
  double c,fa,fb,fc;
  double d,e;
  double tol1,xm;
  double p,q,r,s;

  a = ina;
  b = inb;

  fa = f(a);
  fb = f(b);

  /* Init these to be safe */
  c = d = e = 0.0;

  if (fb*fa > 0)
    {
      fprintf(stderr,"Root wasn't bracketed:\n");
      fprintf(stderr,"  %13.6e -> %13.6e; ",a,fa);
      fprintf(stderr,"  %13.6e -> %13.6e\n",b,fb);
      exit(1);
    }

  fc = fb;

  for (i = 0; i < MAXITER; i++)
    {
      if (fb*fc > 0)
        {
          /* Rename a, b, c, and adjust bounding interval d */
          c = a;
          fc  = fa;
          d = b - a;
          e = d;
        }

      if (fabs(fc) < fabs(fb))
        {
          a = b;
          b = c;
          c = a;
          fa  = fb;
          fb  = fc;
          fc  = fa;
        }

      /* Convergence check */
      tol1  = 2.0 * EPS * fabs(b) + 0.5 * tol;
      xm    = 0.5 * (c - b);

      if (fabs(xm) <= tol1 || fb == 0.0)
        {
          break;
        }

      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
        {
          /* Attempt inverse quadratic interpolation */
          s = fb / fa;
          if (a == c)
            {
              p = 2.0 * xm * s;
              q = 1.0 - s;
            }
          else
            {
              q = fa / fc;
              r = fb / fc;
              p = s * (2.0 * xm * q * (q-r) - (b-a) * (r-1.0));
              q = (q-1.0) * (r-1.0) * (s-1.0);
            }

          /* Check whether in bounds */
          if (p > 0) q = -q;

          p = fabs(p);

          if (2.0 * p < fmin(3.0*xm*q-fabs(tol1*q),fabs(e*q)))
            {
              /* Accept interpolation */
              e = d;
              d = p / q;
            }
          else
            {
              /* Interpolation failed, use bisection */
              d = xm;
              e = d;
            }
        }
      else
        {
          /* Bounds decreasing too slowly, use bisection */
          d = xm;
          e = d;
        }

      /* Move last best guess to a */
      a = b;
      fa  = fb;

      /* Evaluate new trial root */
      if (fabs(d) > tol1)
        {
          b = b + d;
        }
      else
        {
          if (xm < 0) b = b - tol1;
          else        b = b + tol1;
        }

      fb = f(b);
    }

  if (i >= MAXITER)
    {
      fprintf(stderr,"brentRootFinder: exceeding maximum iterations: %d\n",MAXITER);
    }

  return b;
}
