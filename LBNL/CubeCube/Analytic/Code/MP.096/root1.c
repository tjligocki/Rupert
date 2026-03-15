#include <stdlib.h>
#include <stdio.h>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_multiroots.h"

typedef struct params
{
  size_t n;
  size_t nn;
} PARAMS;

int cubecube_f(const gsl_vector* x, void* params, gsl_vector* f)
{
  size_t n  = ((PARAMS *)params)->n;
  size_t nn = ((PARAMS *)params)->nn;

  size_t i,ii,iii;
  size_t j,jj,jjj;
  size_t k;
  size_t r;

  double* xx;
  double* ff;

  xx = (double *)malloc(nn*sizeof(*xx));
  ff = (double *)malloc(nn*sizeof(*ff));

  for (i = 0; i < nn; i++) {
    xx[i] = gsl_vector_get(x,i);
  }

#if 0
fprintf(stderr,"n = %d, nn = %d\n",n,nn);
#endif

  r = 0;
  ii = 0;

  for (i = 0; i < n-1; i++) {
#if 0
fprintf(stderr,"r = %d, ii = %d\n",r,ii);
#endif
    ff[r] = 1.0;
    for (j = i; j < n; j++) {
      ff[r] = ff[r] - xx[ii]*xx[ii];
#if 0
fprintf(stderr,"  - x[%d]*x[%d]\n",ii,ii);
#endif
      ii = ii+1;
    }
    r = r+1;
  }

#if 0
fprintf(stderr,"r = %d\n",r,ii);
#endif
#if 0
fprintf(stderr,"  + x[%d]*x[%d]\n",0,0);
#endif
  ff[r] = xx[0]*xx[0];

  for (i = 1; i < n; i++) {
    ff[r] = ff[r] - xx[i]*xx[i];
#if 0
fprintf(stderr,"  - x[%d]*x[%d]\n",i,i);
#endif
  }

  r = r+1;

  ii = 0;
  iii = n;
  for (i = 0; i < n-2; i++) {
    jj = ii + iii;
    jjj = iii - 1;
    for (j = i+1; j < n-1; j++) {
#if 0
fprintf(stderr,"r = %d\n",r,ii);
#endif
      ff[r] = 0;
      for (k = j; k < n; k++) {
        ff[r] = ff[r] + xx[ii+k-i]*xx[jj+k-j];
#if 0
fprintf(stderr,"  + x[%d]*x[%d]\n",ii+k-i,jj+k-j);
#endif
      }
      
      jj = jj + jjj;
      jjj = jjj - 1;
      
      r = r+1;
    }
    ii = ii + iii;
    iii = iii - 1;
  }

  ii = 0;
  iii = n;
  for (i = 0; i < n-1; i++) {
    ff[r] = 0;
#if 0
fprintf(stderr,"r = %d\n",r,ii);
#endif
    for (j = i; j < n-1; j++) {
      ff[r] = ff[r] + fabs(xx[ii+j-i]);
#if 0
fprintf(stderr,"  + fabs(x[%d])\n",ii+j-i);
#endif
    }
    ff[r] = ff[r] - fabs(xx[nn-1]);
#if 0
fprintf(stderr,"  - fabs(x[%d])\n",nn-1);
#endif
    
    ii = ii + iii;
    iii = iii - 1;
    
    r = r+1;
  }
#if 0
fprintf(stderr,"\n");
#endif

  for (i = 0; i < nn; i++) {
    gsl_vector_set(f,i,ff[i]);
  }

  free(xx);
  free(ff);

  return GSL_SUCCESS;
}

int print_state(size_t iter, gsl_multiroot_fsolver* s, const size_t n)
{
    size_t i;
    double cur,max;

    fprintf(stderr,"iter = %3u    ",iter);

    fprintf(stderr,"side = %23.20f    ",1.0/gsl_vector_get(s->x,n-1));

    max = 0;
    for (i = 0; i < n; i++) {
      cur = fabs(gsl_vector_get(s->f,i));

      if (cur > max) {
        max = cur;
      }
    }

    fprintf(stderr,"max(f) = %10.3e\n",max);
}

main(int argc, char** argv)
{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i,count;
  size_t iter;

  size_t n;
  size_t nn;

  PARAMS p;

  double result;

  gsl_multiroot_function f;

  gsl_vector* x;

  n = 5;

  if (argc > 1) {
    n = atoi(argv[1]);
  }

  nn = n*(n+1)/2;

  p.n  = n;
  p.nn = nn;

  f.f      = &cubecube_f;
  f.n      = nn;
  f.params = &p;

  srand48(time(0));

  x = gsl_vector_alloc(nn);
  T = gsl_multiroot_fsolver_hybrid;

  count = 0;

  do {
    s = gsl_multiroot_fsolver_alloc(T,nn);

    for (i = 0; i < nn; i++) {
      gsl_vector_set(x,i,2.0*drand48()-1.0);
    }

    gsl_multiroot_fsolver_set(s,&f,x);

    iter = 0;

#if 0
    print_state(iter,s,nn);
#endif

    do {
      iter++;

      status = gsl_multiroot_fsolver_iterate(s);

#if 0
      print_state(iter,s,nn);
#endif

      if (status) {
        break;
      }

      status = gsl_multiroot_test_residual(s->f,1e-12);
    } while (status == GSL_CONTINUE && iter < 100);

#if 0
    fprintf(stderr,"%s\n\n",gsl_strerror(status));
#endif

    if (count % 100 == 0) {
      fprintf(stderr,".");
    }

    count++;

    result = gsl_vector_get(s->x,nn-1);

    gsl_multiroot_fsolver_free(s);
  } while (status != GSL_SUCCESS);

  fprintf(stderr,"\nside (%d %d) = %23.20f\n",n-1,n,fabs(1.0/result));

  gsl_vector_free(x);

  exit(0);
}
