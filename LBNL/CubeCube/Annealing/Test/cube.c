#include "define.h"
#include "global.h"

extern char *optarg;
extern int optind,opterr;

main(int argc, char **argv)
{
  int option;
  char paramname[8192];
  PARAMETERS params;
  CHAIN nextchain;
  CHAIN curchain;
  CHAIN bestchain;
  int length;
  REAL emin;
  REAL e0,e1;
  double fail;
  double try;
  double curtry;
  double tick,tock;
  int all;
  int n;
  REAL T;
  REAL sd;
  int count;
  REAL etotal;

  signal(SIGHUP,panic);
  signal(SIGINT,panic);
  signal(SIGQUIT,panic);
  signal(SIGPIPE,panic);
  signal(SIGTERM,panic);
  signal(SIGXCPU,panic);

  strcpy(paramname,"\0");

  initparams(&params);

  while ((option = getopt(argc,argv,"hp:")) != -1) {
    switch (option) {
    case 'h':
      helpparams();
      exit(1);
      break;

    case 'p':
      strcpy(paramname,optarg);
      break;

    case '?':
      fprintf(stderr,"Usage:  %s [-h] [-p paramfile]\n",argv[0]);
      helpparams();
      exit(1);
      break;
    }
  }

  if ((int)strlen(paramname) > 0) {
    readparams(paramname,&params);
  }

  if (params.outdim < params.indim) {
      fprintf(stderr,"Outer cube dimension, %d, must be larger than inner cube dimension, %d\n",params.outdim,params.indim);
  }

  initchain(&curchain,&length,&params);

  newchain(&nextchain,length,params.outdim);
  newchain(&bestchain,length,params.outdim);

  n = 0;

  e0 = energy(&curchain,length,&params);
  emin = e0;
  copychain(&curchain,&bestchain,length,params.outdim);

  if (params.T0 < 0) {
    params.T0 = emin;
  }

  T = params.T0;

  if (params.seed == 0) {
    params.seed = (int)time((time_t *)NULL);
  }

  srand48(params.seed);
  
  etotal = e0;
  count = 1;

  if (params.maxfail == -1) {
    params.maxfail = length / 5.0 * NUMFAIL;
  }

  if (params.maxtick == -1) {
    params.maxtick = length * NUMTRY;
  }

  if (params.maxtock == -1) {
    params.maxtock = params.maxtick * NUMSAVE;
  }

  initpanic(&params,&curchain,&bestchain,&length,&emin,&e0);

  fprintf(stdout,"Code info: none\n\n");
  printparams(stdout,&params);

  printchain(stdout,&curchain,length,&params);
  fprintf(stdout,"\n");
  fflush(stdout);

  all = 0;

  try = params.iteration;

  curtry = params.iteration;
  while (curtry >= params.retry) {
    switch (params.cooling) {
      case EXPONENTIAL:
        T = params.T0/exp(params.rate*params.retry);
        break;

      case LINEAR:
        T = params.T0/(1+params.rate*params.retry);
        break;

      case LOGARITHMIC:
        T = params.T0/log(M_E+params.rate*params.retry);
        break;
    }
    params.T0 = T;
    curtry -= params.retry;
  }

  tick = 0;
  tock = 0;

  fail = 0;

  switch (params.cooling) {
    case EXPONENTIAL:
      T = params.T0/exp(params.rate*curtry);
      break;

    case LINEAR:
      T = params.T0/(1+params.rate*curtry);
      break;

    case LOGARITHMIC:
      T = params.T0/log(M_E+params.rate*curtry);
      break;
  }

  fprintf(stderr,"emin+1: %13.6Le, e+1: %13.6Le, T: %10.4Le\n",emin+1,e0+1,T);

  while (fail < params.maxfail && T >= params.Tmin) {
    CHAIN cur;

    sd = params.sdrate * sqrt(T) + params.sdminimum;

    all += next(&curchain,&nextchain,length,sd,&params,&e1);

    if (approve(e0,e1,T,&params) == 1) {
      fail = 0;

      etotal += e1;
      count++;

      cur = curchain;
      curchain = nextchain;
      nextchain = cur;

      e0 = e1;

      if (e0 < emin) {
        emin = e0;
        copychain(&curchain,&bestchain,length,params.outdim);
      }

      n++;
    } else {
      fail++;
    }

    try++;
    curtry++;
    tick++;
    tock++;

    switch (params.cooling) {
      case EXPONENTIAL:
        T = params.T0/exp(params.rate*curtry);
        break;

      case LINEAR:
        T = params.T0/(1+params.rate*curtry);
        break;

      case LOGARITHMIC:
        T = params.T0/log(M_E+params.rate*curtry);
        break;
    }

    if (curtry == params.retry) {
      params.T0 = T;
      curtry = 0;
    }

    if (tick == params.maxtick || fail == params.maxfail) {
      if (params.movement == LOCAL) {
        fprintf(stderr,"emin+1: %13.6Le, e+1: %13.6Le, T: %10.4Le, no: %7.0lf, try: %12.0lf\n",emin+1,etotal/count + 1,T,fail,try);
      } else {
        fprintf(stderr,"emin+1: %13.6Le, e+1: %13.6Le, T: %10.4Le, no: %7.0lf, try: %12.0lf\n",emin+1,etotal/count + 1,T,fail,try);
      }

      tick = 0;
      all = 0;

      etotal = e0;
      count = 1;
    }

    if (tock == params.maxtock || fail == params.maxfail) {
      printchain(stdout,&curchain,length,&params);
      fprintf(stdout,"\n");
      fflush(stdout);

      tock = 0;
    }

    if ((try/10000) == floor(try/10000)) {
      orthonormalize(curchain.trans,params.outdim);
      energy(&curchain,length,&params);
    }
  }

  printchain(stdout,&bestchain,length,&params);
  fprintf(stdout,"\n");
  fflush(stdout);

  fprintf(stderr,"\n");
  fprintf(stderr,"Last chain: %34.27Le\n",e0);
  fprintf(stderr,"Best chain: %34.27Le (%34.27Le)\n",emin,energy(&bestchain,length,&params));

  fprintf(stderr,"\n");
  check(stderr,&bestchain,length,&params);
  fprintf(stderr,"\n");

  orthonormalize(curchain.trans,params.outdim);
  energy(&bestchain,length,&params);

  printchain(stdout,&bestchain,length,&params);
  fprintf(stdout,"\n");
  fflush(stdout);

  fprintf(stderr,"\n");
  fprintf(stderr,"Last chain: %34.27Le\n",e0);
  fprintf(stderr,"Best chain: %34.27Le (%34.27Le)\n",emin,energy(&bestchain,length,&params));

  fprintf(stderr,"\n");
  check(stderr,&bestchain,length,&params);
  fprintf(stderr,"\n");

  endpanic();

  return(0);
}
