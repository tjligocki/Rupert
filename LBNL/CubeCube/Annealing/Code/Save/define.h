#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/time.h>

#define REAL  long double

#define ABS(x)		(((x) > 0) ? (x) : -(x))

#define DEBUG		0

#define BIASED		0
#define SYMMETRIC	1

#define UNIFORM		0
#define NORMAL		1

#define PHYSICAL	0
#define FREEDMAN	1

#define LOCAL		0
#define MULTI		1
#define GLOBAL		2

#define NUMFAIL		 1000
#define NUMTRY		 1500
#define NUMSAVE		    1

#define DECAY		0.998
#define THRES		0.80

#define GEOMETRIC	0
#define EXPONENTIAL	0
#define LINEAR		1
#define LOGARITHMIC	2

typedef struct chain CHAIN;
typedef struct parameters PARAMETERS;

struct chain {
	REAL **trans;
  REAL *insum;
  REAL maxsum;
};

struct parameters {
	int seed;
	int cooling;
	int indim,outdim;
	REAL T0;
	REAL Tmin;
	REAL rate;
	double iteration;
	double retry;
	double maxfail;
	double maxtick;
	double maxtock;
	int acceptance;
	int distribution;
	int movement;
	REAL sdminimum;
	REAL sdrate;
};
