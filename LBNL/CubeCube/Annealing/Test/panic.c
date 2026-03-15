#include "define.h"
#include "global.h"

static int dopanic = 0;
static PARAMETERS *saveparams;
static CHAIN *savecurchain;
static CHAIN *savebestchain;
static int *savelength;
static REAL *saveemin;
static REAL *savee0;

void initpanic(PARAMETERS *params, CHAIN *curchain, CHAIN *bestchain,
         int *length, REAL *emin, REAL *e0)
{
  saveparams = params;
  savecurchain = curchain;
  savebestchain = bestchain;
  savelength = length;
  saveemin = emin;
  savee0 = e0;

  dopanic = 1;
}

void panic(int signal)
{
  if (dopanic == 1) {
    printchain(stdout,savecurchain,*savelength,saveparams);
    fprintf(stdout,"\n");
    fflush(stdout);

    printchain(stdout,savebestchain,*savelength,saveparams);
    fprintf(stdout,"\n");
    fflush(stdout);

    fprintf(stderr,"\n");
    fprintf(stderr,"Last chain: %34.27Le\n",*savee0);
    fprintf(stderr,"Best chain: %34.27Le (%34.27Le)\n",*saveemin,energy(savebestchain,*savelength,saveparams));
  }

  exit(1);
}

void endpanic()
{
  dopanic = 0;
}
