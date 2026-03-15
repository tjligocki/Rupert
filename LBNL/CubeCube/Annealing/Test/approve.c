#include "define.h"
#include "global.h"

int approve(REAL e0, REAL e1, REAL T, PARAMETERS *params)
{
  REAL delta,power;
  REAL prob,pick;
  int result;

  delta = e1-e0;

  result = 0;

  if (params->acceptance == BIASED) {
    if (delta >= 0) {
      power = delta/T;

      if (power < 700) {
        prob = 1/exp(power);
      } else {
        prob = 0;
      }
      pick = drand48();

      if (pick <= prob) {
        result = 1;
      }
    } else {
      result = 1;
    }
  } else {
    power = delta/T;

    if (power < 700) {
      if (power > -700) {
        prob = 1/(1+exp(power));
      } else {
        prob = 1;
      }
    } else {
      prob = 0;
    }
    pick = drand48();

    if (pick <= prob) {
      result = 1;
    }
  }

  return(result);
}
