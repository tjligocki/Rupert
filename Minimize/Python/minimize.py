#!/usr/bin/env python3

import sys,math
import numpy as np
import scipy
import rupert

def main(argc,argv):
  m = 2
  if argc > 1:
    m = int(argv[1])

  n = 3
  if argc > 2:
    n = int(argv[2])

  rng = np.random.default_rng()

  total = rupert.total(m,n)

  mins = []
  for i in range(1000):
    rot = math.pi * (2.0*rng.random(total) - 1.0)

    #print(rot,'->',end=' ')

    result = scipy.optimize.minimize(func,rot,args=(m,n),method='Nelder-Mead')
    mins.append(result.get('fun'))

    #print(1.0/mins[-1])

  print(1.0/min(mins))
  print()


def func(x,m,n):
  rup = rupert.rupert(m,n,x)

  #print("func:",x,m,n,end=' ')
  #print("%f" % rup.value())
  #print()

  return rup.value()


if __name__ == "__main__":
  main(len(sys.argv),sys.argv)
