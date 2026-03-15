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

  eps = 1e-8
  if argc > 3:
    eps = float(argv[3])

  rng = np.random.default_rng()

  total = rupert.total(m,n)

  #rot = [math.pi / 4.0]
  #rot = [math.pi / 4.0, 0.000000, 4 * math.atan(math.sqrt(5.0 - 2.0*math.sqrt(6.0)))]
  #print(rot,'->',end=' ')
  #print('  ',1.0/func(rot,1,2,eps))
  #print()

  mins = []
  for i in range(1000):
    rot = math.pi * (2.0*rng.random(total) - 1.0)

    # print(rot,'->',end=' ')

    options = {}
    options['disp'] = False
    options['maxiter'] = 10000
    result = scipy.optimize.minimize(func,rot,args=(m,n,eps),method='CG',tol=1e-16,options=options)
    mins.append(result.get('fun'))

    # print(result.get('x'),'  ',1.0/result.get('fun'))

  print(1.0/min(mins))


def func(x,m,n,eps):
  # print("func:",x,m,n)

  rup = rupert.rupert(m,n,x)

  return rup.value_approx(eps)


if __name__ == "__main__":
  main(len(sys.argv),sys.argv)
