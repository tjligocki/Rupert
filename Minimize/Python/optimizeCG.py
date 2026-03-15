#!/usr/bin/env python3

import sys,math,time,timeit,signal
import numpy as np
import scipy
import rupert


def main(argc,argv):
  signal.signal(signal.SIGINT, signal.default_int_handler)

  m = 2
  if argc > 1:
    m = int(argv[1])

  n = 3
  if argc > 2:
    n = int(argv[2])

  eps_start = 0.001
  if argc > 3:
    eps_start = float(argv[3])

  eps_end = 0.001
  if argc > 4:
    eps_end = float(argv[4])

  multi = 100
  if argc > 5:
    multi = int(argv[5])

  rng = np.random.default_rng()

  total = rupert.total(m,n)
  print('%3d %3d  %5d     ' % (m,n,total),end='')

  time_start = time.time()

  try:
    values = []
    for i in range(total*100):
      rot = math.pi * (2.0*rng.random(total) - 1.0)

      value1_prev = 0.0
      #value2_prev = 0.0

      eps = eps_start

      while eps >= eps_end:
        #print('eps: %.2e' % eps,end=' ')

        rv = scipy.optimize.minimize(func,rot,jac=jac,method='CG',args=(m,n,eps),options={'maxiter':100000})

        value1 = 1.0/rv.get('fun')
        grad1  = rv.get('jac')

        rot = rv.get('x')
        #rup = rupert.rupert(m,n,rot,eps)

        #value2 = 1.0/rup.value
        #grad2  = rup.grad

        #print('  ',rv)
        #print()

        #print('  %.16f ' % value1,end='')
        #if value1_prev > 0:
        #  print(' (%.10e)' % (value1-value1_prev),end='')
        #else:
        #  print('                   ',end='')

        #print(' %17.10e' % np.linalg.norm(grad1),end='')
        #print()

        #print('  %.16f ' % value2,end='')
        #if value2_prev > 0:
        #  print(' (%.10e)' % (value2-value2_prev),end='')
        #else:
        #  print('                   ',end='')

        #print(' %17.10e' % np.linalg.norm(grad2),end='')
        #print()

        value1_prev = value1
        #value2_prev = value2

        eps *= 0.1

      values.append(value1)

  except KeyboardInterrupt:
    pass

  time_end = time.time()

  print('%21.16f %21.16f' % (min(values),max(values)),end=' ')
  # print()

  print('%12.6f' % (time_end - time_start))


def func(rot,m,n,eps):
  rup = rupert.rupert(m,n,rot,eps)

  return rup.value


def jac(rot,m,n,eps):
  rup = rupert.rupert(m,n,rot,eps)

  return rup.grad


if __name__ == '__main__':
  main(len(sys.argv),sys.argv)

  sys.exit(0)
