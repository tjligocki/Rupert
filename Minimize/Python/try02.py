#!/usr/bin/env python3

import sys,math
import numpy as np
import scipy
import rupert
import timeit

def func_exact(x,m,n):
  rup = rupert.rupert(m,n,x)

  return rup.value()


def func_approx(x,m,n,eps):
  rup = rupert.rupert(m,n,x)

  return rup.value_approx(eps)


if __name__ == "__main__":
  argc = len(sys.argv)
  argv = sys.argv

  m = 2
  if argc > 1:
    m = int(argv[1])

  n = 3
  if argc > 2:
    n = int(argv[2])

  r = 100000
  if argc > 3:
    r = int(argv[3])

  eps = 1e-5

  rng = np.random.default_rng()

  total = rupert.total(m,n)

  rot = math.pi * (2.0*rng.random(total) - 1.0)

  t = timeit.timeit(
      'func_exact(rot,m,n)',
      'from __main__ import func_exact,rot,m,n',
      number=r)
  print('Time (%d,%d) (%d) exact  %d times: %0.3f' % (m,n,total,r,t))

  t = timeit.timeit(
      'func_approx(rot,m,n,eps)',
      'from __main__ import func_approx,rot,m,n,eps',
      number=r)
  print('Time (%d,%d) (%d) approx %d times: %0.3f' % (m,n,total,r,t))
