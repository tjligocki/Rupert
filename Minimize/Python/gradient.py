#!/usr/bin/env python3

import sys,math
import numpy as np
import scipy
import rupert
import timeit
import signal


def main(argc,argv):
  signal.signal(signal.SIGINT, signal.default_int_handler)

  m = 2
  if argc > 1:
    m = int(argv[1])

  n = 3
  if argc > 2:
    n = int(argv[2])

  eps = 0.001
  if argc > 3:
    eps = float(argv[3])

  iter = 10000
  if argc > 4:
    iter = int(argv[4])

  maxiter = 100000
  if argc > 5:
    maxiter = int(argv[5])

  #rng = np.random.default_rng(4393483)
  rng = np.random.default_rng()

  total = rupert.total(m,n)

  rot = math.pi * (2.0*rng.random(total) - 1.0)

  i = 0
  #tprev = -1
  #t = 0

  v = 0
  grad = np.ones(n)
  t = 1.0

  try:
    while True:
      #if (i-1) % iter == 0:
      #  print('top:',rot,i-1)

      vprev = v

      rup = rupert.rupert(m,n,rot,eps)

      v     = rup.value
      grad  = rup.grad

      ''' 
      if i == 0:
        t = 1.0
      else:
        print(rot-rotprev)
        print(grad-gradprev)
        print(np.linalg.norm(grad-gradprev))

        t = abs(np.matmul((rot-rotprev),(grad-gradprev)))/np.linalg.norm(grad-gradprev)
      '''

      t = eps

      #tprev = t
      #t = scipy.optimize.minimize_scalar(func,args=(m,n,rot,grad,eps)).get('x')

      if i % iter == 0 or t*np.linalg.norm(grad) <= 1e-14 or abs(v-vprev) < 1e-15:
        print('[',end=' ')
        for r in rot:
          print('%13.10f ' % r,end='')
        print(' ]',end=' ')

        print('%7d: %13.10f (%13.6e) ' % (i,v,abs(v-vprev)),end=' ')

        print('[',end=' ')
        for g in grad:
          print('%8.5f ' % g,end='')
        print(' ]',end=' ')

        print('%13.10e ' % (t*np.linalg.norm(grad),),end='')

        print('%13.10e' % t)

      if t*np.linalg.norm(grad) > 1e-14 and i < maxiter and abs(v-vprev) >= 1e-15:
        rotprev = rot
        gradprev = grad

        #if i % iter == 0:
        #  print('up: ',rot,t,grad,(t*grad),'->',end='')

        rot = rot - t*grad

        #if i % iter == 0:
        #  print(rot,i)

        i += 1
      else:
        break

      #if (i-1) % iter == 0:
      #  print('end:',rot,i-1)
  except KeyboardInterrupt:
    pass

  print(1/rup.value,', ',sep='',end=' ')

  print('%13.10e (%13.10e %13.10e) %d' % (t*np.linalg.norm(grad),np.linalg.norm(grad),t,i),end=' ')

  print('{',end='')
  for i in range(total):
    print('%.10f' % (rot[i],),end='')

    if i < total-1:
      print(',',end='')
  print('} ',end='')

  print('{',end='')
  for i in range(total):
    print('%.10e' % (grad[i],),end='')

    if i < total-1:
      print(',',end='')
  print('} ',end='')

  print()


def func(t,m,n,rot,grad,eps):
  rup = rupert.rupert(m,n,rot-t*grad,eps)

  return rup.value


if __name__ == '__main__':
  main(len(sys.argv),sys.argv)

  sys.exit(0)
