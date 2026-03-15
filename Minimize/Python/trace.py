#!/usr/bin/env python3

import sys,math
import numpy as np
import scipy
import rupert

def main(argc,argv):
  n = 2
  if argc > 1:
    n = int(argv[1])

  if False:
    rng = np.random.default_rng()

    tryit = {}
    tryit['x'] = [math.pi/4]
    tryit['fun'] = math.sqrt(2)/2
    traceit(tryit)
    
    for i in range(1):
      x = (2.0*rng.random(n) - 1.0)

      options = {}
      options['maxiter'] = 1000

      result = scipy.optimize.minimize(func,x,method='Nelder-Mead',callback=traceit,tol=1e-14,options=options)
      print(result)

  dx = 2/(n-1)
  eps = 1e-1

  for h in [1e-5,9e-6,8e-6]:
    for i in range(n):
      x = i*dx - 1.0

      print(x,func(x,eps),dfunc([x],eps,h)[0],dabs_approx(x,eps),abs(x))

    print()
    print()


def abs_approx(x,eps):
  return math.sqrt(x**2 + eps**2)


def dabs_approx(x,eps):
  return x/math.sqrt(x**2 + eps**2)


def func(x,eps):
  return abs_approx(x,eps)
  #return max(abs_approx(math.sin(x[0])),abs_approx(math.cos(x[0])))


def dfunc(x,eps,h):
  d = len(x)

  df = d*[0.0]

  #print("dfunc:","%20.16f" % x[0])

  for i in range(d):
    xcur = x[:]
    xcur[i] -= h 
    vlo = func(xcur[i],eps)

    #print("  %20.16f" % xcur[0])

    xcur[i] += 2*h
    vhi = func(xcur[i],eps)

    #print("  %20.16f" % xcur[0])

    #print("  %20.16f" % vlo,"%20.16f" % vhi,"%.16e" % (vhi - vlo))

    df[i] = (vhi - vlo) / (2*h)

  return df


Ncallback = 1
def traceit(intermediate_result):
  global Ncallback

  x = intermediate_result.get('x')
  v = intermediate_result.get('fun')

  df = dfunc(x)

  print("%6d" % Ncallback,x,v,df)

  Ncallback += 1


if __name__ == "__main__":
  main(len(sys.argv),sys.argv)
