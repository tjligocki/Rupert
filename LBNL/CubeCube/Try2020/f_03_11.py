#!/usr/bin/env python3

from math import *

sign = lambda x: -1 if x < 0 else (1 if x > 0 else 0)

####
####

def f(n,t):
  if n == 0:
    return fabs(cos(t[0])*cos(t[1])*cos(t[2])) + fabs(cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9])*sin(t[0]))
  elif n == 1:
    return fabs(cos(t[1])*cos(t[2])*sin(t[0])) + fabs(-cos(t[0])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
  elif n == 2:
    return fabs(cos(t[2])*sin(t[1]))
  elif n == 3:
    return fabs(sin(t[2]))
  elif n == 4:
    return fabs(cos(t[3])*cos(t[4])*cos(t[5])) + fabs(-cos(t[7])*cos(t[8])*cos(t[9])*sin(t[3])*sin(t[6]))
  elif n == 5:
    return fabs(cos(t[4])*cos(t[5])*sin(t[3])) + fabs(cos(t[3])*cos(t[7])*cos(t[8])*cos(t[9])*sin(t[6]))
  elif n == 6:
    return fabs(cos(t[5])*sin(t[4]))
  elif n == 7:
    return fabs(sin(t[5]))
  elif n == 8:
    return fabs(cos(t[8])*cos(t[9])*sin(t[7]))
  elif n == 9:
    return fabs(cos(t[9])*sin(t[8]))
  else:
    return fabs(sin(t[9]))

eps = 1e-10

def df(n,m,t):
  tc = t[:]

  tc[m] -= eps
  fm = f(n,tc)

  tc[m] += 2*eps
  fp = f(n,tc)

  return (fp-fm)/(2*eps)

def fa(t):
  return [f(i,t) for i in range(11)]

def ff(t):
  return max(fa(t))
