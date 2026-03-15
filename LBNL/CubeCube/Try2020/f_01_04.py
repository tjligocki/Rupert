#!/usr/bin/env python3

from math import *

sign = lambda x: -1 if x < 0 else (1 if x > 0 else 0)

def f(n,t):
  if n == 0:
    return fabs(cos(t[0])*cos(t[1])*cos(t[2]))
  elif n == 1:
    return fabs(cos(t[0])*cos(t[1])*sin(t[2]))
  elif n == 2:
    return fabs(cos(t[0])*sin(t[1]))
  else:
    return fabs(sin(t[0]))

eps = 1e-10

def df(n,m,t):
  tc = t[:]

  tc[m] -= eps
  fm = f(n,tc)

  tc[m] += 2*eps
  fp = f(n,tc)

  return (fp-fm)/(2*eps)

def fa(t):
  return [f(i,t) for i in range(4)]

def ff(t):
  return max(fa(t))
