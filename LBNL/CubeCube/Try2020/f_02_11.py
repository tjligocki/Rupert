#!/usr/bin/env python3

from math import *

sign = lambda x: -1 if x < 0 else (1 if x > 0 else 0)

def f(n,t):
  if n == 0:
    return fabs(cos(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4])) + fabs(-sin(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
  elif n == 1:
    return fabs(sin(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4])) + fabs( cos(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
  elif n == 2:
    return fabs(sin(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
  elif n == 3:
    return fabs(sin(t[2])*cos(t[3])*cos(t[4]))
  elif n == 4:
    return fabs(sin(t[3])*cos(t[4]))
  elif n == 5:
    return fabs(sin(t[4]))
  elif n == 6:
    return fabs(sin(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
  elif n == 7:
    return fabs(sin(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
  elif n == 8:
    return fabs(sin(t[7])*cos(t[8])*cos(t[9]))
  elif n == 9:
    return fabs(sin(t[8])*cos(t[9]))
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

"""
def df(n,m,t):
  if n == 0:
    if m == 0:
      return -sin(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4])*sign( cos(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4])) - cos(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9])*sign(-sin(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 1:
      return -cos(t[0])*sin(t[1])*cos(t[2])*cos(t[3])*cos(t[4])*sign(cos(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 2:
      return -cos(t[0])*cos(t[1])*sin(t[2])*cos(t[3])*cos(t[4])*sign(cos(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 3:
      return -cos(t[0])*cos(t[1])*cos(t[2])*sin(t[3])*cos(t[4])*sign(cos(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 4:
      return -cos(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*sin(t[4])*sign(cos(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 5:
      return  sin(t[0])*sin(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9])*sign(-sin(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 6:
      return  sin(t[0])*cos(t[5])*sin(t[6])*cos(t[7])*cos(t[8])*cos(t[9])*sign(-sin(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 7:
      return  sin(t[0])*cos(t[5])*cos(t[6])*sin(t[7])*cos(t[8])*cos(t[9])*sign(-sin(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 8:
      return  sin(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*sin(t[8])*cos(t[9])*sign(-sin(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    else:
      return  sin(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*sin(t[9])*sign(-sin(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
  elif n == 1:
    if m == 0:
      return  cos(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4])*sign( sin(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4])) - sin(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9])*sign( cos(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 1:
      return -sin(t[0])*sin(t[1])*cos(t[2])*cos(t[3])*cos(t[4])*sign(-sin(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 2:
      return -sin(t[0])*cos(t[1])*sin(t[2])*cos(t[3])*cos(t[4])*sign(-sin(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 3:
      return -sin(t[0])*cos(t[1])*cos(t[2])*sin(t[3])*cos(t[4])*sign(-sin(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 4:
      return -sin(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*sin(t[4])*sign(-sin(t[0])*cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 5:
      return -cos(t[0])*sin(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9])*sign( cos(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 6:
      return -cos(t[0])*cos(t[5])*sin(t[6])*cos(t[7])*cos(t[8])*cos(t[9])*sign( cos(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 7:
      return -cos(t[0])*cos(t[5])*cos(t[6])*sin(t[7])*cos(t[8])*cos(t[9])*sign( cos(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 8:
      return -cos(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*sin(t[8])*cos(t[9])*sign( cos(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    else:
      return -cos(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*sin(t[9])*sign( cos(t[0])*cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
  elif n == 2:
    if m == 0:
      return 0.0
    elif m == 1:
      return  cos(t[1])*cos(t[2])*cos(t[3])*cos(t[4])*sign(sin(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 2:
      return -sin(t[1])*sin(t[2])*cos(t[3])*cos(t[4])*sign(sin(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 3:
      return -sin(t[1])*cos(t[2])*sin(t[3])*cos(t[4])*sign(sin(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 4:
      return -sin(t[1])*cos(t[2])*cos(t[3])*sin(t[4])*sign(sin(t[1])*cos(t[2])*cos(t[3])*cos(t[4]))
    elif m == 5:
      return 0.0
    elif m == 6:
      return 0.0
    elif m == 7:
      return 0.0
    elif m == 8:
      return 0.0
    else:
      return 0.0
  elif n == 3:
    if m == 0:
      return 0.0
    elif m == 1:
      return 0.0
    elif m == 2:
      return  cos(t[2])*cos(t[3])*cos(t[4])*sign(sin(t[2])*cos(t[3])*cos(t[4]))
    elif m == 3:
      return -sin(t[2])*sin(t[3])*cos(t[4])*sign(sin(t[2])*cos(t[3])*cos(t[4]))
    elif m == 4:
      return -sin(t[2])*cos(t[3])*sin(t[4])*sign(sin(t[2])*cos(t[3])*cos(t[4]))
    elif m == 5:
      return 0.0
    elif m == 6:
      return 0.0
    elif m == 7:
      return 0.0
    elif m == 8:
      return 0.0
    else:
      return 0.0
  elif n == 4:
    if m == 0:
      return 0.0
    elif m == 1:
      return 0.0
    elif m == 2:
      return 0.0
    elif m == 3:
      return  cos(t[3])*cos(t[4])*sign(sin(t[3])*cos(t[4]))
    elif m == 4:
      return -sin(t[3])*sin(t[4])*sign(sin(t[3])*cos(t[4]))
    elif m == 5:
      return 0.0
    elif m == 6:
      return 0.0
    elif m == 7:
      return 0.0
    elif m == 8:
      return 0.0
    else:
      return 0.0
  elif n == 5:
    if m == 0:
      return 0.0
    elif m == 1:
      return 0.0
    elif m == 2:
      return 0.0
    elif m == 3:
      return 0.0
    elif m == 4:
      return cos(t[4])*sign(sin(t[4]))
    elif m == 5:
      return 0.0
    elif m == 6:
      return 0.0
    elif m == 7:
      return 0.0
    elif m == 8:
      return 0.0
    else:
      return 0.0
  elif n == 6:
    if m == 0:
      return 0.0
    elif m == 1:
      return 0.0
    elif m == 2:
      return 0.0
    elif m == 3:
      return 0.0
    elif m == 4:
      return 0.0
    elif m == 5:
      return  cos(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9])*sign(sin(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 6:
      return -sin(t[5])*sin(t[6])*cos(t[7])*cos(t[8])*cos(t[9])*sign(sin(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 7:
      return -sin(t[5])*cos(t[6])*sin(t[7])*cos(t[8])*cos(t[9])*sign(sin(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 8:
      return -sin(t[5])*cos(t[6])*cos(t[7])*sin(t[8])*cos(t[9])*sign(sin(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    else:
      return -sin(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*sin(t[9])*sign(sin(t[5])*cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
  elif n == 7:
    if m == 0:
      return 0.0
    elif m == 1:
      return 0.0
    elif m == 2:
      return 0.0
    elif m == 3:
      return 0.0
    elif m == 4:
      return 0.0
    elif m == 5:
      return 0.0
    elif m == 6:
      return  cos(t[6])*cos(t[7])*cos(t[8])*cos(t[9])*sign(sin(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 7:
      return -sin(t[6])*sin(t[7])*cos(t[8])*cos(t[9])*sign(sin(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    elif m == 8:
      return -sin(t[6])*cos(t[7])*sin(t[8])*cos(t[9])*sign(sin(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
    else:
      return -sin(t[6])*cos(t[7])*cos(t[8])*sin(t[9])*sign(sin(t[6])*cos(t[7])*cos(t[8])*cos(t[9]))
  elif n == 8:
    if m == 0:
      return 0.0
    elif m == 1:
      return 0.0
    elif m == 2:
      return 0.0
    elif m == 3:
      return 0.0
    elif m == 4:
      return 0.0
    elif m == 5:
      return 0.0
    elif m == 6:
      return 0.0
    elif m == 7:
      return  cos(t[7])*cos(t[8])*cos(t[9])*sign(sin(t[7])*cos(t[8])*cos(t[9]))
    elif m == 8:
      return -sin(t[7])*sin(t[8])*cos(t[9])*sign(sin(t[7])*cos(t[8])*cos(t[9]))
    else:
      return -sin(t[7])*cos(t[8])*sin(t[9])*sign(sin(t[7])*cos(t[8])*cos(t[9]))
  elif n == 9:
    if m == 0:
      return 0.0
    elif m == 1:
      return 0.0
    elif m == 2:
      return 0.0
    elif m == 3:
      return 0.0
    elif m == 4:
      return 0.0
    elif m == 5:
      return 0.0
    elif m == 6:
      return 0.0
    elif m == 7:
      return 0.0
    elif m == 8:
      return  cos(t[8])*cos(t[9])*sign(sin(t[8])*cos(t[9]))
    else:
      return -sin(t[8])*sin(t[9])*sign(sin(t[8])*cos(t[9]))
  else:
    if m == 0:
      return 0.0
    elif m == 1:
      return 0.0
    elif m == 2:
      return 0.0
    elif m == 3:
      return 0.0
    elif m == 4:
      return 0.0
    elif m == 5:
      return 0.0
    elif m == 6:
      return 0.0
    elif m == 7:
      return 0.0
    elif m == 8:
      return 0.0
    else:
      return cos(t[9])*sign(sin(t[9]))
"""

def fa(t):
  return [f(i,t) for i in range(11)]

def ff(t):
  return max(fa(t))
