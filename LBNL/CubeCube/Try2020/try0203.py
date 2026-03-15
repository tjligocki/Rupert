#!/usr/bin/env python3

import sys,random,time
from math import *

sign = lambda x: -1 if x < 0 else (1 if x > 0 else 0)

def value(t1,t2):
  return max(fabs(cos(t1))+fabs(-sin(t1)*cos(t2)),
             fabs(sin(t1))+fabs( cos(t1)*cos(t2)),
             fabs(                       sin(t2)))

def func(t):
  global p1,p2,v1,v2
  return value(p1+t*v1,p2+t*v2)

def value3(t1,t2):
  return (fabs(cos(t1))+fabs(-sin(t1)*cos(t2)),
          fabs(sin(t1))+fabs( cos(t1)*cos(t2)),
          fabs(                       sin(t2)))

def func3(t):
  global p1,p2,v1,v2
  return value3(p1+t*v1,p2+t*v2)

def deri3(t1,t2):
  return ((-sin(t1)*sign(cos(t1))-cos(t1)*cos(t2)*sign(-sin(t1)*cos(t2)),
                                  sin(t1)*sin(t2)*sign(-sin(t1)*cos(t2))),
          ( cos(t1)*sign(sin(t1))-sin(t1)*cos(t2)*sign( cos(t1)*cos(t2)),
                                 -cos(t1)*sin(t2)*sign( cos(t1)*cos(t2))),
          (0                                                            ,
                                          cos(t2)*sign(         sin(t2))))
   

def grad3(t):
  global p1,p2,v1,v2
  return deri3(p1+t*v1,p2+t*v2)

invphi = (sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - sqrt(5)) / 2  # 1 / phi^2

def gss(f, a, b, tol=1e-5):
  """Golden-section search.

  Given a function f with a single local minimum in
  the interval [a,b], gss returns a subset interval
  [c,d] that contains the minimum with d-c <= tol.

  Example:
  >>> f = lambda x: (x-2)**2
  >>> a = 1
  >>> b = 5
  >>> tol = 1e-5
  >>> (c,d) = gss(f, a, b, tol)
  >>> print(c, d)
  1.9999959837979107 2.0000050911830893
  """
  (a, b) = (min(a, b), max(a, b))
  h = b - a
  if h <= tol:
    return (a, b)

  # Required steps to achieve tolerance
  n = int(ceil(log(tol / h) / log(invphi)))

  c = a + invphi2 * h
  d = a + invphi * h
  yc = f(c)
  yd = f(d)

  for k in range(n-1):
    if yc < yd:
      b = d
      d = c
      yd = yc
      h = invphi * h
      c = a + invphi2 * h
      yc = f(c)
    else:
      a = c
      c = d
      yc = yd
      h = invphi * h
      d = a + invphi * h
      yd = f(d)

  if yc < yd:
      return (a, d)
  else:
      return (c, b)


print()

seed = int(time.time()*1000)
print("Seed: ",seed)
print()

random.seed(int(seed))

n = 1000

p1 = random.uniform(0,pi/2)
p2 = random.uniform(0,pi/2)

p1 = 1e-32
p2 = 0.25

v1 = 0.0
v2 = 0.0

m = 0.0

outer = 0
inner = 0

while True:
  print("m: %17.10le" % m)
  print()

  f3 = func3(0.0)
  print("Point: (%20.18f %20.18f)" % (p1,p2))
  print()

  print("Values: %10.8f %10.8f %10.8f" % f3)
  print()

  df3 = (fabs(f3[0]-f3[1]),fabs(f3[1]-f3[2]),fabs(f3[2]-f3[0]))
  print("Diffs:  %17.10e %17.10e %17.10e" % df3)
  print()

  if max(df3) < 1e-12:
    break

  d3 = list(deri3(p1,p2))

  for d in d3:
    print("   D: (%17.10le %17.10le)" % tuple(d))
  print()

  for i in range(3):
    d = list(d3[i])

    ld = sqrt(d[0]*d[0]+d[1]*d[1])

    if ld > 0:
      d[0] = d[0]/ld
      d[1] = d[1]/ld

    d3[i] = d

  for d in d3:
    print("   D: (%17.10le %17.10le)" % tuple(d))
  print()

  fmax = max(f3)

  v1 = 0.0
  v2 = 0.0

  for i in range(3):
    if fabs(f3[i] - fmax) < 1e-12:
      v1 -= d3[i][0] 
      v2 -= d3[i][1] 

  lv = sqrt(v1*v1 + v2*v2)
  print("Length: %17.10le" % lv)
  print()

  if lv > 0.0:
    v1 = v1/lv
    v2 = v2/lv

  print("Vector: (%17.10le %17.10le)" % (v1,v2))
  print()

  v1 = v1/n
  v2 = v2/n

  lastv = func(-2.0)
  curv = func(-1.0) 

  i = 0.0;
  while True:
    nextv = func(i)

    #print("%6d %10.8f %10.8f %10.8f (%10.8f %10.8f)" % (i,lastv,curv,nextv,p1+i*v1,p2+i*v2))

    if lastv >= curv and curv <= nextv:
      print("Min: %4d %10.8f %10.8f %10.8f" % (i,lastv,curv,nextv))
      print()

      (m1,m2) = gss(func,i-2,i,1e-16)

      m = (m1+m2)/2

      mt1 = m*v1
      mt1 = mt1 - floor(mt1/pi)*pi

      mt2 = m*v2
      mt2 = mt2 - floor(mt2/pi)*pi

      break

    lastv = curv
    curv  = nextv

    i += 1

    inner += 1

  p1 = p1 + m*v1
  p2 = p2 + m*v2

  outer += 1

print("Values:")
print("   %20.18f" % f3[0])
print("   %20.18f" % f3[1])
print("   %20.18f" % f3[2])
print()

print("Outer iterations: %d" % outer)
print("Inner iterations: %d" % inner)
print()
