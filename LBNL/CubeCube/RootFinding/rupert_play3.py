#!/usr/bin/env python

import numpy as np
from math import *
import sys,time,random

m = int(sys.argv[1])
n = int(sys.argv[2])

nd = (n*(n-1))/2

random.seed(time.time())

# t = pi/6 * np.ones(nd)
t = random.uniform(0.0,pi) * np.ones((nd,1))

# print 
# print t
# print 

iter = 1

old_dt_max = 100.0

while True:
  print '%d' % iter,
  # print

  total_rot = np.identity(n)

  total_der = []
  for i in xrange(nd):
    total_der.append(np.identity(n))

  d = 0
  for i in xrange(n):
    for j in xrange(i+1,n):
      rot = np.identity(n)
      rot[i,i] =  cos(t[d])
      rot[i,j] =  sin(t[d])
      rot[j,i] = -sin(t[d])
      rot[j,j] =  cos(t[d])

      der = np.identity(n)
      der[i,i] = -sin(t[d])
      der[i,j] =  cos(t[d])
      der[j,i] = -cos(t[d])
      der[j,j] = -sin(t[d])

      # print 'rot',i,j
      # print rot
      # print 'deriv'
      # print der
      # print 

      total_rot = np.dot(rot,total_rot)

      # print 'total_der[0] - before'
      # print total_der[0]
      # print

      for k in xrange(nd):
        if k == d:
          total_der[k] = np.dot(der,total_der[k])
        else:
          total_der[k] = np.dot(rot,total_der[k])

      # print 'total_der[0] - after'
      # print total_der[0]
      # print

      d += 1

  sums = np.zeros((n,1))
  sums_der = np.zeros((n,nd))

  for i in xrange(n):
    sums[i] = 0.0

    for k in xrange(nd):
      sums_der[i,k] = 0.0

    for j in xrange(m):
      sums[i] += abs(total_rot[i,j])

      for k in xrange(nd):
        sums_der[i,k] += total_der[k][i,j]

  # print
  # print total_rot
  # print 

  # for k in xrange(nd):
  #   print k,total_der[k]
  #   print

  # print sums
  # print 

  # print sums_der
  # print

  mf = np.zeros((nd,1))

  for i in xrange(n-1):
    mf[i] = -(sums[i]-sums[i+1])

  mf[n-1] = -total_rot[n-2,0]
  mf[n  ] = -total_rot[n-1,0]
  mf[n+1] = -total_rot[n-1,1]

  fp = np.identity(nd)

  for i in xrange(n-1):
    for j in xrange(nd):
      fp[i,j] = sums_der[i,j]-sums_der[i+1,j]

  for j in xrange(nd):
    fp[n-1,j] = total_der[j][n-2,0]
    fp[n  ,j] = total_der[j][n-1,0]
    fp[n+1,j] = total_der[j][n-1,1]

  # print fp
  # print

  # print mf
  # print

  dt = np.linalg.solve(fp,mf)

  # print
  if abs(dt).max() > 0.10:
    # print '   ',abs(dt).max()
    dt *= 0.10/abs(dt).max()

  if abs(dt).max() < 1e-12:
    # print 1.0/abs(sums).max()
    # print
    break

  print 1.0/abs(sums).max(),
  print abs(dt).max(),

  t += dt

  for i in xrange(nd):
    while t[i] < 0.0:
      t[i] += 2*pi
    while t[i] > 2*pi:
      t[i] -= 2*pi

  print t.min(),t.max()
  # print

  # print dt
  # print
  # print t

  iter += 1

  old_dt_max = abs(dt).max()

print
print
print '%.20lf' % (1.0/abs(sums).max(),)
