#!/usr/bin/env python

import sys,os

from scipy.optimize import *
import numpy as np

from math import *

def rupertMatrix(x,m,n):
  globalRotation = np.eye(n)

  index = 0

  for i in xrange(0,m):
    for j in xrange(i+1,n):
      localRotation = np.eye(n)

      localRotation[i,i] = cos(x[index])
      localRotation[j,j] = cos(x[index])

      localRotation[i,j] = -sin(x[index])
      localRotation[j,i] =  sin(x[index])

      globalRotation = np.matmul(localRotation,globalRotation)

      index += 1

  return globalRotation

def rupert(x,m,n):
  globalRotation = np.eye(n)

  index = 0

  for i in xrange(0,m):
    for j in xrange(i+1,n):
      localRotation = np.eye(n)

      localRotation[i,i] = cos(x[index])
      localRotation[j,j] = cos(x[index])

      localRotation[i,j] = -sin(x[index])
      localRotation[j,i] =  sin(x[index])

      globalRotation = np.matmul(localRotation,globalRotation)

      index += 1

  maxRow = 0.0

  for i in xrange (0,n):
    curRow = 0.0
    for j in xrange (0,m):
      curRow += fabs(globalRotation[i,j])

    maxRow = max(curRow,maxRow)

  return maxRow

m = 1
n = 7

ranges = ((0,pi/2),)*((m*(2*n-(m+1)))/2)
problem = (m,n)
res = brute(rupert,ranges,problem,disp=True)

print res
print 
print rupertMatrix(res,m,n)
print
print 1.0/rupert(res,m,n)
