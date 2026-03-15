# Represent a given Rupert function for (m,n)

import math,numpy

class rupertFunction:

  def __init__(self,m,n):
    self.m = m
    self.n = n
    self.initMatrix()
    return

  def initMatrix(self):
    self.rot = numpy.identity(self.n)
    return

  def rotMatrix(self,i,j,theta):
    t = math.pi*theta/180.0

    s = math.sin(t)
    c = math.cos(t)
    mat = numpy.identity(self.n)
    mat[i][i] =  c
    mat[i][j] = -s
    mat[j][i] =  s
    mat[j][j] =  c

    self.rot = numpy.matmul(self.rot,mat)
    return

  def setMatrix(self,thetas):
    self.initMatrix()
    count = 0
    for i in range(self.n):
      for j in range(i+1,self.n):
        self.rotMatrix(i,j,thetas[count])
        count += 1
        
  def f(self):
    maxSum = 0
    i = 0
    for j in range(self.m):
      maxSum += abs(self.rot[i][j])

    for i in range (1,self.n):
      curSum = 0
      for j in range(self.m):
        curSum += abs(self.rot[i][j])

      if curSum > maxSum:
        maxSum = curSum

    return 1.0/maxSum
