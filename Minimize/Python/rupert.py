import sys,math
import numpy as np

class rupert:
  def __init__(self,m,n,rot,eps):
    self.m = m
    self.n = n

    self.eps = eps

    self.mat = np.identity(self.n)

    tot = total(self.m,self.n)

    gradMats = np.zeros([tot,n,n])
    for i in range(tot):
      gradMats[i] = np.identity(self.n)

    irot = 0

    #print('  rot:')
    #print(rot)
    #print()

    for i in range(0,self.m):
      for j in range(i+1,self.n):
        t = rot[irot]

        s = math.sin(t)
        c = math.cos(t)

        rotmat = np.identity(self.n)

        rotmat[i][i] =  c
        rotmat[i][j] = -s
        rotmat[j][i] =  s
        rotmat[j][j] =  c

        #print('    rotmat:',rotmat)

        self.mat = np.matmul(self.mat,rotmat)

        drotmat = np.zeros((self.n,self.n))

        drotmat[i][i] = -s
        drotmat[i][j] = -c
        drotmat[j][i] =  c
        drotmat[j][j] = -s

        for k in range(tot):
          if k == irot:
            gradMats[k] = np.matmul(gradMats[k],drotmat)
          else:
            gradMats[k] = np.matmul(gradMats[k], rotmat)

        irot += 1

    self.mat = self.mat[:,0:m]
    #print('self.mat:',np.shape(self.mat))
    #print(self.mat)
    #print()

    absMat = absSqrt(self.mat,self.eps)
    #print('absMat:',np.shape(absMat))
    #print(absMat)
    #print()

    sumAbsRows = absMat.sum(axis=1)
    #print('sumAbsRows:',np.shape(sumAbsRows))
    #print(sumAbsRows)
    #print()

    self.value = maxExp(sumAbsRows,self.eps)
    #print('self.value:')
    #print(self.value)
    #print()

    dAbsSqrtMat = dAbsSqrt(self.mat,self.eps)
    #print('dAbsSqrtMat:',np.shape(dAbsSqrtMat))
    #print(dAbsSqrtMat)
    #print()

    dMaxExpRows = dMaxExp(sumAbsRows,self.eps)
    #print('dMaxExpRows:',np.shape(dMaxExpRows))
    #print(dMaxExpRows)
    #print()

    self.grad = np.zeros(tot)
    for k in range(tot):
      #print('gradMats[%d]:' % k,np.shape(gradMats[k][:,0:m]))
      #print(gradMats[k][:,0:m])
      #print()

      self.grad[k] = (dMaxExpRows * (dAbsSqrtMat * gradMats[k][:,0:m]).sum(axis=1)).sum()

    #print('self.grad:',np.shape(self.grad))
    #print(self.grad)
    #print()

    #imax = 0
    #self.grad = (deriv[imax] * gradMats[:,imax,0:m]).sum(axis=1)
    #print('self.grad:',np.shape(self.mat))
    #print(self.grad)
    #print()


def total(m,n):
  return (m * (2*n - m - 1)) // 2


def absSqrt(x,eps):
  return np.sqrt(x**2 + eps**2)


def dAbsSqrt(x,eps):
  return x/np.sqrt(x**2 + eps**2)


def maxExp(x,eps):
  xm = np.max(x)

  max = eps * np.log(np.exp((x-xm)/eps).sum(0)) + xm

  return max


def dMaxExp(x,eps):
  xm = np.max(x)

  max = np.exp((x-xm)/eps)/(np.exp((x-xm)/eps).sum(0))

  return max
