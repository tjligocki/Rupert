import sys,math
import numpy as np

class rupert:
  def __init__(self,rot):
    self.m = 2
    self.n = 3
    self.mat = np.identity(self.n)

    tot = total(m,n)

    self.grad = tot*[0]
    for i in range(tot):
      self.grad[i] = np.identity(self.n)

    irot = 0

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

        self.mat = np.matmul(self.mat,rotmat)

        for k in range(tot):
          if k = irot:
            rotmat[i][i] = -s
            rotmat[i][j] = -c
            rotmat[j][i] =  c
            rotmat[j][j] = -s

          self.grad[k] = np.matmul(self.grad[k],rotmat)

        irot += 1

  def value(self):
    return abs(self.mat[:,0:self.m]).sum(axis=1).max()

  def value_approx(self,eps):
    return vec_max_approx(mat_abs_approx(self.mat[:,0:self.m],eps).sum(axis=1),eps)


def total(m,n):
  return (m * (2*n - m - 1)) // 2


def mat_abs_approx(x,eps):
  return np.sqrt(x**2 + eps**2)


def abs_approx(x,eps):
  return math.sqrt(x*x + eps*eps)


def mat_max_approx(x,eps):
  return np.apply_along_axis(vec_max_approx,0,x,eps)


def vec_max_approx(x,eps):
  max = x[0]

  for v in x[1:]:
    max = max_approx(max,v,eps)

  return max


def max_approx(x,y,eps):
  return (x + y + abs_approx(x - y,eps)) / 2.0
