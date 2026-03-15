import math,random,time,sys
import numpy


class Rupert:
  def __init__(self,m,n):
    self.m = m
    self.n = n

    self.matrix    = numpy.identity(n)
    self.cur_value = self.value()

    self.limit = math.pi

  def __mul__(self,other):
    prod = Rupert(self.m,self.n)
    prod.matrix = numpy.dot(self.matrix,other.matrix)
    return prod

  def value(self):
    return abs(self.matrix[:,0:self.m]).sum(axis=1).max()

  def new_opt(self,N=1000,num_levels=1):
    rotation1 = Random_Rotation(self.n)
    rotation2 = Random_Rotation(self.n)

    # min_value = self.value()
    min_value = self.cur_value
    min_matrix = self

    start1 = -self.limit
    end1   =  self.limit

    start2 = -self.limit
    end2   =  self.limit

    min_t1 = 0.0
    min_t2 = 0.0

    for level in xrange(num_levels):
      dt1 = (end1 - start1) / (N-1)
      dt2 = (end2 - start2) / (N-1)

      for it1 in xrange(N):
        t1 = it1 * dt1 + start1

        rotation1.rotate_matrix(t1)

        new_matrix = self * rotation1

        for it2 in xrange(N):
          t2 = it2 * dt2 + start2

          rotation2.rotate_matrix(t2)

          new_new_matrix = new_matrix * rotation2
          value = new_new_matrix.value()

          if value <= min_value:
            min_t1     = t1
            min_t2     = t2
            min_value  = value
            min_matrix = new_new_matrix

          value = new_matrix.value()

      start1 = min_t1 - dt1
      end11  = min_t1 + dt1

      start2 = min_t2 - dt2
      end12  = min_t2 + dt2

    if min_value != self.cur_value:
      self.limit = (N**num_levels) * abs(max(min_t1,min_t2))
      pass
    else:
      self.limit = self.limit / 1.05
      if self.limit < 1e-16:
        self.limit = math.pi

    self.limit = min(self.limit,math.pi)

    self.matrix    = min_matrix.matrix
    self.cur_value = min_value

    return min_value

  def new_random(self):
    rotation = Random_Rotation(self.n)

    t = random.uniform(-math.pi,math.pi)

    rotation.rotate_matrix(t)

    new_matrix = self * rotation
    value = new_matrix.value()

    self.matrix    = new_matrix.matrix
    self.cur_value = value


class Random_Unit_Vector:
  def __init__(self,n):
    self.n = n

  def unit_vector(self):
    norm = 0.0
    while norm == 0.0:
      vect = numpy.array([random.gauss(0.0,1.0) for i in xrange(self.n)])
      norm = math.sqrt(numpy.dot(vect,vect.conj()))
    vect /= norm
    return vect


class Random_Rotation:
  def __init__(self,n):
    self.n = n

    self.vector1 = Random_Unit_Vector(self.n).unit_vector()

    norm = 0.0
    while norm == 0.0:
      vector2 = Random_Unit_Vector(self.n).unit_vector()
      vector2 = vector2 - numpy.dot(self.vector1,vector2.conj())*self.vector1
      norm = math.sqrt(numpy.dot(vector2,vector2.conj()))

    self.vector2 = vector2/norm

    self.t = 0.0
    self.rotate_matrix(self.t)

  def rotate_matrix(self,t):
    self.t = t

    self.matrix = numpy.identity(self.n)

    ct = math.cos(t)
    st = math.sin(t)

    for i in xrange(self.n):
      for j in xrange(self.n):
        self.matrix[i,j] += (self.vector1[i]*self.vector1[j] + self.vector2[i]*self.vector2[j]) * (ct-1) + (self.vector1[j]*self.vector2[i] - self.vector1[i]*self.vector2[j]) * st


if __name__ == "__main__":
  random.seed(time.time())

  m = 2
  n = 3

  N = 20
  num_levels = 1

  tries  = 2000
  report =  100

  argc = len(sys.argv)

  if argc >= 3:
    m = int(sys.argv[1])
    n = int(sys.argv[2])

  if argc >= 5:
    N          = int(sys.argv[3])
    num_levels = int(sys.argv[4])

  if argc >= 6:
    tries = int(sys.argv[5])

  if argc >= 7:
    report = int(sys.argv[6])

  rm = Rupert(m,n)

  rm.new_random()
  #rm.new_random()
  #rm.new_random()
  #rm.new_random()

  for i in xrange(tries):
    if i % report == 0:
      print "%7d  %30.23le  %30.23le" % (i,rm.limit,1.0/rm.cur_value)

    rm.new_opt(N,num_levels)

  print

  print "Best (%d,%d): %30.23le" % (m,n,1.0/rm.cur_value)
