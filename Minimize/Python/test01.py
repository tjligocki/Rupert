#!/usr/bin/env python3

import sys,math
import numpy as np
import rupert

def main(argc,argv):
  mat1 = np.zeros(9).reshape(3,3)
  mat2 = np.zeros(9).reshape(3,3)
  mat3 = np.zeros(9).reshape(3,3)
  mat4 = np.zeros(9).reshape(3,3)
  mat5 = np.zeros(9).reshape(3,3)

  print("Starting to write Python code to work on the Generalized Rupert Problem.")
  print()
  
  a01 = math.pi/4
  a02 = 0.0
  a12 = 4 * math.atan(math.sqrt(5.0 - 2.0*math.sqrt(6.0)))

  print("a01: %9.6f (%9.4f)" % (a01,math.degrees(a01)))
  print("a02: %9.6f (%9.4f)" % (a02,math.degrees(a02)))
  print("a12: %9.6f (%9.4f)" % (a12,math.degrees(a12)))
  print()

  rot = [a01,a02,a12]

  rup = rupert.rupert(2,3,rot)

  print(rup.mat)
  print()

  print(np.linalg.det(rup.mat))
  print()

  print(rup.value())
  print()


if __name__ == "__main__":
  main(len(sys.argv),sys.argv)
