#!/usr/bin/env python3

import math,sys

def main(argc,argv):
  x = 0.1
  if (argc >= 2):
    x = float(argv[1])

  alpha = 0.01
  if (argc >= 3):
    alpha = float(argv[2])

  absx = ABS(alpha,x)

  print(x,"->",absx)


def ABS(alpha,x):
  return math.sqrt(x*x+alpha*alpha)


def MAX(alpha,x,y):
  return 0.5*(x + y - ABS(alpha,x-y))


if __name__ == "__main__":
  main(len(sys.argv),sys.argv)
