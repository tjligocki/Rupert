#!/usr/bin/env python3

import math,sys

def main(argc,argv):
  alpha = 0.01
  if (argc >= 2):
    alpha = float(argv[1])

  vars = map(float,argv[2:])

  max_vars = LSE(alpha,vars)


def MAX(alpha,vars):
  inv_alpha = 1.0/alpha

  sum = 0.0
  for v in vars:
    sum += math.exp(inv_alpha*float(v))

  aver = alpha * math.log(sum) 

  return aver


def ABS(alpha,x):
  return 2*MAX(alpha,(x,0)) - x


if __name__ == "__main__":
  main(len(sys.argv),sys.argv)
