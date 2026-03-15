#!/usr/bin/env python3

import sys,LSE

def main(argc,argv):
  alpha = 0.01
  if (argc >= 2):
    alpha = float(argv[1])

  xmin = -1.0
  xmax =  1.0
  if (argc >= 4):
    xmin = float(argv[2])
    xmax = float(argv[3])

  nx = 100
  if (argc >= 5):
    nx = int(argv[4])

  scan(alpha,xmin,xmax,nx)


def scan(alpha,xmin,xmax,nx):
  dx = (xmax-xmin) / nx

  for i in range(0,nx+1):
    x = i*dx + xmin

    absx = LSE.ABS(alpha,x)

    print(x,absx)


if __name__ == "__main__":
  main(len(sys.argv),sys.argv)
