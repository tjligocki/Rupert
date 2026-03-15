#!/usr/bin/env python3

import sys,math
import numpy as np
import rupert

def main(argc,argv):
  m = 2
  if argc > 1:
    m = int(argv[1])

  n = 3
  if argc > 2:
    n = int(argv[2])

  rng = np.random.default_rng()

  total = rupert.total(m,n)

  print("%d, %d -> %d" % (m,n,total))
  print()

  arot = math.pi * (2.0*rng.random(total) - 1.0)

  print(arot)
  print()

  numpts = 1000
  delta  = 0.001

  count = 0

  while count < 1000 and delta > 1e-15:
    darot = 2.0*rng.random(total) - 1.0

    norm = np.linalg.norm(darot)
    if norm != 0.0:
      darot /= norm

    rmin = 0.0
    for i in range(-numpts//2,numpts//2 + 1):
      brot = arot + i*delta * darot

      rup = rupert.rupert(m,n,brot)

      r = rup.value()

      if i == -numpts//2 or r < rmin:
        smin = i*delta
        rmin = r

        bminrot = brot.copy()

    if smin != 0.0:
      print("%24.21f (%e) %24.21f (%24.21f), %d" % (smin,delta,rmin,1.0/rmin,count));

      if delta < 1.0:
        delta *= 1.2

      count = 0
    else:
      count += 1

    arot = bminrot.copy()

    delta *= 0.99

  print("%19.16f (%19.16f), %23.16e, " % (rmin,1.0/rmin,delta))


if __name__ == "__main__":
  main(len(sys.argv),sys.argv)
