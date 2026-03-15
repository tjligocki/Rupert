#!/usr/bin/env python3

import sys

for line in sys.stdin:
  parts = line.split()

  if len(parts) > 0:
    m = int(parts[0])
    n = int(parts[1])
    s = float(parts[2])

    print("%4d %4d %23.18f %23.18f" % (m,n,s,s**2))
  else:
    print()
