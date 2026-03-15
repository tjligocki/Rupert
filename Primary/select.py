#!/usr/bin/env python3

import sys

if len(sys.argv) != 3:
  sys.stderr.write("Usage: select.py m r\n")
  sys.exit(1)

m0 = int(sys.argv[1])
r0 = int(sys.argv[2])

for line in sys.stdin:
  parts = line.split()

  if len(parts) > 0:
    m = int(parts[0])
    n = int(parts[1])

    if m == m0 and n % m == r0:
      sys.stdout.write(line)

sys.exit(0)
