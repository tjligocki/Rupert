#!/usr/bin/env python

import sys

print "                                                                  f(m-1,m) - 1"
print "  m    n        f(m,m+1) - 1                  k               9 - ------------"
print "                                                                  f(m,m+1) - 1"
print " ---  ---  ----------------------  ----------------------  ----------------------"

fm1_int_prev  = -1

for line in (sys.stdin):
  parts = line.split()
  m = int(parts[0])
  n = int(parts[1])
  f_str = parts[2]

  fm1_str = "0" + f_str[1:]

  f   = float(f_str)
  fm1 = float(fm1_str)
  fm1_int = int(fm1_str[2:])

  print "%4d %4d %23.16le" % (m,n,fm1),
  print "%23.16le" % ((9**m * fm1_int)/float(10**200)),

  if fm1_int_prev > 0:
    print "%23.16le" % ((9*fm1_int - fm1_int_prev)/float(fm1_int))
  else:
    print

  fm1_int_prev  = fm1_int
