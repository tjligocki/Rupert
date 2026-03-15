#!/usr/bin/env python3

import sys,operator,random,time,importlib
from math import *

sign = lambda x: -1 if x < 0 else (1 if x > 0 else 0)

m = 2
n = 3

if len(sys.argv) >= 3:
  m = int(sys.argv[1])
  n = int(sys.argv[2])

funcfile = "f" + ("_%02d" % m) + ("_%02d" % n) 

funcs = importlib.import_module(funcfile)

myseed = int(time.time()*10000*time.time()) % 1000000000
#myseed = 157694464 # bad (2,3)
#myseed =  83624960 # bad (2,5)
#myseed = 746582528 # bad (2,11)
#myseed = 816931328 # bad (3,4)

print()
print(myseed)

random.seed(myseed)

t = [random.uniform(0,pi/2) for i in range(n-1)]
tinit = t[:]

m = 1.0
count = 0

var = 1e-1
ccount = 0

while m > 1e-14 and count < 50000:
  count += 1

  fl = funcs.fa(t)

  min_index, min_value = min(enumerate(fl), key=operator.itemgetter(1))
  max_index, max_value = max(enumerate(fl), key=operator.itemgetter(1))

  m = (max_value - min_value)/2

  while 2*m < var:
    print("%8.2e %6d (%6d)" % (var,count,count-ccount))
    var /= 10
    ccount = count

  dt = [-funcs.df(max_index,i,t)+funcs.df(min_index,i,t) for i in range(n-1)]

  """
  print()
  print("  ",t)
  print("  ",fl)
  print()
  for i in range(n):
    print("  ",end="")
    for j in range(n-1):
      print(funcs.df(i,j,t),end=", ")
    print()
  print()
  print("  ",min_index,min_value)
  print("  ",max_index,max_value)
  print()
  print("  ",m)
  print()
  print("  ",dt)
  print()
  """

  t = [t[i]+m*dt[i] for i in range(n-1)]

print()
print(tinit)
print(funcs.fa(tinit))
print(1.0/funcs.ff(tinit))
print()
print(t)
print(funcs.fa(t))
print(1.0/funcs.ff(t))
print()
print(2*m)
print(count)
print()
