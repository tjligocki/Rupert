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

print()
print(myseed)

random.seed(myseed)

t = [random.uniform(0,pi/2) for i in range(n-1)]
tinit = t[:]

m = 1.0
count = 0

var = 1e-1
ccount = 0

while m > 1e-15 and count < 500000:
  count += 1

  fl = funcs.fa(t)

  min_index, min_value = min(enumerate(fl), key=operator.itemgetter(1))
  max_index, max_value = max(enumerate(fl), key=operator.itemgetter(1))

  m = (max_value - min_value)
  aver_value = sum(fl)/len(fl)

  while m < var:
    print("%8.2e %6d (%6d)" % (var,count,count-ccount))
    var /= 10
    ccount = count

  dt = [0.0]*(n-1)

  ci = random.randint(0,n-1)
  for i in [ci]:
    norm = 0.0
    for j in range(n-1):
      norm += funcs.df(i,j,t)**2
    norm = sqrt(norm)

    if norm == 0.0:
      norm = 1.0

    scale = random.uniform(0.1,10)

    for j in range(n-1):
      dt[j] += (aver_value-funcs.f(i,t))*funcs.df(i,j,t)/(scale*norm)

  """
  print()
  print("  ",t)
  print("  ",fl)
  print("  ",ci)
  print()
  print("  ",end="")
  for i in range(n):
    print("  ",end="")
    for j in range(n-1):
      print(funcs.df(i,j,t),end=", ")
    print()
  print()
  for i in range(n):
    print("  ",aver_value-f(i,t))
  print()
  print("  ",min_index,min_value)
  print("  ",max_index,max_value)
  print()
  print("  ",m)
  print()
  print("  ",dt)
  print()
  """

  t = [t[i]+dt[i] for i in range(n-1)]

print()
print(tinit)
print(funcs.fa(tinit))
print(1.0/funcs.ff(tinit))
print()
print(t)
print(funcs.fa(t))
print(1.0/funcs.ff(t))
print()
print(m)
print(count)
print()
