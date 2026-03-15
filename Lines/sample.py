#!/usr/bin/env python3

import rupertFunction

f = rupertFunction.rupertFunction(2,3)

deg1 = 180
nsteps1 = 90
dt1 = deg1/nsteps1

deg2 = 360
nsteps2 = 180
dt2 = deg2/nsteps2

deg3 = 0
nsteps3 = 1
dt3 = deg3/nsteps3

for it1 in range(nsteps1+1):
  t1 = it1*dt1
  for it2 in range(nsteps2+1):
    t2 = it2*dt2
    for it3 in range(nsteps3):
      t3 = it3*dt3
      f.setMatrix((t1,t2,t3))
      print(t1,t2,t3,f.f())
  print()
