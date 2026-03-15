#!/usr/bin/env python3

import rupertFunction

f = rupertFunction.rupertFunction(2,3)

print("(m,n) = (",f.m,",",f.n,")")
print()

print("Initial matrix:")
print(f.rot)
print()
print(f.f())
print()

print("New matrix:")
f.setMatrix((45,60,0))
print(f.rot)
print()
print(f.f())
print()
