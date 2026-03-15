from math import *
import sys,copy

def radians(d):
  return d*pi/180.0

def rotate(iv,t,i1,i2):
  ov = []
  for ip in iv:
    op = copy.copy(ip)
    op[i1] = cos(t)*ip[i1] - sin(t)*ip[i2]
    op[i2] = sin(t)*ip[i1] + cos(t)*ip[i2]
    ov.append(op)

  return ov

v = [[ 0.00, 0.00, 0.00],
     [ 1.00, 0.00, 0.00],
     [ 0.00,-1.00, 0.00],
     [ 0.00, 0.00, 1.00],
     [ 1.00,-1.00, 0.00],
     [ 0.00,-1.00, 1.00],
     [ 1.00, 0.00, 1.00],
     [ 1.00,-1.00, 1.00],
     [ 0.25,-1.00, 0.00],
     [ 0.00,-0.75, 1.00],
     [ 1.00,-0.25, 0.00],
     [ 0.75, 0.00, 1.00]]

s = 3.0

ox = 1.5
oy = 2.0

v = rotate(v,radians(90),1,2)
v = rotate(v,radians(-70),0,2)
v = rotate(v,radians(30),1,2)

es = [[ 1, 2],[ 1, 3],[ 1, 4],
      [ 2, 5],[ 2, 7],
      [ 3, 5],[ 3, 6],
      [ 4, 6],[ 4, 7],
      [ 5, 8],[ 6, 8],[ 7, 8],
      [ 9,10],[ 9,11],[10,12],[11,12]]

fv = v
for p in fv:
  p[0] += ox
  p[1] += oy

print "#FIG 3.2  Produced by xfig version 3.2.5b"
print "Landscape"
print "Center"
print "Inches"
print "Letter" 
print "100.00"
print "Single"
print "-2"
print "1200 2"

pt = 2400
i = 1
for e in es:
  v1 = e[0]-1
  v2 = e[1]-1
  if i <= 12:
    print "2 1 0 2  0 7 50 -1 -1 1.000 0 0 -1 0 0 2"
  else:
    print "2 1 1 2 20 7 55 -1 -1 6.000 0 0 -1 0 0 2"
  print "  %5d %5d %5d %5d" % (int(pt*fv[v1][0]+0.5),int(pt*fv[v1][1]+0.5),int(pt*fv[v2][0]+0.5),int(pt*fv[v2][1]+0.5))

  i += 1
