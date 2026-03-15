from pprint import pprint
import operator,math

def trans_p(p,t):
  pout = []
  for c in p:
    cout = []
    for r in t:
      cout.append(reduce(operator.add,[e[0]*e[1] for e in zip(r,c)]))
    pout.append(cout)

  return pout

def v_from_p(p0,p1):
  return [e[0] - e[1] for e in zip(p0,p1)]

def len_v(v):
  return math.sqrt(reduce(operator.add,[c*c for c in v]))

def scale_v(a,v):
  return [a*c for c in v]

def add_v(v0,v1):
  return [e[0] + e[1] for e in zip(v0,v1)]

def sub_v(v0,v1):
  return [e[0] - e[1] for e in zip(v0,v1)]

m = 2
n = 4

p0 = []
p0.append([-0.50,-0.50,0.00,0.00])
p0.append([ 0.50,-0.50,0.00,0.00])
p0.append([-0.50, 0.50,0.00,0.00])
p0.append([ 0.50, 0.50,0.00,0.00])

t = []
t.append([ 0.707107, 0.000000,-0.361083, 0.607963])
t.append([ 0.707107, 0.000000, 0.361083,-0.607963])
t.append([ 0.000000, 0.707107,-0.607963,-0.361083])
t.append([ 0.000000, 0.707107, 0.607963, 0.361083])

p1 = trans_p(p0,t)
pm = max(max(p1))

p2 = []
for c1 in p1:
  c2 = scale_v(0.5/pm,c1)
  p2.append(c2)

p = p2

v1 = v_from_p(p[1],p[0])
v2 = v_from_p(p[2],p[0])

print 'set xrange [-0.5:1.5]'
print 'set yrange [-0.5:1.5]'
print 'set trange [-0.5:1.5]'

print 'set size ratio -1'

print 'set parametric'

plot_str = 'plot'
for i in xrange(0,4):
  if v2[i] != 0:
    m = -v1[i]/v2[i]
    b = (-0.5-p[0][i])/v2[i]
    if b >= 0:
      print '%s t,%f*t+%f lw 2 lc 1' % (plot_str,m,b)
    else:
      print '%s t,%f*t-%f lw 2 lc 1' % (plot_str,m,-b)
  else:
    b = (-0.5-p[0][i])/v1[i]
    print '%s %f,t lw 2 lc 1' % (plot_str,b)

  plot_str = 'replot'

  if v2[i] != 0:
    m = -v1[i]/v2[i]
    b = (0.5-p[0][i])/v2[i]
    if b >= 0:
      print '%s t,%f*t+%f lw 2 lc 1' % (plot_str,m,b)
    else:
      print '%s t,%f*t-%f lw 2 lc 1' % (plot_str,m,-b)
  else:
    b = (0.5-p[0][i])/v1[i]
    print '%s %f,t lw 2 lc 1' % (plot_str,b)

print 'set style arrow 1 nohead lw 2 front'

print 'set arrow from 0,0 to 1,0 as 1'
print 'set arrow from 1,0 to 1,1 as 1'
print 'set arrow from 1,1 to 0,1 as 1'
print 'set arrow from 0,1 to 0,0 as 1'
print 'replot NaN,NaN' 

