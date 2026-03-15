import sys

lines = sys.stdin.readlines()

num_lines = len(lines)

n = num_lines - 4

for i in xrange(1,n+1):
  columns = lines[i].split()
  for c in columns[1:-1]:
    print "%9.6f" % float(c),
  print
