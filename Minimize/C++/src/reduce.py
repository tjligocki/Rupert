#!/usr/bin/env python3

import sys

def main(argc,argv):
  m = 2
  if (argc >= 2):
    m = int(argv[1])

  n = 3
  if (argc >= 3):
    n = int(argv[1])

  for line in sys.stdin:
    line = line.replace("("," ").replace(","," ").replace(")"," ")
    parts = line.split()

    parts = map(float,parts)

    for part in parts:
      print("%13.8lf" % (part,),end = " ")
    print()


if __name__ == "__main__":
  main(len(sys.argv),sys.argv)
