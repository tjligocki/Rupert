#!/usr/bin/env python3

import sys,time


def main(argc,argv):
  if argc <= 2:
    sys.stderr.write("Usage: results_compare reference_file result_file")
    sys.exit(1)

  ref_filename = argv[1]
  res_filename = argv[2]

  ref_file = open(ref_filename,'r')
  res_file = open(res_filename,'r')

  res_dict = {}
  for res_line in res_file:
    res_parts = res_line.split()

    res_m = int(res_parts[0])
    res_n = int(res_parts[1])

    res_dim = int(res_parts[2])

    res_val = float(res_parts[4])

    res_dict[(res_m,res_n)] = [res_dim,res_val,False]

  for ref_line in ref_file:
    ref_parts = ref_line.split()

    if len(ref_parts) >= 3:
      ref_m = int(ref_parts[0])
      ref_n = int(ref_parts[1])

      ref_val = float(ref_parts[2])

      if (ref_m,ref_n) in res_dict:
        res_list = res_dict.get((ref_m,ref_n))
        res_list[2] = ref_val
        res_dict[(ref_m,ref_n)] = res_list

  for e in res_dict:
    l = res_dict.get(e)

    if l[2] != False:
      print('%3d %3d  %21.16f  %21.16f  %11.4e' % (e[0],e[1],l[2],l[1],l[1] - l[2]))
    else:
      print(e[0],e[1],'case not found in reference')


if __name__ == '__main__':
  main(len(sys.argv),sys.argv)

  sys.exit(0)
