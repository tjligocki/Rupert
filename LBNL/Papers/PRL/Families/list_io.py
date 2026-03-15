def read_key(file):
  list_key = {}
  for line in file:
    parts = line.split()
    if len(parts) == 0:
      break
    else:
      list_key[parts[0]] = parts[1:]

  return list_key

def write_key(list_key,file):
  list_indices = list_key.keys()
  list_indices.sort()

  for index in list_indices:
    file.write("%s " % index)

    for item in list_key[index]:
      file.write("%s " % item)

    file.write("\n")

  file.write("\n")

def read_body(file):
  list_pairs = {}

  maxm = 0
  maxn = 0

  for line in file:
    parts = line.split()

    if len(parts) > 0:
      m = int(parts[0])
      n = int(parts[1])

      if m > maxm:
        maxm = m

      if n > maxn:
        maxn = n

      info = []
      if len(parts) > 2:
        info = parts[2].split(",")

      list_pairs[(m,n)] = set(info)

  return (maxm,maxn,list_pairs)

def set_to_str(in_set):
  items = sorted(in_set)

  out_str = ""
  nitems = len(items)
  if nitems > 0:
    for n in xrange(nitems-1):
      out_str += "%s," % items[n]
    out_str += "%s" % items[nitems-1]

  return out_str


def write_body(maxm,maxn,list_pairs,file):
  for n in xrange(1,maxn+1):
    for m in xrange(1,min(n+1,maxm+1)):
      if (m,n) in list_pairs:
        file.write("%3d %3d   %s\n" % (m,n,set_to_str(list_pairs[(m,n)])))


    file.write("\n")
