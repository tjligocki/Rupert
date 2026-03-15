BEGIN {
  curM = 0;
  prevN = 0;

  countMissing = 0;
}
NF == 0 {
  print;

  total = 100-m;
  countPresent = total - countMissing;
  printf("%3d missing (%3d%%), %3d present (%3d%%) - %3d total\n",countMissing,(100.0*countMissing)/total,countPresent,(100.0*countPresent)/total,total);

  curM = 0;
  prevN = 0;

  countMissing = 0;

  print;
  print;

  next;
}
{
  m = $1;
  n = $2;
  value = $3;

  if (m != curM) {
    curM = m;
    prevN = n;
  }

  for (pn = prevN+1; pn < n; pn++) {
    countMissing++;
    printf("%3d %3d\n",m,pn);
  }
  printf("%3d %3d    %s\n",m,n,value);

  prevN = n;

  next;
}
