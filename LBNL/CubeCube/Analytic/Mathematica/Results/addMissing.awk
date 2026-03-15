BEGIN {
  curN = 0;
  prevM = 0;

  countMissing = 0;
}
NF == 0 {
  print;

  total = n-1;
  countPresent = total - countMissing;
  printf("%3d missing (%3d%%), %3d present (%3d%%) - %3d total\n",countMissing,(100.0*countMissing)/total,countPresent,(100.0*countPresent)/total,total);

  curN = 0;
  prevM = 0;

  countMissing = 0;

  print;
  print;

  next;
}
{
  m = $1;
  n = $2;
  value = $3;

  if (n != curN) {
    curN = n;
    prevM = 0;
  }

  for (pm = prevM+1; pm < m; pm++) {
    countMissing++;
    printf("%3d %3d\n",pm,n);
  }
  printf("%3d %3d    %s\n",m,n,value);

  prevM = m;

  next;
}
