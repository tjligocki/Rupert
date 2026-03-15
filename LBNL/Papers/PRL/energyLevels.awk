BEGIN {
  maxM = 0;
  maxN = 0;
}
/^$/ {
  next;
}
{
  m = $1;
  n = $2;
  v = $3;

  if (m > maxM) {
    maxM = m;
  }

  if (n > maxN) {
    maxN = n;
  }

  table[m "," n] = v;

  next;
}
END {
  for (m = 1; m <= maxM; m++) {
    for (mod = 0; mod < m; mod++) {
      file = sprintf("Split/m%03dmod%03d",m,mod);
      for (n = m+mod; n <= maxN; n += m) {
        if (table[m "," n] > 0) {
          printf("%3d %3d    %s\n",m,n,table[m "," n]) >> file;
        }
      }
    }
  }
}
