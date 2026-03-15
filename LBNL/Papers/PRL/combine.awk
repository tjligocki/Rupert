BEGIN {
  num = 0;

  maxM = 0;
  maxN = 0;
}
/^$/ {
  next;
}
{
  m[num] = $1;
  n[num] = $2;
  v[num] = $3;

  if (m[num] > maxM) {
    maxM = m[num];
  }

  if (n[num] > maxN) {
    maxN = n[num];
  }

  num += 1;

  next;
}
END {
  for (i = 0; i < num; i++)
  {
    h = m[i] "," n[i];

    if (g[h] == 0) {
      g[h] = v[i];
    }
    else
    {
      use = v[i];
      mj = length(g[h]);

      if (length(v[i]) < mj) {
        use = g[h];
        mj = length(v[i]);
      }

      for (j = 1; j <= mj; j++)
      {
        curG = substr(g[h],j,1);
        curV = substr(v[i],j,1);

        if (curG > curV) {
          use = g[h];
          break;
        }
        else if (curG < curV) {
          use = v[i];
          break;
        }
      }

      g[h] = use;
    }
  }

  for (mm = 1; mm <= maxM; mm++) {
    for (nn = 1; nn <= maxN; nn++) {
      h = mm "," nn;
      if (g[h] != 0) {
        printf("%3d %3d    %s\n",mm,nn,g[h]);
      }
    }

    print;
  }
}
