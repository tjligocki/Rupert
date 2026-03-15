BEGIN {
  count = 0;
}
NF == 1 {
  curfile = $1;
}
NF == 3 {
  if (max[$1 "," $2] == 0) {
    key[count] = $1 "," $2;
    count++;
  }
  if ($3 > max[$1 "," $2]) {
    max[$1 "," $2] = $3;
    file[$1 "," $2] = curfile;
    p[$1 "," $2] = $1;
    q[$1 "," $2] = $2;
  }
}
END {
  for (i = 0; i < count; i++) {
    k = key[i];
    printf("%3d %3d    %50.47lf    %-45s\n",p[k],q[k],max[k],file[k]);
  }
}
