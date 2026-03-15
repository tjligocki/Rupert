BEGIN {
  count = 0;
  max = 0;
  maxcount = -1;
}
NF == 4 {
  m = $2;
  n = $3;

  result[count] = $4;
  file[count] = substr($1,1,length($1)-1);

  if ($4 > max) {
    max = $4;
    maxcount = count;
  }

  count++;
}
END {
  for (i = 0; i < count; i++) {
    printf("%3d %3d     %s    %13.6e  %s\n",m,n,result[i],(max - result[i])/max,file[i]);
  }
}
