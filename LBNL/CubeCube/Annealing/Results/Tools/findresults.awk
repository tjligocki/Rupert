NF == 3 {
  file=substr($1,1,length($1)-1);
  printf("%3d %3d     %s    %s\n",$2,$3,"-------------------",file);
}
NF == 4 {
  file=substr($1,1,length($1)-1);
  printf("%3d %3d     %s    %s\n",$2,$3,$4,file);
}
