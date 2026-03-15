{
  if (last != $1) {
    print "";
    last = $1;
  }
  printf("%3d %3d    %50.47lf    %-45s\n",$1,$2,$3,$4);
}
