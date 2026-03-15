BEGIN {
}
{
  if ($3 == "0.e0") {
    printf("%3d %3d     %19.17e\n",$1,$2,0.0);
  } else if (substr($3,1,1) == "0" && substr($3,3,1) != "0") {
    sci = index($3,"e-");
    power = substr($3,sci+2) + 1;
    new = substr($3,3,1) "." substr($3,4,sci-4) "e-" sprintf("%02d",power);
    printf("%3d %3d     %s\n",$1,$2,new);
  } else {
    print $0;
  }
}
