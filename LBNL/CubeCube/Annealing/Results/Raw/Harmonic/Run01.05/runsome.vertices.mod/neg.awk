NF < 2 {
  print $0;
  next;
}
NF >= 2 {
  printf("%3d ",$1);
  for (i = 2; i < NF; i++) {
    if (substr($i,1,1) == "-") {
      printf(" %s ",substr($i,2));
    } else {
      printf("-%s ",$i);
    }
  }
  printf(" %s \n",$NF);
}
