NF <= 1 {
  print $0;
  next;
}
{
  comp = NF-2;
  printf("%3d",$1);
  for (cur = 2; cur <= comp+1; cur++) {
    if (cur == i+1) {
      theOne = $(j+1);
    } else 
    if (cur == j+1) {
      theOne = $(i+1);
    } else
    {
      theOne = $cur;
    }

    if (substr(theOne,1,1) == "-") {
      printf(" %s",theOne);
    } else {
      printf("  %s",theOne);
    }
  }
  printf("  %s\n",$NF);
  next;
}
