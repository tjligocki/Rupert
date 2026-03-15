{
  comp = NF;
  start = 0;
  max = 0.0;
  for (i = 1; i <= comp; i++) {
    if ($i > 0) {
      cur = $i;
    } else {
      cur = -$i;
    }
    if (cur > max) {
      max = cur;
      if ($i < 0) {
        neg = 1;
      } else {
        neg = 0;
      }
    }
  }
  for (i = 1; i <= comp; i++) {
    if ($i == 0.0) {
      if (substr($i,1,1) == "-") {
        printf(" %s ",substr($i,2));
      } else {
        printf(" %s ",$i);
      }
    } else {
      if (neg == 1) {
        if (substr($i,1,1) == "-") {
          printf(" %s ",substr($i,2));
        } else {
          printf("-%s ",$i);
        }
      } else {
        if (substr($i,1,1) == "-") {
          printf("%s ",$i);
        } else {
          printf(" %s ",$i);
        }
      }
    }
  }
  printf("\n");
}
