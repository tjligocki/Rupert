{
  m = $1;
  n = $2;

  num = n - m;

  for (div = 1; div <= num; div++) {
    if ((div - num + m) % num == 0) {
      break;
    }
  }

  m1 = (div - num + m) / num + 1;
  m2 = m1 + 1;

  if (div < num) {
    print m2-1;
  } else {
    print m1-1;
  }
}
