NF == 3 && $1 ~ /[0-9][0-9]*/ && $2 ~ /[0-9][0-9]*/ {
  mm = $1;
  nn = $2;

  for (i = mm; i >= 1; i--) {
    if ((nn % i == 0) && (mm % i == 0)) {
      m = mm / i;
      n = nn / i;
      break;
    }
  }

  if (m == 1) {
    error = $3 - sqrt(n);
    if (error < 0) error = -error;

    printf("%2d %2d   %13.6e\n",mm,nn,error);
  } else
  if (m == 2) {
    if (n % 2 == 0) {
      error = $3 - sqrt(n/2.0);
    } else {
      error = $3 - sqrt((n-1)/2.0 + 1.0/8.0);
    }

    if (error < 0) error = -error;

    printf("%2d %2d   %13.6e\n",mm,nn,error);
  } else {
    print $0;
  }

  next;
}
{
  print $0;
}
