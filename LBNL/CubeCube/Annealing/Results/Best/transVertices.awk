BEGIN {
  i = 1;
}
<<<<<<< .mine
NF == n+2 {
  for (j = 2; j <= n+1; j++) {
    t[i "," j-1] = $j;
||||||| .r1717
{
  for (j = 1; j <= n; j++) {
    t[i "," j] = $j;
=======
NR > 1 && NR <= n+1 {
  for (j = 1; j <= n; j++) {
    t[i "," j] = $(j+1);
>>>>>>> .r1778
  }
  i++;
}
END {
  for (i = 1; i <= n; i++) {
    if (i <= m) {
      v[i] = -0.5;
    } else {
      v[i] =  0;
    }
  }

  while (1) {
    max = 0.0;

    for (i = 1; i <= n; i++) {
      w[i] = 0.0;
      for (j = 1; j <= n; j++) {
        w[i] = w[i] + t[i "," j] * v[j];
      }

      if (abs(w[i]) > max) {
        max = abs(w[i]);
      }
    }

    max = 2*max;

    for (i = 1; i <= n; i++) {
<<<<<<< .mine
      printf("%15.12lf ",w[i]/max);
||||||| .r1717
      printf("%30.27lf ",w[i]/max);
=======
      printf("%15.12lf ",w[i]/max/2.0);
>>>>>>> .r1778
    }
<<<<<<< .mine
    printf("--- ");
    for (i = 1; i <= m; i++) {
      printf("%15.12lf ",v[i]);
    }
||||||| .r1717
=======
    printf(" - ");
    for (i = 1; i <= m; i++) {
      printf("%15.12lf ",v[i]);
    }
>>>>>>> .r1778
    printf("\n");

    i = 1;
    while (i <= m && v[i] == 0.5) {
      v[i] = -0.5;
      i++;
    }

    if (i > m) {
      break;
    }

    v[i] = 0.5;
  }
}

function abs(x) {
  if (x < 0) {
    return -x;
  } else {
    return x;
  }
}
