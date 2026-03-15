BEGIN {
  cur = 0;
  num = 0;
}
{
  for (i = 1; i <= NF; i++) {
    x[cur,i-1] = $i;
  }

  cur++;
  num = NF;
}
END {
  numcorn = lshift(1,indim);
  for (i = 0; i < numcorn; i++) {
    for (j = 0; j < indim; j++) {
      if (and(rshift(i,j),1) == 0) {
        corn[i,j] = -1.0;
      } else {
        corn[i,j] =  1.0;
      }
    }

    for (j = indim; j < cur; j++) {
      corn[i,j] = 0.0;
    }
  }

  max = 0.0;

  for (i = 0; i < numcorn; i++) {
    for (j = 0; j < cur; j++) {
      trans[i,j] = 0.0;
      for (k = 0; k < cur; k++) {
        trans[i,j] += x[j,k]*corn[i,k];
      }

      if (trans[i,j] > max) max = trans[i,j];
    }
  }

  printf("%19.16f\n\n",1.0/max);

  for (i = 0; i < numcorn; i++) {
    for (j = 0; j < cur; j++) {
      printf("%19.16f ",trans[i,j]/max);
    }
    printf("\n");
  }
}
