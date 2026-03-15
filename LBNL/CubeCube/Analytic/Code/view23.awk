NF == 5 {
  a11 = a21;
  a12 = a22;
  a13 = a23;

  a21 = a31;
  a22 = a32;
  a23 = a33;

  a31 = $1;
  a32 = $2;
  a33 = $3;
}
END {
  sum1 = 0;
  if (a11 > 0) {
    sum1 += a11;
  } else {
    sum1 += -a11;
  }
  if (a12 > 0) {
    sum1 += a12;
  } else {
    sum1 += -a12;
  }
  
  sum2 = 0;
  if (a21 > 0) {
    sum2 += a21;
  } else {
    sum2 += -a21;
  }
  if (a22 > 0) {
    sum2 += a22;
  } else {
    sum2 += -a22;
  }
  
  sum3 = 0;
  if (a31 > 0) {
    sum3 += a31;
  } else {
    sum3 += -a31;
  }
  if (a32 > 0) {
    sum3 += a32;
  } else {
    sum3 += -a32;
  }

  max = sum1;
  if (sum2 > max) {
    max = sum2;
  }
  if (sum3 > max) {
    max = sum3;
  }
  
  print a11,a12,a13,a21,a22,a23,a31,a32,a33,max;
}
