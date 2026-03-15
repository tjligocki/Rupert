BEGIN {
  prevM = -1;
  prevN = -1;
  prevV = "";
}
$2 != prevN {
  if (prevN != -1) {
    printf("%3d %3d     %s\n",prevM,prevN,prevV);
  }

  print "";

  prevM = $1;
  prevN = $2;
  prevV = $3;

  next;
}
$2 == prevN {
  if (prevN != -1) {
    prevPower = substr(prevV,index(prevV,"e-")+2);
    power = substr($3,index($3,"e-")+2);

    if (prevPower == power) {
      for (i = 1; i < index(prevV,"e-"); i++) {
        if (substr(prevV,i,1) != substr($3,i,1)) {
          break;
        }
      }

      if (substr(prevV,i,1) < substr($3,i,1)) {
        prevV = substr(prevV,1,i-1) " " substr(prevV,i);
      }
    }

    printf("%3d %3d     %s\n",prevM,prevN,prevV);
  }

  prevM = $1;
  prevN = $2;
  prevV = $3;

  next;
}
END {
  if (prevN != -1) {
    printf("%3d %3d     %s\n",prevM,prevN,prevV);
  }

  print "";
}
