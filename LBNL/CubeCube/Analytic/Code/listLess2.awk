BEGIN {
  prevN = -1;
  count = 0;

  maxLen = -1;
  maxIndex = -1;
}
$2 != prevN {
  if (count > 0) {
    for (i = 0; i < count; i++) {
      minLen = maxLen;
      minIndex = maxIndex;
      for (j = 0; j < count; j++) {
        if (length(value[j]) <= minLen && value[j] != "") {
          minLen = length(value[j]);
          minIndex = j;
        }
      }

      print lines[minIndex];
      value[minIndex] = "";
    }
  }

  count = 0;
  lines[count] = $0;
  value[count] = $3;

  maxLen = length(value[count]);
  maxIndex = count;

  count++;

  prevN = $2;

  next;
}
$2 == prevN {
  lines[count] = $0;
  value[count] = $3;
  if (length(value[count]) > maxLen) {
    maxLen = length(value[count]);
    maxIndex = count;
  }
  count++;

  next;
}
END {
  if (count > 0) {
    for (i = 0; i < count; i++) {
      minLen = maxLen;
      minIndex = maxIndex;
      for (j = 0; j < count; j++) {
        if (length(value[j]) <= minLen && value[j] != "") {
          minLen = length(value[j]);
          minIndex = j;
        }
      }

      print lines[minIndex];
      value[minIndex] = "";
    }
  }
}
