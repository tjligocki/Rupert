BEGIN {
  prev = 1;
}
$1 != prev {
  print "";
  print $0;

  prev = $1;

  next;
}
$1 == prev {
  print $0;
  next;
}
