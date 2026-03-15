BEGIN {
  on = 0;
}
/{{t99[^0-9]/ {
  print $0;
  on = 1;
  next;
}
/{{t100[^0-9]/ {
  on = 0;
  exit;
}
on == 1 {
  print $0;
  next;
}
