BEGIN {
  on = 0;
}
/{{t49[^0-9]/ {
  print $0;
  on = 1;
  next;
}
/{{t50[^0-9]/ {
  on = 0;
  exit;
}
on == 1 {
  print $0;
  next;
}
