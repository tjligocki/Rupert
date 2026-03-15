{
  if (NR == 1) {
    print $0;
  } else {
    if ($2 != last) {
      print "";
    }
    print $0;
  }

  last = $2;
}
