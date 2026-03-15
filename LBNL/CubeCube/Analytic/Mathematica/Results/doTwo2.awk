BEGIN {
  state = 0;
}
/{{t/ && state == 0 {
  state = 1;
  next;
}
/{t/ && state == 1 {
  var = gensub("{","","g",$1);
  var = gensub(",","","g",var);
  printf("%-3s = ",var);
  val = substr($2,1,length($2)-1);
  state = 2;
  next;
}
state == 2 {
  where = index($1,"`");
  if (where == 0) {
    val = val substr($1,1,length($1)-1);
  } else {
    val = val substr($1,1,where-1);
    printf("%s\n",val);
    state = 1;
  }
}
