#!/bin/csh -f
set mold=""
foreach i ($1/run.*)
  set file=`basename $i` 
  set n=$file:e
  set part=$file:r
  set m=$part:e
  if ("$mold" != "$m") then
    echo ""
    echo "$m" | awk '{printf("%3d\n",$1);}'
    set mold="$m"
  endif
  check2 $m $n < $i
end
