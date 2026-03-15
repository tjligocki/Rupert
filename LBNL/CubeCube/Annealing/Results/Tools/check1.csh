#!/bin/csh -f
foreach i ($1/run.*)
  set file=`basename $i` 
  set n=$file:e
  set part=$file:r
  set m=$part:e
  echo $m $n | awk '{printf("%3d %3d  ",$1,$2);}'
  check1 $m $n < $i
end
