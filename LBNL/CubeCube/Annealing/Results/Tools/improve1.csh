#!/bin/csh -f
foreach i ($1/run.*)
  set file=`basename $i` 
  set n=$file:e
  set part=$file:r
  set m=$part:e
  echo $file
  improve1 $m $n < $i >! $2/$file
end
