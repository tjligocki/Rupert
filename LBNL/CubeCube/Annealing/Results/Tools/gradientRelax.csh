#!/bin/csh -f
set nn=002
while ($nn < 50)
  foreach i ($1/run.*.$nn)
    set file=`basename $i` 
    set n=$file:e
    set part=$file:r
    set m=$part:e
    echo $file
    set t=`expr $n + 4`
    gradientRelax $m $n < $i | tail -$t >! $2/$file
  end
  @ nn++
  set nn=`echo $nn | awk '{printf("%03d",$1);}'`
end
