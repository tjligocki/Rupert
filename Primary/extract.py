#!/bin/tcsh

foreach m (05 06 07 08 09 10 11 12 13 14 15)
  set r = 0
  while ($r < $m)
    set rr = `echo $r | awk '{printf("%02d",$1);}'`
    echo $m $r $rr
    ./select.py $m $r < ResultsCombinedWithSquares.txt > m${m}_${rr}.txt
    @ r++
  end
end
