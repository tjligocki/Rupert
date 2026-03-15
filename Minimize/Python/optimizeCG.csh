#!/bin/tcsh -f
set ms=$1
set me=$2
set ns=$3
set ne=$4

set m=$ms
while ($m <= $me)
  set n=$ns
  if ($ns < `expr $m + 1`) then
    set n=`expr $m + 1`
  endif
  while ($n <= $ne)
    optimizeCG.py $m $n 1e-3 1e-16
    @ n++
  end
  @ m++
end
