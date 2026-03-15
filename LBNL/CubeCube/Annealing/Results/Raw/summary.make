#!/bin/csh -f
# Process the results
/bin/rm -f summary.accum
set count = 0
foreach i (summary.list $*)
  echo $i
  awk -f summary.awk1 count=$count < $i >> summary.accum
  @ count++
end
sort +1n -2 +0n -1 +4r -5 +3n -4 < summary.accum > summary.sort
awk -f summary.awk2 < summary.sort > summary.body
# Make a LaTeX document
cat summary.head summary.body summary.tail > summary.tex
# Start with a LaTeX document and produce a DVI and PostScript version.
set target=summary
latex $target.tex
latex $target.tex
latex $target.tex
dvips -f -D 600 -t landscape -p 2 $target.dvi > $target.ps
