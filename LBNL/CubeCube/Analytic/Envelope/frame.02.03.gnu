set xrange [-0.5:1.5]
set yrange [-0.5:1.5]
set trange [-0.5:1.5]
set size ratio -1
set parametric
plot t,-3.000000*t-0.000002 lw 2 lc 1
replot t,-3.000000*t+4.000002 lw 2 lc 1
replot t,3.000000*t+1.000002 lw 2 lc 1
replot t,3.000000*t-3.000002 lw 2 lc 1
replot t,-0.000000*t-0.000000 lw 2 lc 1
replot t,-0.000000*t+1.000000 lw 2 lc 1
set style arrow 1 nohead lw 2 front
set arrow from 0,0 to 1,0 as 1
set arrow from 1,0 to 1,1 as 1
set arrow from 1,1 to 0,1 as 1
set arrow from 0,1 to 0,0 as 1
replot NaN,NaN
