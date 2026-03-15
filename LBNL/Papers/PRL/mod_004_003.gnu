# set xrange [1:15]
set yrange [-0.51:0.01]
set key on

plot   "Split/m008mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 8"
replot "Split/m012mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 12"
replot "Split/m016mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 16"
