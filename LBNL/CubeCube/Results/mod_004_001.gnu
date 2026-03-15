set xrange [1:15]
set yrange [-0.51:0.01]
set key on

plot   "Split/m006mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 6"
replot "Split/m010mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 10"
replot "Split/m014mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 14"
replot "Split/m018mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 18"
