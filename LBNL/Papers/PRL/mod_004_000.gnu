# set xrange [1:15]
set yrange [-0.51:0.01]
set key on

plot   "Split/m005mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 5"
replot "Split/m009mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 9"
replot "Split/m013mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 13"
replot "Split/m017mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 17"
