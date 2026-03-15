# set xrange [1:15]
set yrange [-0.51:0.01]
set key on

plot   "Split/m004mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 4"
replot "Split/m007mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 7"
replot "Split/m010mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 10"
replot "Split/m013mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 13"
replot "Split/m016mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 16"
replot "Split/m019mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 19"
