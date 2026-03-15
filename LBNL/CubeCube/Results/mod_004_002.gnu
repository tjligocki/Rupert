set xrange [1:15]
set yrange [-0.51:0.01]
set key on

plot   "Split/m007mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 7"
replot "Split/m011mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 11"
replot "Split/m015mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 15"
replot "Split/m019mod004" using ($2/$1):($3*$3 - $2/$1) title "m = 19"
