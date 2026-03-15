set xrange [1:15]
set yrange [-0.51:0.01]
set key on

plot   "Split/m006mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 6"
replot "Split/m009mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 9"
replot "Split/m012mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 12"
replot "Split/m015mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 15"
replot "Split/m018mod003" using ($2/$1):($3*$3 - $2/$1) title "m = 18"
