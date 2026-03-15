set xrange [1:15]
set yrange [-0.51:0.01]
set key on

plot   "Split/m004mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 4"
replot "Split/m006mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 6"
replot "Split/m008mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 8"
replot "Split/m010mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 10"
replot "Split/m012mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 12"
replot "Split/m014mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 14"
replot "Split/m016mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 16"
replot "Split/m018mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 18"
