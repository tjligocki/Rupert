set xrange [1:15]
set yrange [-0.51:0.01]
set key on

plot   "Split/m001mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 1"
replot "Split/m002mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 2"
replot "Split/m003mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 3"
replot "Split/m004mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 4"
replot "Split/m005mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 5"
replot "Split/m006mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 6"
replot "Split/m007mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 7"
replot "Split/m008mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 8"
replot "Split/m009mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 9"
replot "Split/m010mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 10"
replot "Split/m011mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 11"
replot "Split/m012mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 12"
replot "Split/m013mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 13"
replot "Split/m014mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 14"
replot "Split/m015mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 15"
replot "Split/m016mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 16"
replot "Split/m017mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 17"
replot "Split/m018mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 18"
replot "Split/m019mod000" using ($2/$1):($3*$3 - $2/$1) title "m = 19"
