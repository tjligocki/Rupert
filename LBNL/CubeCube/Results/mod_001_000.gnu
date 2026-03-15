set xrange [1:15]
set yrange [-0.51:0.01]
set key on

plot   "Split/m002mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 2"
replot "Split/m003mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 3"
replot "Split/m004mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 4"
replot "Split/m005mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 5"
replot "Split/m006mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 6"
replot "Split/m007mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 7"
replot "Split/m008mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 8"
replot "Split/m009mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 9"
replot "Split/m010mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 10"
replot "Split/m011mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 11"
replot "Split/m012mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 12"
replot "Split/m013mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 13"
replot "Split/m014mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 14"
replot "Split/m015mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 15"
replot "Split/m016mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 16"
replot "Split/m017mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 17"
replot "Split/m018mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 18"
replot "Split/m019mod001" using ($2/$1):($3*$3 - $2/$1) title "m = 19"
