# set xrange [1:15]
set yrange [-0.51:0.01]
set key on

plot   "Split/m003mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 3"
replot "Split/m005mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 5"
replot "Split/m007mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 7"
replot "Split/m009mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 9"
replot "Split/m011mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 11"
replot "Split/m013mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 13"
replot "Split/m015mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 15"
replot "Split/m017mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 17"
replot "Split/m019mod002" using ($2/$1):($3*$3 - $2/$1) title "m = 19"
