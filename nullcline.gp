#! /usr/bin/gnuplot
filename="IntpltdFr.dat"
filename1="IntpltdFt.dat"
set xrange [0:0.4]
set yrange [0:0.4]
set size ratio 1
set contour base 
set view map
unset surface
set cntrparam levels discrete 0
splot filename u ($1):($2):($3) w l ,\
filename1 u ($1):($2):($3) w l
pause -1 "press [Enter] key or [OK] button to quit"