function monitor_header()
{
cat - << EOF 
reset
max(x,y)=(x>y?x:y)
set datafile fortran
EOF
}

function monitor_z_am()
{
   local calc=$1
cat - << EOF 

# Plot axial angular momentum
# collect statistics
stats '$calc.am' u 2:3 prefix "am" nooutput
r=max(max(-am_min_x,am_max_x),max(-am_min_y,am_max_y))
am_stddev=max(am_stddev_x,am_stddev_y)

# Plot axial angular momentum
set xlabel 'time'
set ylabel 'L_z'
plot '$calc.am' u 1:4 w lp notitle
EOF
}

function monitor_h_am()
{
   local calc=$1
cat - << EOF 

# Plot horizontal angular momentum
# collect statistics
stats '$calc.am' u 2:3 prefix "am" nooutput
r=max(max(-am_min_x,am_max_x),max(-am_min_y,am_max_y))
am_stddev=max(am_stddev_x,am_stddev_y)
set size ratio -1
set zeroaxis
unset key
set xlabel 'L_x'
set ylabel 'L_y'
set palette defined (0 "white", 0.5 "light-green", 1 "dark-green", 1.5 "dark-yellow", 2 "orange", 2.5 "red", 3 "dark-red", 3.5 "black")
set xrange[-r:r]
set yrange[-r:r]
set y2range [-1:1]
plot '$calc.am' u 2:3:1 w p pt 7 ps 2 lw 2 lc palette, \
     '' u (sqrt(\$2**2+\$3**2)):(am_stddev/am_records) axes x1y2 w l lw 3 lc 0 smooth kdensity
unset zeroaxis
EOF
}

function monitor_KE()
{
   local calc=$1
cat - << EOF 

# Plot kinetic energies
set xrange[*:*]
set yrange[*:*]
set size 1,1
set xlabel 'time'
set ylabel 'Energy'
set key on
set key below
set logscale y
plot  '$calc.ek' using 1:2  w l lw 2 lc 0 t 'Total kinetic e.', \
      '$calc.ek' using 1:3  w lp lw 2 lc 1 pt 5 t 'zonal pol e-s', \
      '$calc.ek' using 1:4  w lp lw 2 lc 1 pt 6 t 'zonal tor e-s', \
      '$calc.ek' using 1:7  w lp lw 2 lc 2 pt 5 t 'zonal pol e-a', \
      '$calc.ek' using 1:8  w lp lw 2 lc 2 pt 6 t 'zonal tor e-a', \
      '$calc.ek' using 1:5  w lp lw 2 lc 3 pt 5 t 'non-zonal pol e-s', \
      '$calc.ek' using 1:6  w lp lw 2 lc 3 pt 6 t 'non-zonal tor e-s', \
      '$calc.ek' using 1:9  w lp lw 2 lc 4 pt 5 t 'non-zonal pol e-a', \
      '$calc.ek' using 1:10 w lp lw 2 lc 4 pt 6 t 'non-zonal tor e-a'

EOF
}

function monitor_Nusselt()
{
   local calc=$1
cat - << EOF 

# Plot the Nusselt number at the boundaries.
set xlabel 'time'
set ylabel 'Nu'
plot '$calc.nu' u 1:2 w l t 'Nu @ ICB', \
     '$calc.nu' u 1:3 w l t 'Nu @ CMB
EOF
}

function monitor_new_screen()
{
   local panel=$1

cat - << EOF 
set terminal x11 $panel
EOF
}

function monitor_repeat()
{
   local interval=$1
cat - << EOF 

pause $interval
reread
EOF
}

