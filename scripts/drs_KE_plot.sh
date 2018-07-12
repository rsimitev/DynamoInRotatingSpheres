#!/bin/bash

function usage()
{
   echo -e "Plots the flow energy components.\n"
   echo -e "Usage:"
   echo -e "\t`basename $0` <runname>"
   echo -e "\t`basename $0` -h|--help"
   echo -e "\t <runname> is the unnumbered  basename of the files containing the data."
   exit 1
}
if [[ ${#} != 1 ]]
then
   usage
elif [[ $1 == '-h' || $1 == '--help' ]]
then
   usage
fi
runName=$1

# Rebuilds the time-series of whatever needs to be rebuilt
rebuildTimeSeries.sh $runName

# Create the plots folder if required
mkdir -p l--plots--l

cat -  << EOF > l--plots--l/plot_KE.gpi
reset
set datafile fortran
set key below
set logscale y
set xlabel 'time'
set ylabel 'Energy'
set output 'l--plots--l/KineticEnergy.eps'
set terminal postscript eps enhanced color lw 2 16 size 10 cm, 7.7 cm
plot \
      '$runName.ek' using 1:2  w l lw 2 lc 0 t 'Total KE', \
      '$runName.ek' using 1:3  w l lw 2 lc 1 t 'zonal pol e-s', \
      '$runName.ek' using 1:4  w l lw 2 lc 1 t 'zonal tor e-s', \
      '$runName.ek' using 1:7  w l lw 2 lc 2 t 'zonal pol e-a', \
      '$runName.ek' using 1:8  w l lw 2 lc 2 t 'zonal tor e-a', \
      '$runName.ek' using 1:5  w l lw 2 lc 3 t 'non-zonal pol e-s', \
      '$runName.ek' using 1:6  w l lw 2 lc 3 t 'non-zonal tor e-s', \
      '$runName.ek' using 1:9  w l lw 2 lc 4 t 'non-zonal pol e-a', \
      '$runName.ek' using 1:10 w l lw 2 lc 4 t 'non-zonal tor e-a'
show output
EOF

gnuplot l--plots--l/plot_KE.gpi 
