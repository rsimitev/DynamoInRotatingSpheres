#!/bin/bash
function usage()
{
   echo -e "Plots the m-spectra of all quantities for the highest numbered run."
   echo -e "Spectra are normalised to the higest value.\n"
   echo -e "Usage:"
   echo -e "\t`basename $0` <runname>"
   echo -e "\t`basename $0` -h|--help"
   echo -e "\t <runname> is the unnumbered  basename"
   echo -e "\t or the numbered  basename of the files containing the data."
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
# Create the plots folder if required
mkdir -p l--plots--l

runName=${runName%.}
# If we already have a number, use that
if [[ -f $runName.par ]]
then
   latest=$runName
else
   # Otherwhise, find the latest state in this folder
   latest=${runName}.`find . -iname $runName'.*.par' | sed -e 's/\.\/'$runName'.//' | sort -n | tail -n 1|sed -e 's/.par//'`
fi

cat - << EOF > l--plots--l/$latest.norm.mspec.gpi
set datafile fortran
plot '$latest.mspec' u 1:2 w l
tempmax=GPVAL_DATA_Y_MAX
plot '$latest.mspec' u 1:3 w l
compmax=GPVAL_DATA_Y_MAX
plot '$latest.mspec' u 1:4 w l
flowmax=GPVAL_DATA_Y_MAX
set xrange [0:20]
set xlabel 'm'
set yrange [1.0e-7:1]
set logscale y
set ylabel 'Normalised power'
set terminal postscript eps enhanced color lw 2 16 size 10 cm, 7.7 cm
set output 'l--plots--l/$latest.norm.mspec.eps'
plot '$latest.mspec' u 1:(\$2/tempmax) w lp lw 2 t 'Temperature', \
     '$latest.mspec' u 1:(\$3/compmax) w lp lw 2 t 'Composition', \
     '$latest.mspec' u 1:(\$4/flowmax) w lp lw 2 t 'Flow'
show output
EOF
gnuplot l--plots--l/$latest.norm.mspec.gpi
