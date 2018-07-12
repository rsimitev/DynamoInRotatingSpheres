#!/bin/bash
# Prints the usage and exits
function usage()
{
   echo -e "Rebuilds the global time series from the individual runs."
   echo -e 'Usage:' 
   echo -e "\t$1 <runname>"
   echo -e "\t$1 -h|--help"
   echo '<runname> is the unnumbered basename of the files containing the data.'
   exit
}

# Tests a file for NaN's
function testForNaNs()
{
grep -qi '\bNaN\b' $1 && true || false
}


if [[ $# != 1 || $1 == '-h' || $1 == '--help' ]]
then
   usage $0
fi

runName=$1

extensions="eb ek nu koeb koeu dissB dissu"
for ext in ${extensions}
do
   seriesFile=$runName.$ext
   allFiles=`ls -f $runName.[0-9]*.$ext 2> /dev/null`
   # If we have no  files with this extension, continue
   test -n "${allFiles}" || continue
   # Extract the run numbers available
   numbers=`echo ${allFiles} | sed -e "s/$runName\.//g" -e "s/\.$ext//g" -e's/ /\n/g' | sort -n`
   echo ${numbers}
   test -z "${numbers}" && exit 2
   cat /dev/null > $seriesFile
   for n in $numbers
   do
      seriesFileN="$runName.$n.$ext"
      if [[ -f $seriesFileN ]]
      then
         # Test for NaN's before continuing
         if testForNaNs $seriesFileN 
         then
            echo "$seriesFileN contains NaN's. Stoping here!"
            exit 1
         fi
         cat $seriesFileN >> $seriesFile
      fi
   done
done

