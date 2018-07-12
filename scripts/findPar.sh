#!/bin/bash
function usage()
{
	exit 1
}
DATADIR=${1:-.}
par=$2
if [[ $# == 2 ]]
then
   parVal=$2
   if [[ $par == "BC" ]]
   then
      find -H $DATADIR -type d -links 2 | grep -E "_t$parVal|_b$parVal" | sed -e 's/\/l--plots--l//g'| sort -u 
   else
      find -H $DATADIR -type d -links 2 | grep "$par=$parVal" | sed -e 's/\/l--plots--l//g'| sort -u 
   fi
elif [[ $# == 1 ]]
then
   echo "Possible values for $par are:"
   if [[ $par == "BC" ]]
   then
      echo -e "FT\nFdT\nFC\nFdC\nFU\nFdU"
   else
      find $DATADIR -type d -links 2 | grep $par | sed -e 's@.*'$par'=\([0-9\.e]*\).*@\1@' | sort -u 
   fi
   echo "Use $0 <parname> <parval>"
   echo "to see a list of all computations with the given parameter and value."
else
   echo -en "Use\n $0 <parname> <parval>\n"
   echo " to see a list of all computations with the given parameter and value."
   echo -en "Use\n $0 <parname>\n"
   echo " to see a list of all values for the given parameter."
   echo "Possible parameters are:"
   echo -e "Ta\neta\nPt\nPc\nRt\nRc\nBC"
fi

