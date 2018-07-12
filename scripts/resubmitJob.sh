#!/bin/bash

function usage()
{
   echo "Resubmits jobs to a queueing system or a simple shell."
   echo "Usage: `basename $0` <jobscript> {pbs|shell}"
   echo
   echo "At the moment only submiting to pbs ot to a shell is supported."
   exit 1
}

if [[ $# != 2 ]]
then
	usage
fi
jobscript=$1
jobtype=$2
if [[ "x$jobtype" != "xpbs" && "x$jobtype" != "xshell" ]]
then
	usage
fi
cd `dirname ${jobscript}`
if [[ -e drs.lock ]]
then
   echo "Found a lock file."
   echo -n "Should I remove it to ensure the job will run? (Y/n) "
   read -n 1 a
   test "x$a" = "x" -o "x$a" = "xY" -o "x$a" = "xy" && rm drs.lock || exit 2
fi
if [[ -e ${jobscript}.stop ]]
then
   echo "Found a stop file."
   echo -n "Should I remove it to ensure the job will run? (Y/n) "
   read -n 1 a
   test "x$a" = "x" -o "x$a" = "xY" -o "x$a" = "xy" && rm ${jobscript}.stop || exit 3
fi
echo 
num=`cat ${jobscript}.num||echo 1`
echo "Error path is ${jobscript}.err.${num}"
echo "Output path is ${jobscript}.out.${num}"

if [[ $jobtype == "pbs" ]]
then
   qsub -e ${jobscript}.err.${num} -o ${jobscript}.out.${num} ${jobscript}
else
   ./${jobscript} 1> ${jobscript}.out.${num} 2> ${jobscript}.err.${num} & 
fi

