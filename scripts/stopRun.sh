#!/bin/bash
function showUsage()
{
   echo "Usage: `basename $0` <jobscript> [-f][-n]"
   exit 1
}

function softKill()
{
  # Create a stop file for the job itself
  touch ${jobscript}.stop
  # Remove the lock file so the program stops at the end of this step.
  rm -f drs.lock
  echo "Job will stop cleanly."
}

if [[ $# > 2 || $# < 1 ]]
then
   showUsage
fi
jobscript=$1
cd `dirname ${jobscript}`

# If we have no options just do a soft kill
if [[ $# == 1 ]]
then
   softKill
# If we have a -f option then search for the job that is running this jobscript and kill it.
elif [[ $# == 2 ]]
then
   allJobs=`qstat | awk '/'$(hostname -s)'/ {print $1}'`
   jobdir=`pwd`
   for job in $allJobs
   do
      # The xml output is easyer to parse
      workdirXML=`qstat -fx $job | xmllint --format -|grep Output_Path`
      # remove xml tags and hostname
      workdir=`echo $workdirXML|sed -e's/<\?Output_Path>//g' -e 's/.*://'`
      # Extract the dirname
      workdir=`dirname $workdir`
      if [[ $jobdir == $workdir ]]
      then
         if [[ $2 == "-f" ]]
         then
            qdel $job
            echo "Killed job $job."
            rm -f ${jobscript}.stop drs.lock
            exit 0
         elif [[ $2 == "-n" ]]
         then
            echo "Would have killed job $job."
            exit 0
         else
            echo "Invalid option $2"
            exit 2
	 fi
      fi
   done
fi
