#!/bin/bash

if [[ $# != 1 ]]
then
   echo "Removes derived quantities leaving only the state files."
   echo "Usage:"
   echo " `basename $0` <base name>"
fi
base=$1
extensions="adv am cfl dissu ek koeu lspec mspec nspec nu t uaz u_mid ur uzon"
for ext in $extensions
do
   rm -vf $base.$ext
done
