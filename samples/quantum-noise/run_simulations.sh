#!/bin/bash
#------------------------------------------------------------------
# Runs the four quantum noise samples and prints the mean
# magnetisation length of each (running mean at the end of the run).
# Usage: ./run_simulations.sh [path-to-vampire-serial]
#------------------------------------------------------------------
VAMPIRE=${1:-../../vampire-serial}
if [ ! -x "$VAMPIRE" ]; then
   echo "Error: vampire binary not found at $VAMPIRE"
   exit 1
fi
for case in input-classical input input-quantum input-open-system; do
   dir=run-${case}
   mkdir -p $dir
   cp Ni.mat $dir/
   cp $case $dir/input
   ( cd $dir && $VAMPIRE > /dev/null )
   # column 2 is the running mean of |m| (output:mean-magnetisation-length)
   mean=$(grep -v '^#' $dir/output | tail -1 | awk '{printf "%.4f", $2}')
   echo "$case: mean |m| = $mean"
done
