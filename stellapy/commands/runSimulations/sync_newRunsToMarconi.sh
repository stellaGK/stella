#!/bin/bash

# ======================================
# Script to upload new runs to marconi
# ======================================

# Location of the runs on the local computer
here=$NEWRUNS

# Location of the runs on MARCONI
there=$RUNS_MARCONI

# Output message to show the process
echo "----------- START SYNCHRONIZATION NEWRUNS --------------"
echo "Syncronization between the NEWRUNS and RUNS src directories:"
echo " "$here
echo " "$there

# Upload the new runs to marconi (small l because we're copying symlinks that only exist on marconi, it doesn't exist locally which is intended)
rsync -avzh -l --progress --exclude '*_ref' --exclude '*.dat' --exclude '*.pdf' $here $there

echo " "

exit
