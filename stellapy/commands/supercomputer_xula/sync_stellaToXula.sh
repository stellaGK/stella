#!/bin/bash

# ==================================================
# Script to synchronize the stella version that is on to local computer to the stella on marconi
# ==================================================

# Location of the stella code on the local computer
here=$STELLA
code='stella'

# Location of the stella code on MARCONI
there="e5368@xula01.ciemat.es:/mnt/lustre/home/e5368/stella/"

# Output message to show the process
echo "Syncronization between the " $code " src directories:"
echo " "$here
echo " "$there

# Syncronize the stella versions
rsync -avzh -L --progress --exclude '.git' --exclude stellapy/stella_dirs.py $here $there

exit
