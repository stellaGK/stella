#!/bin/bash

# ===============================================================================================
# Script to synchronize the stella version that is on to local computer to the stella on marconi
# ================================================================================================

# Location of the stella code on the local computer
here=$STELLA
code='stella'

# Location of the stella code on MARCONI
there=$STELLA_MARCONI

# Output message to show the process
echo "-------- START SYNCHRONIZATION STELLA CODE -------------"
echo "Syncronization between the " $code " src directories:"
echo " "$here
echo " "$there

# Syncronize the stella versions
rsync -avzh -L --progress --exclude '.git' --exclude '*/config/research/*' --exclude '*/config/experiments/*' --exclude '*/config/config.ini' --exclude '__pycache__' --exclude 'linearmap*' --exclude '.project' --exclude '*.aux' --exclude '*.bbl' --exclude '*.blg' --exclude '*.log' --exclude '*.gz' --exclude '*.toc' $here $there

exit
