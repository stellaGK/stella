#!/bin/bash  

#########################################################
# Reduce the size of the "*.out.nc" and "wout*" file in each subfolder of $RUNS
#########################################################

# Only perform this action on marconi
if [ $HOME = "/marconi/home/userexternal/hthienpo" ]; then

# Load the python environment
source ~/my_pyth/bin/activate

# Location of the bashscripts
BASHSCRIPTS="/marconi/home/userexternal/hthienpo/stella/stellapy/commands"

# Don't flood the command prompt
echo " "
python3 $BASHSCRIPTS/utils/toggle_verbose.py --off

# Run the reduce_size_netcdf and reduce_size_wout commands in each subfolder
for D in `find $SCRATCH -mindepth 1 -maxdepth 10 -type d`
do
    if ! [[ $D == *"olddd"* || $D == *"finished"* || $D == *"restart"* || $D == *"failed"* || $D == *"extra"* || 
            $D == *"not_important"* || $D == *"figures"* || $D == *"stellacode"* ]]; then
    echo "------------------------------------------------------"
    echo "Reduce the size of the netcdf and wout file inside $D:"
    echo "------------------------------------------------------"
    python3 $BASHSCRIPTS/reduceSizeFiles/reduce_sizeNetcdf.py --folder $D #--dontwriteagain
    python3 $BASHSCRIPTS/reduceSizeFiles/reduce_sizeWout.py --folder $D #--dontwriteagain
    echo
    fi
done

# Turn the verbose back on
echo " "
python3 $BASHSCRIPTS/utils/toggle_verbose.py --on
fi


