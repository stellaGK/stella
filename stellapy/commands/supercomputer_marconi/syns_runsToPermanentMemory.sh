#!/bin/bash                                                                                                                                              

# ==================================================                                                                                                     
# Syncronize runs from $SCRATCH to $RUNS                                                                                   
# ==================================================                                                                                                     

# Move the reduced data from the $SCRATCH folder to the $FINISHED folder
SOURCE=$SCRATCH
DESTINATION=$FINISHED
echo "Syncronization between the directories:"
echo " "$SOURCE
echo " "$DESTINATION

# Syncronize the stella versions                                                                                                                         
rsync -avzh -L --progress --exclude 'not_important*' --exclude 'Finished' --exclude 'stella' --exclude '*.nc' --exclude 'restart*' --exclude 'restart' --exclude '*.out.nc'  --exclude '*.error' --exclude '*.vmec_geo' --exclude '*.jacob' --exclude '*.dat' --exclude '*.scratch' --exclude 'mat' $SOURCE $DESTINATION

# Save the simulations in the temporary FINISHED folder to the permanent memory of the RUNS folder
SOURCE=$FINISHED
DESTINATION=$RUNS
echo "Syncronization between the directories:"
echo " "$SOURCE
echo " "$DESTINATION

# Syncronize the stella versions                                                                                                                         
rsync -avzh -L --progress --exclude 'not_important*' --exclude 'Finished' --exclude 'stella' --exclude 'stella*' --exclude '*.nc' --exclude 'restart*' --exclude 'restart' --exclude '*.out.nc'  --exclude '*.error' --exclude '*.vmec_geo' --exclude '*.jacob' --exclude '*.dat' --exclude '*.scratch' --exclude 'mat' $SOURCE $DESTINATION

exit

