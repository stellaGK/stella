#!/bin/bash                                                                                                                                             \
                                                                                                                                                         
                                                                                                                                                        \
HERE=$PWD

for (( i = 1 ; i <= 9; i++ ))
do
    c=$(printf '%02d\n' $i)
    THERE=$HERE'/run'$c

    if [ ! -d "$THERE" ]; then
        mkdir $THERE
        rsync -avzh -L --progress --exclude 'restart*' --exclude 'run*' --exclude 'stella' --exclude 'wout*' $HERE/  $THERE
        rm .*
        rm *.error
        rm *.err
        rm *.vmec_geo
        rm *.final_fields
        rm *.out.h5
        rm *.scratch
        rm *.geometry
        rm *.h5
        rm *.stella.*
        rm *.species.*
        rm slurm*
        rm "species.input"
        rm "vmec.geo"
        # The following files will be appended
        #rm *.out.nc
        #rm *.omega
        #rm *.fluxes
        #rm *.out
        echo "Run directory created: "$THERE
        exit
    fi
done
