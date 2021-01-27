#!/bin/bash                                                                                                                                             \
                                                                                                                                                         
# script to load or shift the modules and enviroment variables for different codes (euterpe and stella)                                                 \
                                                                                                                                                         
echo "GK_SYSTEM = "$GK_SYSTEM

module purge     # To purge compiler and libraries modules                                                                                               
source ~/.bashrc # To load again those that need to be set for other matters (e.g. python)                                                               

if [ "$GK_SYSTEM" == "marconi" ]; then

    if [ "$1" == "euterpe" ]; then
        module load profile/advanced
        module load intel/pe-xe-2018--binary
        module load intelmpi/2018--binary
        module load szip/2.1--gnu--6.1.0
        module load zlib/1.2.8--gnu--6.1.0
        module load hdf5/1.8.18--intelmpi--2018--binary
        module load petsc/3.8.3--intelmpi--2018--binary
        module load mkl
        export EUTERPE_FCOMP=intel16
        export EUTERPE_PETSC=petsc38
        export EUTERPE_OMPI=impi51
        export EUTERPE_SPARSELIB=PETSC

    elif [ "$1" == "stella" ]; then
        # Modules for the old stella version (before November 2020)
        module load intel/pe-xe-2018--binary
        module load intelmpi/2018--binary
        module load zlib/1.2.8--gnu--6.1.0
        module load szip/2.1--gnu--6.1.0
        module load hdf5/1.10.4--intel--pe-xe-2018--binary
        module load netcdf/4.6.1--intel--pe-xe-2018--binary
        module load netcdff/4.4.4--intel--pe-xe-2018--binary
        module load fftw/3.3.7--intelmpi--2018--binary
        module load lapack
        module load blas
        export SFINCS_DIR=/marconi/home/userexternal/jgarciar/stella/version3
        export NETCDF_LIB_DIR=/cineca/prod/opt/libraries/netcdff/4.4.4/intel--pe-xe-2018--binary/lib
        export NETCDF_INC_DIR=/cineca/prod/opt/libraries/netcdff/4.4.4/intel--pe-xe-2018--binary/include
        export FFTW_LIB_DIR=$FFTW_DIR/lib
        export FFTW_INC_DIR=$FFTW_INC
        export MAKEFLAGS=-IMakefiles

    elif [ "$1" == "stellaGithub" ]; then
        # Modules for the new stella version (after November 2020)
        module purge
        module load env-skl
        module load intel/pe-xe-2018--binary 
        module load intelmpi/2018--binary
        module load szip/2.1--gnu--6.1.0
        module load zlib/1.2.8--gnu--6.1.0
        module load hdf5/1.10.4--intel--pe-xe-2018--binary
        module load netcdf/4.6.1--intel--pe-xe-2018--binary
        module load netcdff/4.4.4--intel--pe-xe-2018--binary
        module load fftw/3.3.7--intelmpi--2018--binary 
        module load mkl/2018--binary 
        module load scalapack
        module load blas  
        export SFINCS_DIR=/marconi/home/userexternal/jgarciar/stella/version3
        export NETCDF_LIB_DIR=/cineca/prod/opt/libraries/netcdff/4.4.4/intel--pe-xe-2018--binary/lib
        export NETCDF_INC_DIR=/cineca/prod/opt/libraries/netcdff/4.4.4/intel--pe-xe-2018--binary/include
        export FFTW_LIB_DIR=$FFTW_DIR/lib
        export FFTW_INC_DIR=$FFTW_INC
        export MAKEFLAGS=-IMakefiles

    elif [ "$1" == "sfincs" ]; then
        module load intel/pe-xe-2017--binary
        module load intelmpi/2017--binary
        module load petsc/3.7.5--intelmpi--2017--binary
        module load netcdf/4.4.1--intel--pe-xe-2017--binary
        module load netcdff/4.4.4--intel--pe-xe-2017--binary
        module load fftw/3.3.4--intelmpi--2017--binary
        module load lapack
        module load blas
        module load szip/2.1--gnu--6.1.0
        module load zlib/1.2.8--gnu--6.1.0
        module load hdf5/1.8.17--intel--pe-xe-2017--binary
        export NETCDF_LIB_DIR=/cineca/prod/opt/libraries/netcdff/4.4.4/intel--pe-xe-2017--binary/lib
        export NETCDF_INC_DIR=/cineca/prod/opt/libraries/netcdff/4.4.4/intel--pe-xe-2017--binary/include
    fi

elif [ "$GK_SYSTEM" == "mn4" ]; then

    if [ "$1" == "euterpe" ]; then
        module load impi petsc/3.7.6-real fftw/3.3.6 fabric/1.4.2 netcdf
        export FC=mpiifort
    elif [ "$1" == "stella" ]; then
        module load fftw/3.3.6
        module load impi/2017.6
        module load netcdf/4.4.1.1
        export NETCDF_INC_DIR=$NETCDF_INC
        export NETCDF_LIB_DIR=$NETCDF_LIB
        export FFTW_INC_DIR=-I$FFTW_INCL
        export FFTW_LIB_DIR=$FFTW_DIR/lib
        export MAKEFLAGS=-IMakefiles
        export PETSC_INC
        export FC=ifort
    fi

elif [ "$GK_SYSTEM" == "" ]; then
    echo "GK_SYSTEM is empty."
fi

echo "Modules and evironment set for "$1"."

