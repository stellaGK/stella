####################################################################
#                                                                  #
#       Makefile for the stella gyrokinetic turbulence code        # 
#                                                                  #
####################################################################
#
# Makefile written by Bill Dorland and Ryusuke Numata
# Cleaned-up and modified by Hanne Thienpondt on 07/2024
#
# In this Makefile the flags for the Fortran and C compilers are set.
#
# In the ~/.source.sh or ~/.bashrc file of your computer define, e.g., 
#     export STELLA_SYSTEM = 'marconi'
# which matches stella/Makefiles/Makefile.$STELLA_SYSTEM
#
# * Available Compilers (tested on limited hosts)
#   (must be Fortran 95 Standard compliant)
#
# Intel ifort
# GNU's gfortran and g95
# IBM XL Fortran xlf90
# PathScale Compiler Suite pathf90
# The Portland Group pgf90
# NAGWare f95 (v5.1)
# Lahey/Fujitsu Fortran lf95
# 
# * Frequently Tested Hosts, Systems
#
# Standard Linux
# Standard Mac OS X with MacPorts
# Franklin at NERSC and Jaguar at NCCS (Cray XT4 with PGI)
# Bassi at NERSC (IBM Power 5 with IBM XL Fortran)
# Ranger (Sun Constellation Linux Cluster with Intel)
####################################################################

# Set MAKEFILE_PATH to the location of this Makefile.
# Using 'abspath' will not resolve symbolic links
MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST))) 
PROJECT_DIR := $(realpath $(dir $(MAKEFILE_PATH))/..)
COMPILATION_DIR := $(PROJECT_DIR)/COMPILATION

####################################################################
#                          COMPILER FLAGS                          #
####################################################################

include $(COMPILATION_DIR)/Makefile.compilerflags

####################################################################
#                              STELLA                              #
####################################################################

include $(COMPILATION_DIR)/Makefile.stella

