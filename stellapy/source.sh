#===============================================================================
#                            Stellapy source file                              #
#===============================================================================
# Edit the file below so that each file path is correct, which will ensure     #
# that the variables required to run stellapy are loaded. To load this file,   #
# add the following line to your ~/.alias file which is loaded on start-up.    #
#      source /home/user/STELLA/stella/stellapy/source.sh                      #
#===============================================================================

#===============================================================================
#                             Stella directories                               #
#===============================================================================
export STELLA='/home/user/STELLA/stella/'
export STELLAPY='/home/user/STELLA/stella/stellapy/'
export MARCONI_STELLA='user@login.marconi.cineca.it:/marconi/home/userexternal/user/stella/' 

#===============================================================================
#                             Runs directories                                 #
#===============================================================================
export RUNS='/home/user/STELLA/RUNS/'
export NEWRUNS='/home/user/STELLA/NEWRUNS/'
export MARCONI_RUNS='user@login.marconi.cineca.it:/marconi_scratch/userexternal/user/SCRATCH/'

#===============================================================================
#                 Make sure python loads the stellapy package                  #
#===============================================================================
export PYTHONSTARTUP=$STELLAPY/utils/commandprompt/python_startup.py

#===============================================================================
#                          Load the stellapy commands                          #
#===============================================================================
source $STELLAPY/commands.sh 

#===============================================================================
#                              Supercomputer logins                            #
#===============================================================================
if [[ $HOME != *"marconi"* ]]; then
    alias marconi='ssh -X user@login.marconi.cineca.it' 
fi

#===============================================================================
#                                Marconi commands                              #
#===============================================================================
if [[ $HOME == *"marconi"* ]]; then
    alias emacs='emacs -nw'
    alias qsub='sbatch'
    alias qdel='scancel'
    alias quota='saldo -b --skl'
fi



