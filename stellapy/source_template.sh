#=====================================================================================
#                               Stellapy source file                                 #
#=====================================================================================
# Edit the file below so that each path is correct, which will ensure that the 
# variables required to run stellapy are loaded. To load this file, add the    
# following line to your ~/.alias or .bashrc file which is loaded on start-up. 
#      source /.../stellapy/source.sh                                          
#=====================================================================================


#=====================================================================================
#                                 Stellapy directory                                 #
#=====================================================================================
# Location of the stellapy directory, used to access os.environ.get('STELLAPY') 
# within python in order to import the stellapy modules.
#=====================================================================================
export STELLAPY='MANDATORY/stellapy/'

#=====================================================================================
#                             Load the stellapy commands                             #
#=====================================================================================
# The stellapy commands are defined in "commands.sh", source this file in order
# to have the commands loaded into your command prompt. Personal commands can 
# be defined in a subfolder <User>, as mentioned in the stellapy.ini file.
#=====================================================================================
source $STELLAPY/commands.sh
#source $STELLAPY/OPTIONAL/commands.sh

#=====================================================================================
#                    Make sure python loads the stellapy package                     #
#=====================================================================================
# Uncomment the following export if you'd like to make the stellapy package
# automatically available in the '>> python3' interactive prompt. Make sure
# that this export doesn't clash with previously defined "PYTHONSTARTUP" exports.
#=====================================================================================
#export PYTHONSTARTUP=$STELLAPY/utils/commandprompt/python_startup.py

#=====================================================================================
#                             Other useful short commands                            #
#=====================================================================================
# Commands commonly used when working with stella, comment/uncomment as you see fit.
#=====================================================================================
alias touchall='find . -type f -exec touch {} \;'
alias diskusage='du -h --max-depth=1 | sort -hr'
alias filesizes='ls -Slha --block-size=mB -d input*'
alias filesizesr='ls -Srlha --block-size=mB -d input*'
alias findSaturatedFiles='find . -type f -name "*out.nc.t*"'
alias findSymLink='readlink -f'
alias remove_excessOutputFiles='rm -f $(find . -maxdepth 3 -type f -name "*.out" ! -name "*dt*" ! -name "slurm*"); rm *.fluxes; rm *.omega; rm */*.omega; rm */*.fluxes; rm */*/*.omega; rm */*/*.fluxes'

#=====================================================================================
#                                   Marconi commands                                 #
#=====================================================================================
# Commands commonly used when working with stella, comment/uncomment as you see fit.
#=====================================================================================
if [[ $HOME == *"marconi"* ]]; then
    alias emacs='emacs -nw'
    alias qsub='sbatch'
    alias qdel='scancel'
    alias quota='saldo -b --skl'
    #export USER='OPTIONAL'
    #alias run='squeue -u $USER --format "        %45j %10i %10T %6D %8M   %9L     %15a   %10Q"'
    #alias run_default='squeue -u $USER'
    #alias run_submissiondata=' squeue -u $USER --format "        %10i %10T %6D %8M   %8L    %10r %20V   %20j"'
    #alias run_directory='squeue -u $USER --format "        %7i  %5D  %8M   %110Z"'
    alias run_everyone='squeue  --format "        %10u %20j %10i %10T %6D %8M   %9L     %15a   %10Q"'
    alias stella="source $STELLAPY/supercomputer/load_modulesToCompileStella.sh stella"  
fi

#=====================================================================================
#                                 Supercomputer logins                               #
#=====================================================================================
# On the local computer, load the ssh logins, comment/uncomment as you see fit.
#=====================================================================================
#if [[ $HOME != *"marconi"* ]]; then
    #alias marconi='ssh -X OPTIONAL@login.marconi.cineca.it'
    #alias xula='ssh -X OPTIONAL@xula01.ciemat.es'
    #alias seracre='ssh -X OPTIONAL@seracre.ciemat.es'
    #alias marenostrum='ssh -Y OPTIONAL@mn1.bsc.es'
#fi
