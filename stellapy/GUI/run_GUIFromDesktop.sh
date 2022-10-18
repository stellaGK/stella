#!/bin/sh

# Make sure $STELLAPY is loaded 
. ./source.sh

# Run the GUI by executing the python script, print crash logs to stella_GUI.log 
rm $STELLAPY/GUI/stella_GUI.log
/usr/bin/python3 $STELLAPY/GUI/stella_GUI.py >> $STELLAPY/GUI/stella_GUI.log 2>&1
