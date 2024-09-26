#!/bin/bash
cd "$(dirname "$0")"

# create the header
cp namelist.hdr ../namelists.md

for i in `cat order.txt`; do
   cat $i.nl >> ../namelists.md
   echo -en "\n\n" >> ../namelists.md
done
