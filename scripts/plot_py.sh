##!/usr/bin/bash
source ./parameters.sh
runnumber=`cut -d ',' -f2 runnumber`  

python plotter/detPlot.py $shotnumber $time $runnumber $detector
