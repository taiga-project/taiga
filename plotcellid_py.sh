##!/usr/bin/bash
source ./parameters.sh
runnumber=`cut -d ',' -f2 runnumber`  

python plotter/detPlot_cellid.py $shotnumber $time $runnumber $detector
