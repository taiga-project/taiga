##!/usr/bin/bash
source ./parameters.sh
runnumber=`cut -d ',' -f2 runnumber`

python plotter/detector.py $shotnumber $time $runnumber $detector
