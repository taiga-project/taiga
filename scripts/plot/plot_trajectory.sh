##!/usr/bin/bash
source ./parameters.sh
runnumber=`cut -d ',' -f2 runnumber`
shotAndTime=$shotnumber\_$time

python3 plotter/traj_plotter.py $shotAndTime $runnumber
