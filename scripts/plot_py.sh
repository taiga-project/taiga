##!/usr/bin/bash
source ./parameters.sh
runnumber=`cut -d ',' -f2 runnumber`  

python3 plotter/traj_plotter.py $shotAndTime $runnumber
python3 plotter/ion_profile.py $shotAndTime $runnumber
python plotter/detector_plane.py $shotnumber $time $runnumber $detector
python plotter/detector.py $shotnumber $time $runnumber $detector
