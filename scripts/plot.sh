##!/usr/bin/bash
source ./parameters.sh
runnumber=`cut -d ',' -f2 runnumber`  

matlabscript="cd plotter, try, detPlot('"$shotnumber"_"$time"','$runnumber', "$energy", "$angle", '_all', 0.18), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code: $?"
