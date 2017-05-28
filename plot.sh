##!/usr/bin/bash
source ./parameters.sh
runnumber=`cut -d ',' -f2 runnumber`  

matlabscript="cd plotter,try, detPlot('$shotnumber_$time','$runnumber', "$energy", "$angle"), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code: $?"

matlabscript="cd plotter,try, detPlot('$shotnumber_$time','$runnumber', "$energy", "$angle", '_spx'), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code: $?"
