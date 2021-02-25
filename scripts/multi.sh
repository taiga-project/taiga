##!/usr/bin/bash
source ./parameters.sh
runnumber=`cut -d ',' -f2 runnumber`  

matlabscript="cd plotter, try, multirun_angled_det('"$shotnumber"_"$time"'), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code: $?"

matlabscript="cd plotter, try, detPlot('"$shotnumber"_"$time"','all', "$energy", "$angle", '_all', 0.22), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code: $?"
