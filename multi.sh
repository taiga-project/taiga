##!/usr/bin/bash
source ./parameters.sh
runnumber=`cut -d ',' -f2 runnumber`  

matlabscript="cd plotter, multirun_angled_det('"$shotnumber"_"$time"')"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code: $?"
