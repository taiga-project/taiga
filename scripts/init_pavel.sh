##!/usr/bin/bash
source ./parameters.sh
matlabscript="cd preproc,try, pavel_profile_to_taiga('$shotnumber',"$time"), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code (0): $?"
