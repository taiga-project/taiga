##!/usr/bin/bash
source ./parameters.sh

./sync.sh  $shotnumber

matlabscript="cd preproc,try, cdb_reader('$shotnumber',$time), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code: $?"
