##!/usr/bin/bash
source ./parameters.sh

matlabscript="cd preproc,try, cdb_reader('$shotnumber',$time, $electric_field_module, '$electric_field_value'), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code: $?"
