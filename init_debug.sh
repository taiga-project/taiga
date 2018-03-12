##!/usr/bin/bash
source ./parameters.sh

./sync.sh  $shotnumber

matlabscript="cd preproc, cdb_reader('$shotnumber',$time, $electric_field_module, '$electric_field_value')"
eval '$matlab -nodesktop -r "$matlabscript"'
