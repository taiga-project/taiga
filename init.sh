##!/usr/bin/bash
source ./parameters.sh

./sync.sh  $shotnumber

matlabscript="cd preproc,try, cdb_reader('$shotnumber',$time), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
if (( $? == 0 ))
then
	echo "Initialisation was successfull"
else
	echo "Error in initialisation. Please run init_debug.sh"
fi
