##!/usr/bin/bash

if [ -z $1 ]
then
    parameter_file="./parameters.sh"
else
    parameter_file=$1
fi

source $parameter_file

./sync.sh  $shotnumber

case "$matlab" in
"")
    matlab="matlab"
;;
esac

matlabscript="cd preproc,try, cdb_reader('$shotnumber',$time, $electric_field_module, '$electric_field_value', '$magnetic_field_value'), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
if (( $? == 0 ))
then
    echo "Initialisation was successful"
else
    echo "Error in initialisation. Please run init_debug.sh"
fi

echo $matlabscript>"input/fieldGrid/"$shotnumber"_"$time"/matlab.txt"
