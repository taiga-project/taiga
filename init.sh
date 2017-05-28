##!/usr/bin/bash
source ./parameters.sh

./sync.sh  $shotnumber

matlabscript="cd preproc,try, cdb_reader('$shotnumber',$time), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code: $?"

cp "input/renate110/flux_compass"$shotnumber"_"$time".txt" renate_110/data
cp "input/renate110/nt_compass"$shotnumber"_"$time".txt" renate_110/data

idl -rt=renate_110/am_renate110taiga.sav  -args $shotnumber $time $beammatter $energy $angle 200 2.3 'compass'

cp 'renate_110/output/'$(ls -1t renate_110/output/| head -1) "input/renate110/pop_compass"$shotnumber"_"$time".txt"

matlabscript="cd preproc,try, renate110_to_taiga('$shotnumber',"$time", "$angle"), catch, exit(1), end, exit(0);"
eval '$matlab -nodesktop -r "$matlabscript"'
echo "matlab exit code: $?"
