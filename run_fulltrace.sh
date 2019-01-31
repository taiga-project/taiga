##!/usr/bin/bash
source ./parameters.sh

# set runnumber
oldnum=`cut -d ',' -f2 runnumber`  
newnum=`expr $oldnum + 1`
sed -i "s/$oldnum\$/$newnum/g" runnumber 

mkdir "results/"$shotnumber"_"$time"/"$newnum
cp parameters.sh "results/"$shotnumber"_"$time"/"$newnum

./taiga.exe --fulltrace
