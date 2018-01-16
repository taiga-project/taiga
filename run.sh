##!/usr/bin/bash
source ./parameters.sh

# set runnumber
oldnum=`cut -d ',' -f2 runnumber`  
newnum=`expr $oldnum + 1`
sed -i "s/$oldnum\$/$newnum/g" runnumber 

#./taiga.exe $shotnumber"_"$time $particles $angle
#            <------ 1 -------> <- 2 ->  <-- 3 ---> <- 4 -> <-5->  <-- 6 --> <--- 7 ---> <-- 8 -->
./taiga.exe $shotnumber"_"$time $newnum $beammatter $energy $angle $diameter $detector_R $particles
