##!/usr/bin/bash
source ./parameters.sh

#./taiga.exe $shotnumber"_"$time $particles $angle
#            <------ 1 ------->  <-- 2 ---> <- 3 -> <-4->  <-- 5 --> <-- 6 -->
./taiga.exe $shotnumber"_"$time $beammatter $energy $angle $diameter $particles
