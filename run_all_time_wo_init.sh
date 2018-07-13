##!/usr/bin/bash
source ./parameters.sh

time_slices=($(cat input/tsProf/$shotnumber.time))

for time in "${time_slices[@]}";  do

	echo $time s

        sed -i "s/.*time=.*/time='$time'/" parameters.sh
	./init_renate.sh

	./run_renate.sh
	./plot.sh
done
