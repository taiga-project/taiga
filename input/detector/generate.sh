#!/bin/bash
maxwidth=10.5
width=1.2
gap=1.0
filename="final/detx"

points=$(bc <<< "2*$maxwidth/($width+$gap)")
echo $points
rm $filename
for i in `seq 0 $points`;
do
    left=$(bc <<< "-$maxwidth+$i*($width+$gap)")
    right=$(bc <<< "-$maxwidth+$i*($width+$gap)+$width")
    echo -e "$left\t$right" >> $filename
done
