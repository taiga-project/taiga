#!/bin/bash
maxwidth=20.0
width=1.0
gap=0.0
filename="detx"

points=$(bc <<< "2*$maxwidth/($width+$gap)-1")
rm $filename
for i in `seq 0 $points`;
do
    left=$(bc <<< "-$maxwidth+$i*($width+$gap)")
    right=$(bc <<< "-$maxwidth+$i*($width+$gap)+$width")
    echo -e "$left\t$right" >> $filename
done
