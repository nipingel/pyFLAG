#!/bin/bash

## change variables
projID='AGBT16B_400'
session1='12'
session2='13'
obj='NGC6946'
ra=308.71801
dec=60.153915
for i in {0..7}
do
        bm=$i
	if [ $bm -lt 7 ]
	then
        	gbtgridder "../../../../AGBT16B_400_12/calibratedSpectra/NGC6946/scaled/chanBlank/"$projID"_"$session1"_"$obj"_Beam"$bm"_edge_ss_chanBlank_shift.fits" "../../../../AGBT16B_400_13/calibratedSpectra/NGC6946/scaled/chanBlank/"$projID"_"$session2"_"$obj"_Beam"$bm"_edge_ss_chanBlank.fits" --output=$projID"_"$session1"_"$session2"_"$obj"_Beam"$bm -k gaussbessel --size 128 128 --pixelwidth 105 --restfreq 1420.40575 --noweight --noline --nocont	
	fi
	if [ $bm -eq 7 ] 
	then
		gbtgridder "../../../../AGBT16B_400_12/calibratedSpectra/NGC6946/scaled/chanBlank/"$projID"*.fits" "../../../../AGBT16B_400_13/calibratedSpectra/NGC6946/scaled/chanBlank/"$projID"*.fits" --output=$projID"_"$session1"_"$session2"_"$obj"_Comb" -k gaussbessel --size 128 128 --pixelwidth 105 --restfreq 1420.40575 --noweight --noline --nocont
	fi
done
