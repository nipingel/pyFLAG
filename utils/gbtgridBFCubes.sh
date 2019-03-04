#!/bin/bash

## change variables
projID='AGBT16B_400'
session='12'
obj='NGC6946'
ra=308.71801
dec=60.153915
for i in {0..7}
do
        bm=$i
	if [ $bm -lt 7 ]
	then
        	gbtgridder $projID"_"$session"_"$obj"_Beam"$bm"_edge_ss.fits" --output=$projID"_"$session"_"$obj"_Beam"$bm -k gaussbessel --size 128 128 --pixelwidth 105 --restfreq 1420.40575 --noweight --noline --nocont	
	fi
	if [ $bm -eq 7 ] 
	then
		gbtgridder $projID"_"$session"_"$obj"_Beam*_edge_ss.fits" --output=$projID"_"$session"_"$obj"_Comb" -k gaussbessel --mapcenter $ra $dec --size 128 128 --pixelwidth 105 --restfreq 1420.40575 --noweight --noline --nocont
	fi
done
