#!/bin/bash
# short script to regrid PAF cubes to single pixel spatial/spectral scale

dataDir="/Users/npingel/Desktop/Research/FLAG/data/AGBT16B_400/"
singlePixelCube="N6946_ALL_GBTMAP_SM.mir"
PafPath=dataDir+"AGBT16B_400_12+13/CUBES/Combined"

NGC6946_BEAM0_CUBE_BASE_XX.FITS
for bm in {0..7}
do 
	for polNum in {1..1}
	do
		## set polarization to process
		if [ $polNum -eq 0 ]
		then
			pol='II'
		fi
		if [ $polNum -eq 1 ]
		then
			pol='XX'
		fi
		if [ $polNum -eq 2 ]
		then
			pol='YY'
		fi
		## if bm is under 7, process the individual beams
		if [ $bm -lt 7 ]
		then
			fits in="NGC6946_BEAM"$bm"_CUBE_BASE_"$pol".FITS" op=xyin out="NGC6946_BEAM"$bm"_CUBE_BASE_"$pol".mir"
			regrid in="NGC6946_BEAM"$bm"_CUBE_BASE_"$pol".mir" out="NGC6946_BEAM"$bm"_CUBE_BASE_"$pol".regrid" axes="1,2,3" tin="../../../SinglePixel/"$singlePixelCube
			fits in="NGC6946_BEAM"$bm"_CUBE_BASE_"$pol".regrid" op=xyout out="NGC6946_BEAM"$bm"_CUBE_BASE_"$pol"_SPREGRID.FITS"


		fi
		## if bm is equal to seven process the combined cube
		if [ $bm -eq 7 ]
		then
			fits in="NGC6946_ALL_CUBE_BASE_"$pol".FITS" op=xyin out="NGC6946_ALL_CUBE_BASE_"$pol".mir"
			regrid in="NGC6946_ALL_CUBE_BASE_"$pol".mir" out="NGC6946_ALL_CUBE_BASE_"$pol".regrid" axes="1,2,3" tin="../../../SinglePixel/"$singlePixelCube
			fits in="NGC6946_ALL_CUBE_BASE_"$pol".regrid" op=xyout out="NGC6946_ALL_CUBE_BASE_"$pol"_SPREGRID.FITS"
		fi
	done
done

rm -rf *.mir
rm -rf *.regrid
