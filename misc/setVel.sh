#!/bin/bash 

## short shell script to convert cubes from frequency to velocity using miriad
flagDataRoot='AGBT16B_400_12+13_NGC6946_Beam'
combName='AGBT16B_400_12_13_NGC6946_Comb'

for bm in {0..7}
do
	## if not combined cube
	if [ "$bm" -lt 7 ]
	then
		## read in FLAG cubes
        	fits in=$flagDataRoot$bm'_cube.fits' op=xyin out=$flagDataRoot$bm'_cube.mir'

        	## set correct velocity frame in FLAG data cube
        	velsw in=$flagDataRoot$bm'_cube.mir' axis='FELO'
	
		## delete old FITS files
		rm $flagDataRoot$bm'_cube.fits'

		## output new FITS files
		fits in=$flagDataRoot$bm'_cube.mir' out=$flagDataRoot$bm'_cube.fits' op=xyout
	fi
	if [ "$bm" == 7 ]
	then
		## read in FLAG cubes
                fits in=$combName'_cube.fits' op=xyin out=$combName'_cube.mir'

                ## set correct velocity frame in FLAG data cube
                velsw in=$combName'_cube.mir' axis='FELO'

                ## delete old FITS files
                rm $combName'_cube.fits'

                ## output new FITS files
                fits in=$combName'_cube.mir' out=$combName_cube'.fits' op=xyout
	fi
done
rm -rf $flagDataRoot'_cube.mir'
rm -rf $combName'_cube.mir'
