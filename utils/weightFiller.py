"""
5/29/17
Unpacks binary files containing beamforming weights from FLAG (phased array feed for the GBT). One command line argument is required
which is the directory to the binary weight files. The second argument is either 'FRB', 'Spectral', or 'Calibration'
Usage: 
ipython weightFiller.py full_path_to_binaries CorrMode
Example:
ipython weightFiller.py Users/npingel/Desktop/Research/data/GBT/GBT16B_400/AGBT16B_400_04/weight_files/ 'Calibration'
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "nipingel@mix.wvu.edu"
__status__ = "Production"
"""

## imports 
import struct
import numpy as np
import sys
import glob

from astropy.io import fits


## get directory as command line argument. 
binDir = sys.argv[1]

## global values 
floatBytes = 4
totalElements = 64
numBeams = 7
pols = 2
complexPairs = 2 
numFreqChans = 25

## progress bar
def progressBar(value, endvalue,bar_length=20):

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent of weight files processed: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()

# TODO: input exceptions

binList = glob.glob(binDir + '*.bfw')

for i in range(0,len(binList)):
	progressBar(i,len(binList))
	file = binList[i]
	fileSplit = file.split('/')
	bank = file[-5]
        with open(file, 'rb') as f:
		data=f.read()

	## array to hold weights
	weightArr = np.zeros([14,64*25*2],dtype='float32') ## rows: beams, columns: elements*frequency channels*complex pairs
	
	## total number of bytes in data payload
	totBytes = floatBytes*totalElements*numFreqChans*numBeams*pols*complexPairs
	
	## data payload in binary files (for a PAF calibration scan with course channels) are ordered such that we have correlations for beam0XX, freq0, element 0, beam0X, freq0, element1...
	## beam0X, freq0, element 1 ... beam 0X, freq1, element63, beam0X, freq1, element 0 ... beam0X, freq1, element 63, ... beam0X, freq 24, element 63 
	## ... beam1X, freq0, element0, ... beam6X, freq 24, element 63, beam0Y, freq 0, element0 ... beam6Y, freq 24, element 63

	## based on the above structure, unpack the binary such that we split the data into individual polarizations and beams. First 7 rows are 
	## beams 0-6, XX Pol; last 7 rows are beams 0-6, YY Pol. 
	bytesInBeam = floatBytes*totalElements*numFreqChans*complexPairs
	for beam in range(0,14):
		weightIdx = 0 
		for idx in range(0,bytesInBeam,4):
			absIdx = (beam*bytesInBeam)+idx 
			weightArr[beam,weightIdx] = struct.unpack('f',data[absIdx:absIdx+4])[0]
			weightIdx+=1

	## now, let's unpack the meta data. This portion of the binary file will always have the same number of bytes (4160) with format:
	## Beam 0 offsets (2 floats [arcmin])
	## Beam 1 offsets (2 floats [arcmin])
	## Beam 2 offsets (2 floats [arcmin])
	## Beam 3 offsets (2 floats [arcmin])
	## Beam 4 offsets (2 floats [arcmin])
	## Beam 5 offsets (2 floats [arcmin])
	## Beam 6 offsets (2 floats [arcmin])
	## calibration Set filename [64 char]
	## beamforming algorithm [64 char]
	## X-Engine ID (unit64)

	##Metadata
	startByte = totBytes
	## array to store beam offsets. Rows: Az and El offset, columns (beams)
	offSetArr = np.zeros([2,7],dtype='float32')
	for off in range(0,7):
	    offsetIdx = startByte+(off*8)
            offSetArr[0, off] = struct.unpack('f',data[offsetIdx:offsetIdx+4])[0]
	    offSetArr[1, off] = struct.unpack('f',data[offsetIdx+4:offsetIdx+8])[0]

	## update starting byte to get filenames and BF algorithm
	charStartByte = startByte+(14*4)
	calibFilename = struct.unpack('64c',data[charStartByte:charStartByte+64])
	beamformMethod = struct.unpack('64c',data[charStartByte+64:charStartByte+128])

	## decode the filename and method to strings
	calibFilename_Str = ''
	beamformMethod_Str = ''
	for byteIdx in range(len(calibFilename)):
		val = calibFilename[byteIdx].decode('utf-8')
		if val != '\x00':
			calibFilename_Str+=val
	for byteIdx in range(len(beamformMethod)):
		val = beamformMethod[byteIdx].decode('utf-8')
		if val != '\x00':
			beamformMethod_Str+=val
	## last eight bytes gives X-engine ID
	xEngineID = struct.unpack('Q',data[-8:])[0]

	##Create FITS file
	prihdr = fits.Header()
	##TODO: add header parameters
	prihdr = fits.Header()
	prihdr.set('CALFILE',calibFilename_Str, 'Calibration Set Filename')
	prihdr.set('BEAMFORM',beamformMethod_Str,'Beamforming algorithm used')
	prihdr.set('XENGINE',xEngineID,'X-Engine ID')
	priHdu = fits.PrimaryHDU(header=prihdr)


	##Make table 
	col1 = fits.Column(name='Beam0X', format='1E', array=weightArr[0,:])
	col2 = fits.Column(name='Beam1X', format='1E', array=weightArr[1,:])
	col3 = fits.Column(name='Beam2X', format='1E', array=weightArr[2,:])
	col4 = fits.Column(name='Beam3X', format='1E', array=weightArr[3,:])
	col5 = fits.Column(name='Beam4X', format='1E', array=weightArr[4,:])
	col6 = fits.Column(name='Beam5X', format='1E', array=weightArr[5,:])
	col7 = fits.Column(name='Beam6X', format='1E', array=weightArr[6,:])
	col8 = fits.Column(name='Beam0Y', format='1E', array=weightArr[7,:])
	col9 = fits.Column(name='Beam1Y', format='1E', array=weightArr[8,:])
	col10 = fits.Column(name='Beam2Y', format='1E', array=weightArr[9,:])
	col11 = fits.Column(name='Beam3Y', format='1E', array=weightArr[10,:])
	col12 = fits.Column(name='Beam4Y', format='1E', array=weightArr[11,:])
	col13 = fits.Column(name='Beam5Y', format='1E', array=weightArr[12,:])
	col14 = fits.Column(name='Beam6Y', format='1E', array=weightArr[13,:])
	col15 = fits.Column(name='BeamOff_AZ', format='1E', array=offSetArr[0,:])
	col16 = fits.Column(name='BeamOff_EL', format='1E', array=offSetArr[1,:])
	cols = fits.ColDefs([col1, col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16])
	tbHdu = fits.new_table(cols)
	tbHdu.header.set('EXTNAME','BF Weights', 'name of this binary table extension')
	thdulist = fits.HDUList([priHdu, tbHdu])
	
	## generate filename from split path
        weightFITSName = binDir + '/weight_files/w_' + fileSplit[-3] + '_' + bank +'.FITS'
	thdulist.writeto(weightFITSName)
## when finished, print newline character
print('\n')
