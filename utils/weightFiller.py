#!/opt/local/bin/ipython
"""
08/19/20
Unpacks binary files containing beamforming weights from FLAG (phased array feed for the GBT). The only user input is the path to
the directory containing the raw binary weight files. The weights will be output to a directory called: weight_files
User Inputs:
-p --path - <required> path to weight binary files

Usage: 
ipython weightFiller.py full_path_to_binaries
Example:
ipython weightFiller.py Users/npingel/Desktop/Research/data/GBT/GBT16B_400/AGBT16B_400_04/ 
__email__ = "Nickolas.Pingel@anu.edu.au"
__status__ = "Production"
"""

## imports 
import struct
import numpy as np
import os
import sys
import glob
import argparse
from astropy.io import fits

## make beam dictionaries to map from BYU to WVU convention (e.g. BYU 0 -> WVU 1)
wvuBeamDict = {0:1, 1:2, 2:6, 3:0, 4:3, 5:5, 6:4}
byuBeamDict = {1:0, 2:1, 6:2, 0:3, 3:4, 5:5, 4:6}

## get directory as command line argument. 
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", help = "<required> path to weight binary files (e.g., ./AGBT16B_400_NGC6946)", required = True)
args, unknown = parser.parse_known_args()

binDir = args.path

## global values 
floatBytes = 4
totalElements = 64
#numBeams = 1028
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

## ERROR HANDLING
## test the validity of the user provided data path
if not os.path.exists(binDir):
	print('Path to data directory is not valid; exiting...')
	sys.exit(1)

## check to see whether a backslash is required to add to path
finalChar = binDir[-1]
if not finalChar == '/':
	binDir+='/'

binList = glob.glob(binDir + '*.bin')

## if fitsList is empty, the provided time stamp is wrong
if len(binList) == 0:
	print('No weight files found in' + binDir + '; exiting...')
	sys.exit(1)

binList.sort()
binList = binList[0:20]
print('Using weight files: ')
for weightFile in binList:
	print(weightFile)

## create the weight_files directory to store the output FITS files. 
## check for existence
dirExist = os.path.isdir('./weight_files')
if not os.path.isdir('./weight_files'):
	os.makedirs('./weight_files')
	## inform user
	print('Creating weight_files directory to store weight FITS files...')
for i in range(0,len(binList)):
	progressBar(i,len(binList))
	file = binList[i]
	fileSplit = file.split('/')
	bank = file[-5]
	with open(file, 'rb') as f:
		data=f.read()

	## array to hold weights
	weightArr = np.zeros([numBeams*2,64*25*2],dtype='float32') ## rows: beams, columns: elements*frequency channels*complex pairs
	
	## total number of bytes in data payload
	totBytes = floatBytes*totalElements*numFreqChans*numBeams*pols*complexPairs
	
	## data payload in binary files (for a PAF calibration scan with course channels) are ordered such that we have correlations for beam0X, freq0, element 0, beam0X, freq0, element1...
	## beam0X, freq0, element 63, beam 0X, freq1, element0 ... beam0X, freq1, element 63 ... beam0X, freq24, element 63, beam1X, freq 0, element 0 ... beam1X, freq0, element0 ...
	## beam6X, freq 24, element 63, beam0Y, freq 0, element0 ... beam6Y, freq 24, element 63

	## based on the above structure, unpack the binary such that we split the data into individual polarizations and beams. First 7 rows are 
	## beams 0-6, XX Pol; last 7 rows are beams 0-6, YY Pol. Convert from BYU to WVU beam name convention.  
	bytesInBeam = floatBytes*totalElements*numFreqChans*complexPairs
	for plNum in range(0, 2): ## 0 for XX, 1 for YY
		polIdx = plNum * bytesInBeam * numBeams
		for byuBeam in range(0, numBeams):
			weightIdx = 0 
			for idx in range(0, bytesInBeam, 4):
				absIdx =  polIdx + (byuBeam * bytesInBeam) + idx 
				wvuBeam = wvuBeamDict[byuBeam]
				weightArr[wvuBeam+plNum*7, weightIdx] = struct.unpack('f',data[absIdx:absIdx + 4])[0]
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
    
    ## offsets are stored like el0, xel0, el1, xel1, ... 
	## array to store beam offsets. Rows: el and xel offset, columns (beams)
	offSetArr = np.zeros([2, numBeams],dtype='float32')

	for wvuBeam in range(0, numBeams):
		byuBeam = byuBeamDict[wvuBeam]
		offsetIdx = startByte+(byuBeam*8)
		offSetArr[0, wvuBeam] = struct.unpack('f',data[offsetIdx:offsetIdx+4])[0] ## in deg
		offSetArr[1, wvuBeam] = struct.unpack('f',data[offsetIdx+4:offsetIdx+8])[0]

	## update starting byte to get filenames and BF algorithm
	charStartByte = startByte+(numBeams*2*4)
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
	prihdr.set('EXTNAME', 'BF Info', 'name of this binary table extension')
	priHdu = fits.PrimaryHDU(header=prihdr)


	##Make table
     
	colList = []
	for col in range(0, numBeams):
		colList.append(fits.Column(name = 'Beam' + np.str(col) + 'X', format = '1E', array = weightArr[col,:]))
	for col in range(0, numBeams):
		colList.append(fits.Column(name = 'Beam' + np.str(col) + 'Y', format = '1E', array = weightArr[col + numBeams, :]))
	colList.append(fits.Column(name = 'BeamOff_XEL', format = '1E', array = offSetArr[1, :]))
	colList.append(fits.Column(name = 'BeamOff_EL', format = '1E', array = offSetArr[0, :]))

	tbHdu = fits.BinTableHDU.from_columns(colList)
	tbHdu.header.set('EXTNAME','Beam Weights and Offsets', 'name of this binary table extension')
	thdulist = fits.HDUList([priHdu, tbHdu])
	
	## generate filename from weights file name
	fileName = fileSplit[-1]
	weightFITSName = 'weight_files/' + fileName[:-4] +'.FITS'
	thdulist.writeto(weightFITSName)
## when finished, print newline character
print('\n')
