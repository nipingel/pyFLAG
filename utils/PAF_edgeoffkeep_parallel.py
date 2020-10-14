"""
02/13/20
A python wrapper to submit multiple processing jobs that will calibrate all FLAG SDFITS files.

User Inputs:

-p --path - <required> path to SDFITS files; only requires file prefix (e.g., ./AGBT16B_400_NGC6946)
-o --outFile - <required> name of output file prefixes (will have _edge.fits appended)
-s --sourceName - <required> name of source required to identify relevant scans
-x --sefdX - <required> list of XX SEFD values for each beam
-y --sefdY - <required> list of YY SEFD values for each beam 
-r --chanRange - <required> list of channels over which to fit a polynomial and subtract baseline
-f --order - <required> order of polynomial
-l --raDecLimits - <optional> list of long/latitude coordinates (in deg) that define a box containing edge emission. 
                   Scans with these coordinates will use opposte edge of map to construct off spectrum
-c --factor - <optional> list of scaling factors
-b --blank - <optional> set to True if smoothing is to be applied to blanked data
-n --noSmooth - <optional> set to True if smoothing is NOT to be applied to smoothed data"
-m --beamList -<optional> list of beams ot process; defauls from 0-6. MUST BE EQUAL TO NUMBER OF sefdX and sefdY values 

__email__ = "Nickolas.Pingel@anu.edu.au"
__status__ = "Production"
"""

## imports
import argparse
import sys
import os
import glob
import time
from multiprocessing import Pool

## define function to call PAF_edgeoffkeep_indv.pro to calibrate each beam
def calibrateBeam(fileName, srcStr, outFileName, sefdXVal, sefdYVal, order, chanRangeList, raDecLimitsList, fact, beam):
	## construct file name strings
	if args.blank == True:
		fileNameStr = '%s_Beam%s_blank.fits' % (fileName, beam) 
	elif args.noSmooth == True:
		fileNameStr = '%s_Beam%s.fits' % (fileName, beam)
	elif (args.noSmooth == True) and (args.blank == True):
		fileNameStr = '%s_Beam%s_blank.fits' % (fileName, beam)
	elif (args.noSmooth == False) and (args.blank == True):
		fileNameStr = '%s_Beam%s_ss_blank.fits' % (fileName, beam)
	else:
		fileNameStr = '%s_Beam%s_ss.fits' % (fileName, beam) 
	outFileStr = fileNameStr.replace('Beam%s' % beam, 'Beam%s_edge' % beam)
	
	## construct channel range, sefd, raDecLimit strings
	chanRangeStr = " ".join(chanRangeList)
	raDecLimitsStr = " ".join(raDecLimitsList)
       
	## make call to shell
	os.system("gbtidl -e 'edgeoffkeep_indv' -args %s %s %s %s %s %s %s %s %s" % (fileNameStr, srcStr, outFileStr, sefdXVal, sefdYVal, order, fact, chanRangeStr, raDecLimitsStr))

	print('Finished calibrating  %s' % fileNameStr)

## parse user arguments
parser = argparse.ArgumentParser()

parser.add_argument("-p", "--path", help = "<required> path to SDFITS files; only requires file prefix (e.g., ./AGBT16B_400_NGC6946)", required = True)
parser.add_argument("-o", "--outFile", help = "<required> name of output file prefixes (will have _edge.fits appended)", required = True)
parser.add_argument("-s", "--sourceName", help= "<required> name of source required to identify relevant scans", required = True)
parser.add_argument("-x", "--sefdX", help = "<required> list of XX SEFD values for each beam", required = True, nargs = "+")
parser.add_argument("-y", "--sefdY", help = "<required> list of YY SEFD values for each beam", required = True, nargs = "+")
parser.add_argument("-r", "--chanRange", help = "<required> list of channels over which to fit a polynomial and subtract baseline", required = True, nargs = "+")
parser.add_argument("-f", "--order", help = "<optional> order of polynomial")
parser.add_argument("-l", "--raDecLimits", help = "<optional> list of long/latitude coordinates (in deg) that define a box containing edge emission. Scans with these coordinates will use opposte edge of map to construct off spectrum", nargs = "+")
parser.add_argument("-c", "--factor", help = "<optional> list of scaling factors", nargs = "+")
parser.add_argument("-b", "--blank", help = "<optional> set to True if calibration is to be applied to blanked data", required = False, default = False, type = bool)
parser.add_argument("-n", "--noSmooth", help = "<optional> set to True if smoothing is NOT to be applied to smoothed data", required = False, default = True, type = bool)
parser.add_argument("-m", "--beamList", help = "<optional> list of beams ot process; defauls from 0-6. MUST BE EQUAL TO NUMBER OF sefdX and sefdY values", nargs = "+")

args, unknown = parser.parse_known_args()

filePrefix = args.path
outFile = args.outFile
sname = args.sourceName
sefdX = args.sefdX
sefdY = args.sefdY
order = args.order
chanRange = args.chanRange
raDecLimits = args.raDecLimits

## construct beam list
if args.beamList:
	beamList = args.beamList
if not args.beamList:
	beamList = [str(i) for i in range(0, 7)]
if args.raDecLimits:
	raDecLimits = args.raDecLimits
	if len(raDecLimits) != 4:
		print('Must provide four coordinates (in deg) for raDecLimits')
		sys.exit(1)
if not args.raDecLimits:
	raDecLimits = '0 0 0 0'
if args.factor:
	factorList = args.factor
if not args.factor:
	factorList = ['1.0', '1.0', '1.0', '1.0', '1.0', '1.0', '1.0']

## check that sefd values are equal to the number of beams requisted. If not, exit
if len(beamList) != len(sefdX) or len(beamList) != len(sefdY):
	print('Number of beams must be equal to number of SEFD values givin')
	sys.exit(1)
## pass to processing pool
p = Pool()

## pack up the iterable for function
filePrefixList = [filePrefix for i in range(0, len(beamList))]
outFileList  = [outFile for i in range(0, len(beamList))]
snameList = [sname for i in range(0, len(beamList))]
orderList = [order for i in range(0, len(beamList))]
chanRangeList = [chanRange for i in range(0, len(beamList))]
raDecLimitsList = [raDecLimits for i in range(0, len(beamList))]

## prepare input to processing method
processList = list(zip(filePrefixList, snameList, outFileList, sefdX, sefdY, orderList, chanRangeList, raDecLimitsList, factorList, beamList))

p.starmap(calibrateBeam, processList)

p.close()
p.join()

