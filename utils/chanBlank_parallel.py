"""
02/02/20
A python wrapper to submit multiple processing jobs that will blank channels affected by PFB 'scalloping'
in a series of FLAG SDFITS files.

User Inputs:

-p --path - <required> path to SDFITS files; only requires file prefix (e.g., ./AGBT16B_400_NGC6946_Beam)
-o --outFile - <required> name of output file prefixes (will have _blank.fits appended)
-s --sourceName - <required> name of source required to identify relevant scans
-m --beamList -<optional> list of beams ot process; defauls from 0-6. 

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

## define function to call chanBlank_indv.pro for each beam
def smoothBeam(fileName, srcStr, outFileName, beam):
	
	## construct file name str
	fileNameStr = '%s_Beam%s.fits' % (fileName, beam) 
	outFileStr = fileNameStr.replace('.fits', '_blank.fits')
	## make call to shell
	os.system("gbtidl -e 'chanBlank_indv' -args %s %s %s" % (fileNameStr, srcStr, outFileStr))

	print('Finished blanking %s' % fileNameStr)


## unpack user arguments
parser = argparse.ArgumentParser()

parser.add_argument("-p", "--path", help = "<required> path to SDFITS files; only requires file prefix (e.g., ./AGBT16B_400_NGC6946_Beam)", required = True)
parser.add_argument("-o", "--outFile", help = "<required> name of output file prefixes (will have _ss.fits appended)", required = True)
parser.add_argument("-s", "--sourceName", help= "<required> name of source required to identify relevant scans", required = True)
parser.add_argument("-m", "--beamList",  help = "<optional> list of beams ot process; defauls from 0-6.", nargs = '+')

args, unknown = parser.parse_known_args()

filePrefix = args.path
outFile = args.outFile
sname = args.sourceName

if args.beamList:
	beamList = args.beamList
if not args.beamList:
	beamList = [str(i) for i in range(0, 7)]

# pass to processing pool
p = Pool()

## pack up the iterable for function
filePrefixList = [filePrefix for i in range(0, len(beamList))]
outFileList  = [outFile for i in range(0, len(beamList))]
snameList = [sname for i in range(0, len(beamList))]
processList = list(zip(filePrefixList, snameList, outFileList, beamList)

p.starmap(blankBeam, processList)

p.close()
p.join()

