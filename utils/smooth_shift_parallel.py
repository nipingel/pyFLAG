"""
01/27/20
A python wrapper to submit multiple processing jobs that will smooth all FLAG beam SDFITS files to a user specified
resolution in parallel.

User Inputs:

-p --path - <required> path to SDFITS files; only requires file prefix (e.g., ./AGBT16B_400_NGC6946_Beam)
-o --outFile - <required> name of output file prefixes (will have _ss.fits appended)
-s --sourceName - <required> name of source required to identify relevant scans
-k --kernel - <required> list of kernel parameters (e.g., 0.28 1 1 0.28)
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

## define function to call smooth_shift_indv.pro for each beam
def smoothBeam(fileName, srcStr, outFileName, beam, kern):
	## construct file name str
	fileNameStr = '%s_Beam%s.fits' % (fileName, beam) 
	outFileStr = fileNameStr.replace('.fits', '_ss.fits')
	kernStr = " ".join(kern)
	## make call to shell
	os.system("gbtidl -e 'smooth_shift_indv' -args %s %s %s %s" % (fileNameStr, srcStr, outFileStr, kernStr))

	print('Finished smoothing for %s' % fileNameStr)

## TODO: unpack user arguments
parser = argparse.ArgumentParser()

parser.add_argument("-p", "--path", help = "<required> path to SDFITS files; only requires file prefix (e.g., ./AGBT16B_400_NGC6946_Beam)", required = True)
parser.add_argument("-o", "--outFile", help = "<required> name of output file prefixes (will have _ss.fits appended)")
parser.add_argument("-s", "--sourceName", help= "<required> name of source required to identify relevant scans")
parser.add_argument("-k", "--kernel", help = "<required> list of kernel parameters (e.g., 0.28 1 1 0.28)", nargs = '+')
parser.add_argument("-m", "--beamList",  help = "<optional> list of beams ot process; defauls from 0-6.", nargs = '+')

args, unknown = parser.parse_known_args()

filePrefix = args.path
outFile = args.outFile
sname = args.sourceName
kernel = args.kernel

if args.beamList:
	beamList = args.beamList
if not args.beamList:
	beamList = [str(i) for i in range(0, 7)]

## pass to processing pool
p = Pool()

## pack up the iterable for function
filePrefixList = [filePrefix for i in range(0, len(beamList))]
outFileList  = [outFile for i in range(0, len(beamList))]
snameList = [sname for i in range(0, len(beamList))]
kernelList = [kernel for i in range(0, len(beamList))]
processList = list(zip(filePrefixList, snameList, outFileList, beamList, kernelList))

result = p.starmap(smoothBeam, processList)

p.close()
p.join()