"""
01/27/20
A python wrapper to submit multiple processing jobs that will smooth all FLAG beam SDFITS files to a user specified
resolution in parallel.

User Inputs:

-p --path - <required> path to SDFITS files; only requires file prefix (e.g., ./AGBT16B_400_NGC6946_Beam)
-o --outFile - <required> name of output file prefixes (will have _ss.fits appended)
-s --sourceName - <required> name of source required to identify relevant scans
-k --kernel - <required> list of kernel parameters (e.g., 0.28 1 1 0.28)

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

## TODO: unpack user arguments
filePrefix = args.path
outFile = args.outFile
sname = args.sourceName

## TODO: get files 

## pass to processing pool
p = Pool()
result = p.starmap_async()