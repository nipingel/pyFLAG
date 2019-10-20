# -*- coding: utf-8 -*-
"""
This script generates a 1D vector where each element's value is the index of correlation pair in the GPU output
data buffer. For example, if the fourth correlation in the "map" vector has a value of 20, it means the fourth correlation in the
FITS file corresponds to the 20th correlation pair in the GPU data buffer. 
Usage:
ipython gpuToCovarMatrix.py
Example:
ipython gpuToCovarMatrix.py
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "Nickolas.Pingel@anu.edu.au"
__status__ = "Production"
"""

#imports
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

##The 1D output data buffer from the GPU is a 1D vector of type complex64 with a length of 2112*NUM_CHANS,
##NUM_CHANS differes depending on mode. The order of correlation pairs (real, imaginary) is mode independent. 

##We actually rearrange each 2112 correlation pair block (one frequency channel) of the GPU output data buffer 
##back into it's 2D representation. That way, going to FISHFITS:
# The data order in the 'cross_accum' buffer is ch1Xch1_freq1_real,
# ch1Xch1_freq1_imag, ch1Xch2_freq1_real... ch1Xch20_freq1_imag,
# ch2Xch1_freq1_real...ch2Xch20_freq1_imag, ch3Xch1_freq1_real...
# ch1Xch1_freq2_real, ch1Xch1_freq2_imag...ch20Xch20_freq2_imag,
# ch1Xch1_freq3_real...ch20Xch20_freqN_imag
## is simply grabbing each column while decreasing the number of elements grabbed in each subsequent column by one to avoid zeros
## and redundant correlations. 

##Since each chunk of four correlation pairs constitute a block in the big correlation matrix, we can simply fill out a 2D array
corrIdxMatrix = np.zeros([40,40], dtype='complex64')
gpuIdxVector = np.arange(0,2112,1)
lastBlcIdx = 0
colOffset = 1
for row in range(0,40,2):
    frstBlcIdx = lastBlcIdx
    lastBlcIdx = frstBlcIdx+(colOffset*4)
    fullDataBlock = gpuIdxVector[frstBlcIdx:lastBlcIdx]
    subColOffset = 0
    for col in range(0,row+2):
        if (col % 2 == 0):
           frstSubBlkIdx = subColOffset*4
           lastSubBlkIdx = frstSubBlkIdx+4
           subBlock = fullDataBlock[frstSubBlkIdx:lastSubBlkIdx]
           corrIdxMatrix[row,col] = subBlock[0]
           corrIdxMatrix[row,col+1] = subBlock[1]
           corrIdxMatrix[row+1, col] = subBlock[2]
           corrIdxMatrix[row+1,col+1] = subBlock[3]
           subColOffset+=1
    colOffset+=1

for i in range(0,40):
  print('Data Channel: %s, Index: %d' % (i+1, corrIdxMatrix[i,i].real))
##Now, we must append subsequent columns to get correct FISHFITS order (while decreasing the 
##column index each time to avoid redundant correlation pairs. 
mapVector = []
strtOffset = 0
for col in range(0,40):
    mapVector.extend(np.real(corrIdxMatrix[strtOffset:40,col]))
    strtOffset+=1
np.savetxt('gpuToNativeMap.dat',mapVector,fmt='%d',delimiter='\n')
    
##Let's test the order by constructing the correlation pair (real+imag) for indices:
corrPairMatrix = np.zeros([40,40], dtype='complex64')
rowStrtOffset =  0
val = np.complex(0,0)
for col in range(0,40):
    for row in range(rowStrtOffset,40):
        val = np.complex(row+1,col+1)
        corrPairMatrix[row,col] = val
    rowStrtOffset+=1

