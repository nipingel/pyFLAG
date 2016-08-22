# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 20:41:39 2016

@author: npingel
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


#hdu = fits.open('/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/RawData/2016_07_29_11:59:12A_Mod.fits')
hdu = fits.open('/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_04/RawData/2016_07_30_23:23:45F.fits')
data = hdu[1].data['DATA']
##Select time series of desired frequency channel (e.g. channel 204 is the 10th channel from BANKA)
chanStrt = 2112*7
chanEnd = chanStrt+2112
#chanStrt = 0
#chanEnd = 2112
freqChanTimeSeries = data[:,chanStrt:chanEnd]

tHot = 45.7
tCold = 6

##Now, contruct a 2D matrix in which the real element of the row,column entry is the index of the corresponding complex pair 
##in the 1D GPU output data vector
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
    
##Construct acutal corr Matrix from GBT
hotCorrMatrix = np.zeros([40,40], dtype='complex64')
coldCorrMatrix = np.zeros([40,40], dtype='complex64')
rowStrtOffset =  0
for col in range(0,40):
    for row in range(rowStrtOffset,40):
        idx = corrIdxMatrix[row,col].real
        hotCorrMatrix[row,col] = freqChanTimeSeries[1630,np.int(idx)]
        coldCorrMatrix[row,col] = freqChanTimeSeries[4004,np.int(idx)]
    rowStrtOffset+=1

    
Rhot = hotCorrMatrix+hotCorrMatrix.conj().T-np.diag(hotCorrMatrix)
Rcold = coldCorrMatrix+coldCorrMatrix.conj().T-np.diag(coldCorrMatrix)

hotPowerArr = np.zeros([40],dtype='float32')
coldPowerArr =  np.zeros([40],dtype='float32')

for i in range(0,40):
    hotPowerArr[i] = hotCorrMatrix[i,i].real
    coldPowerArr[i] = coldCorrMatrix[i,i].real
    
Y_factorArr = hotPowerArr/coldPowerArr
TsysArr = (tHot - (Y_factorArr*tCold))/(Y_factorArr-1)

TsysArrDerived = tHot/(Y_factorArr-1)
"""
plt.figure()
plt.plot(TsysArrDerived,color='blue')
plt.plot(TsysArr,color='red')
plt.show()
"""
sigList = []
for i in range(0,52800,2112):
    freqChanTimeSeries = data[:,i:i+2112]
    print('Analyizing channel: '+np.str(i/2112+1))
    maxVal = np.nanmax(freqChanTimeSeries[:,0].real)
    meanVal = np.nanmean(freqChanTimeSeries[4500:4800,0].real)
    signal = np.real(maxVal-meanVal)
    sigList.append(signal)
sigArr = np.array(sigList, dtype='float32')
maxChan = np.array(np.where(sigArr == np.max(sigArr)))
maxChan = maxChan.item(0)
freqChan = data[:,maxChan*2112:maxChan*2112+2112]
Y = np.nanmax(freqChan[:,0].real)/np.nanmean(freqChan[4500:4800,0].real)
tHot = 45.4
tCold = 3.3*np.exp(-0.01)/np.cos(np.deg2rad(14.6))
Tsys = (tHot-(Y*tCold))/(Y-1)

print("Frequency Channel w/ Max Signal:"+np.str(maxChan))
print("Y-factor for Freq chan and central element: "+np.str(Y))
print("Tsys: "+np.str(Tsys))
print("Peak Signal "+np.str(sigArr[maxChan])) 
plt.figure()
freqChannels = np.linspace(1,25,25)
plt.plot(freqChannels,sigList)
plt.show()

"""
data = hdu[1].data['DATA']
row = data[120]
startIdx = 0
corrMatrix = np.zeros([40,40],dtype='complex64')
colStrtIdx = 0
rowStrtIdx = 0
offset=40
row=row[0:820]
for i in range(0,40):
    column = row[rowStrtIdx:rowStrtIdx+offset]    
    corrMatrix[colStrtIdx:40,i] = column
    rowStrtIdx+=offset    
    colStrtIdx+=1
    offset-=1

    

R = corrMatrix+corrMatrix.conj().T-np.diag(corrMatrix)
        
      
    

plt.figure()
plt.imshow(np.abs(R))
plt.show()
"""
