# -*- coding: utf-8 -*-
"""
Created on Mon May 23 13:40:47 2016
This is the main function for the FLAG spectral line filler. It imports the I/O functionality of astropy,
numpy, and matplotlib. It first gathers the necessary metadata from the GBT ancillary FITS files (e.g. Antenna, GO). 
It collates these metadata into a new FITS classification 'BfSpec', which will be GBTIDL friendly

@author: npingel
"""
#imports
from astropy.io import fits
import numpy as np
import sys
import os
import glob
from modules.metaDataModule import MetaDataModule
from modules.beamformingModule import BeamformingModule
import matplotlib.pyplot as pyplot

##global variables
numGPU = 2
numTotalThreads = numGPU*2

##TODO:Extend sorting from 4 FITS files to 20 FITS files.
def xid0(dataBuff,data,ints):
    endChan = 0
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan] = data[i:i+5]
        endChan+=95
def xid1(dataBuff,data,ints):
    endChan = 5
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan] = data[i:i+5]
        endChan+=95
 
def xid2(dataBuff,data,ints):
    endChan = 10
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan] = data[i:i+5]
        endChan+=95 
 
def xid3(dataBuff,data,ints):
    endChan = 15
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan] = data[i:i+5]
        endChan+=95   

##function to determine number of objects observed in session
def numObjs():        
    goFits = glob.glob('/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/GO/*.fits') ##TODO: change to work on flag03/lustre
    objList = []
    fitsList = []
    itr = 0.
    for goName in goFits:
       goHDU = fits.open(goName)
       if itr == 0:
           objList.append(goHDU[0].header['OBJECT'])
           fitsList.append([])
           fitsList[-1].append(goName[-24:])
           itr+=1
       else:
           obj = goHDU[0].header['OBJECT']
           if objList[-1] != obj:
               objList.append(goHDU[0].header['OBJECT']) 
               fitsList.append([])
           fitsList[-1].extend([goName[-24:]])
    return objList, fitsList
def main():
    bf = BeamformingModule()
    banks = {"A" : xid0,
             "B" : xid3,
             "C" : xid1,
             "D" : xid2,
             }
##TODO: include intLen parameter
"""
## change working directory to project dir assumed to be first and only command-line argument. 
    os.chdir(str(sys.argv[1]))    
    print('Changing working directory to: '+str(sys.argv[1]))
    print('Building Primary HDU...')
    #priHeader = md.contstructPriHDUHeader   
    os.chdir('/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/RawData/')
    objList,fitsList = numObjs()    
    ##TODO: put back in: for objs in range(0,len(objList)):
    for objs in range(1,2):
        numInts =4806 ##TODO: grab number of integrations from header?          
        fileList = fitsList[objs]  
        fileList = ['2016_07_29_10:46:06A.fits','2016_07_29_10:46:06B.fits','2016_07_29_10:46:06C.fits','2016_07_29_10:46:06D.fits'] ##TODO: make general
        for beam in range(0,7): ##TODO: change back to 7 
            globalDataBuff_X = np.zeros([int(len(fileList)/len(banks)),numInts,25*20]) ##TODO: make general based on CovMode
            globalDataBuff_Y = np.zeros([int(len(fileList)/len(banks)),numInts,25*20]) ##TODO: make general based on CovMode
            cnt = 0
            fileIdx = -1
            for file in fileList:        
                print('Beamforming correlations in: '+file[-25:]+', Beam: '+np.str(beam)) 
                if cnt % numTotalThreads == 0:
                    dataBuff_X = np.zeros([numInts,25*20])   ##TODO: make mode independent  
                    dataBuff_Y = np.zeros([numInts,25*20])
                    fileIdx+=1
                bank = file[-6]
                bank = bank[0]##TODO: grab from header? 
                if bank == 'A':
                    xid = 0
                elif bank == 'B':
                    xid = 3
                elif bank == 'C':
                    xid = 1
                elif bank == 'D':
                    xid = 2
                xData,yData,nints = bf.getSpectralArray(file,beam,xid)       
                for ints in range(0,nints):
                    banks.get(bank)(dataBuff_X,xData[ints,:],ints)
                    banks.get(bank)(dataBuff_Y,yData[ints,:],ints)
                globalDataBuff_X[fileIdx,:,:] = dataBuff_X
                globalDataBuff_Y[fileIdx,:,:] = dataBuff_Y
                cnt+=1
                ##TODO:REMOVE ONCE METADATA MODULE IS TESTED
            fits.writeto('/Users/npingel/Desktop/BFSpec_NGC6946_beam'+np.str(beam)+'_X.fits',globalDataBuff_X)
            fits.writeto('/Users/npingel/Desktop/BFSpec_NGC6946_beam'+np.str(beam)+'_Y.fits',globalDataBuff_Y)
"""
hdu = fits.open('/Users/npingel/Desktop/BFSpec_NGC6946_beam0_X.fits')
globalDataBuff_X = hdu[0].data
hdu = fits.open('/Users/npingel/Desktop/BFSpec_NGC6946_beam0_Y.fits')
beam = 0
intLen = .25 ##seconds
numInts = 4806
globalDataBuff_Y = hdu[0].data
fileList = ['2016_07_29_10:46:06A.fits','2016_07_29_10:46:06B.fits','2016_07_29_10:46:06C.fits','2016_07_29_10:46:06D.fits']
md = MetaDataModule(fileList,numInts,globalDataBuff_X,globalDataBuff_Y,beam,intLen,numTotalThreads)
md.constuctBinTableHeader()
"""    
    freqChans = np.linspace(1,500,500)
    
    
    pyplot.ylabel('Flux Density [Jy]', fontsize=18)
    pyplot.xlabel('Freq. Channels', fontsize=18)
    pyplot.title('Time Ave BF Spectrum of NGC4646 (Beam 0)', fontsize=18)
    #pyplot.plot(freqChans,np.mean(dataBuff_X, axis=0), color='red', label='X-Pol',linewidth=1.5)
    pyplot.plot(freqChans,np.mean(dataBuff_Y, axis=0), color='blue', label='Y-Pol',linewidth=1.5)
    pyplot.legend(loc=0,prop={'size':14},fontsize='large')
    pyplot.savefig('/Users/npingel/Desktop/AveBFSpec_NGC6946_Beam0.pdf', bbox_inches='tight')
    pyplot.show()
    ##TODO:construct actual binTbl
    ##Get objects observed and their corresponding GO FITS file names.
#    objList,fitsList = numObjs()
#    numInts=50 ##TODO: update and how to figure this out during data processing? Basically number of integrations*numPol
#    for objs in range(0,len(objList)):
#        md = MetaDataModule(fitsList[objs],numInts)
       #md.constuctBinTableHeader()
"""        
    
        
    

#run main function
if __name__=="__main__":
    main()
