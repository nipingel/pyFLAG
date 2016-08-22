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
def bankA(dataBuff,data,ints):
    endChan = 0
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan] = data[i:i+5]
        endChan+=95
def bankB(dataBuff,data,ints):
    endChan = 5
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan] = data[i:i+5]
        endChan+=95
 
def bankC(dataBuff,data,ints):
    endChan = 10
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan] = data[i:i+5]
        endChan+=95 
 
def bankD(dataBuff,data,ints):
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
    banks = {"A" : bankA,
             "B" : bankB,
             "C" : bankC,
             "D" : bankD,
             }
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
        for beam in range(0,7):
            globalDataBuff_X = np.zeros([len(fileList[0]),numInts,25*20]) ##TODO: make general based on CovMode
            globalDataBuff_Y = np.zeros([len(fileList[0]),numInts,25*20]) ##TODO: make general based on CovMode
            cnt = 0
            fileIdx = 0
            for file in fileList:        
                print('Beamforming correlations in: '+file[-25:]+', Beam: '+np.str(beam)) 
                if cnt % numTotalThreads == 0:
                    dataBuff_X = np.zeros([numInts,25*20])   ##TODO: make mode independent  
                    dataBuff_Y = np.zeros([numInts,25*20])        
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
                cnt+=1
                fileIdx+=1
                globalDataBuff_X[fileIdx,:,:] = dataBuff_X
                globalDataBuff_Y[fileIdx,:,:] = dataBuff_Y
            md = MetaDataModule(fitsList[objs],numInts,globalDataBuff_X,globalDataBuff_Y,beam)
            md.constuctBinTableHeader()
        
    freqChans = np.linspace(1,500,500)
    
    
    pyplot.ylabel('Flux Density [Jy]', fontsize=18)
    pyplot.xlabel('Freq. Channels', fontsize=18)
    pyplot.title('Leo A', fontsize=18)
    pyplot.plot(freqChans,dataBuff_X[0], color='red', label='X-Pol',linewidth=1.5)
    pyplot.plot(freqChans,dataBuff_Y[0], color='blue', label='Y-Pol',linewidth=1.5)
    pyplot.legend(loc=0,prop={'size':14},fontsize='large')
#    pyplot.savefig('/Users/npingel/Desktop/FLAG_LeoA_ExtendedFISHFITS.pdf', bbox_inches='tight')
    pyplot.show()
    ##TODO:construct actual binTbl
    ##Get objects observed and their corresponding GO FITS file names.
    objList,fitsList = numObjs()
    numInts=50 ##TODO: update and how to figure this out during data processing? Basically number of integrations*numPol
    for objs in range(0,len(objList)):
        md = MetaDataModule(fitsList[objs],numInts)
        md.constuctBinTableHeader()
        
    
        
    

#run main function
if __name__=="__main__":
    main()
