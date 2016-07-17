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
def bankA(dataBuff,data,ints,beam):
    endChan = 0
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan,beam] = data[i:i+5]
        endChan+=95
def bankB(dataBuff,data,ints,beam):
    endChan = 5
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan,beam] = data[i:i+5]
        endChan+=95
 
def bankC(dataBuff,data,ints,beam):
    endChan = 10
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan,beam] = data[i:i+5]
        endChan+=95 
 
def bankD(dataBuff,data,ints,beam):
    endChan = 15
    for i in range(0,25,5):
        startChan = endChan
        endChan = startChan+5       
        dataBuff[ints,startChan:endChan,beam] = data[i:i+5]
        endChan+=95   

##function to determine number of objects observed in session
def numObjs():        
    goFits = glob.glob('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/GO/*.fits') ##TODO: change to work on flag03/lustre
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
    fitsNames = glob.glob('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/*.fits')
    objList,fitsList = numObjs()    
    for objs in range(0,len(objList)):
        cnt = 0
        numInts =  10        
        file = fitsList[objs]
        numInts = 10 ##TODO: grab number of integrations from header?         
        for beam in range(0,7):         
            print('Beamforming correlations in: '+file[-25:]+', Beam: '+np.str(beam)) 
            if cnt % numTotalThreads == 0:
                dataBuff_X = np.zeros([numInts,25*20,7])   ##TODO: make mode independent  
                dataBuff_Y = np.zeros([numInts,25*20,7])        
            xData,yData,nints = bf.getSpectralArray(file,beam)
            bank = file[-6]
            bank = bank[0]##TODO: grab from header?        
            for ints in range(0,nints):
                banks.get(bank)(dataBuff_X,xData[ints,:],ints,beam)
                banks.get(bank)(dataBuff_Y,yData[ints,:],ints,beam)
        cnt+=1
        md = MetaDataModule(fitsList[objs],numInts,dataBuff_X,dataBuff_Y)
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
