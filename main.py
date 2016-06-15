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

##TODO: How to distinguish between GPU threads, e.g. threads from 1st GPU contianing channels 0-4,100-104, and 5-9,105-109
##TODO:Extend sorting from 4 FITS files to 20 FITS files.
##TODO: how to reset the global data buffers?
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
    

def main():
    bf = BeamformingModule()
    md = MetaDataModule()
    dataBuff_X = np.zeros([10,25*20])   ##TODO: make mode independent  
    dataBuff_Y = np.zeros([10,25*20])   ##TODO: grab number of integrations from header?
    banks = {"A" : bankA,
             "B" : bankB,
             "C" : bankC,
             "D" : bankD,
}
## change working directory to project dir assumed to be first and only command-line argument. 
    os.chdir(str(sys.argv[1]))    
    print('Changing working directory to: '+str(sys.argv[1]))
    print('Building Primary HDU...')
    priHeader = md.contstructPriHDUHeader   
    for fitsName in glob.glob('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/*.fits'):    
        print('Beamforming correlations in: '+fitsName[-25:])
        xData,yData,nints = bf.getSpectralArray(fitsName)
        bank = fitsName[-6]
        bank = bank[0]##TODO: grab from header?        
        for ints in range(0,nints):
            banks.get(bank)(dataBuff_X,xData[ints,:],ints)
            banks.get(bank)(dataBuff_Y,yData[ints,:],ints)
        ##TODO:construct actual binTbl
    
    
        
    

#run main function
if __name__=="__main__":
    main()
