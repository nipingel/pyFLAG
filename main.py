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

##TODO: Makeshift case code to sort beamformed channels to be contiguous
options = {A : bankA,
                B : bankB ,
                C : bankC,
                D : bankD,
}
 
def bankA():
    print "You typed zero.\n"
 
def bankB():
    print "n is a perfect square\n"
 
def bankC():
    print "n is an even number\n"
 
def bankD():
    print "n is a prime number\n"
    

def main():
    bf = BeamformingModule()
    md = MetaDataModule()
    
## change working directory to project dir assumed to be first and only command-line argument. 
    os.chdir(str(sys.argv[1]))    
    print('Changing working directory to: '+str(sys.argv[1]))
    print('Building Primary HDU...')
    priHeader = md.contstructPriHDUHeader
##TODO:Decide on data processing algorithm     
    for fitsName in glob.glob('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/*.fits'):    
        dataBuff = np.zeros([25*20])   ##TODO: make mode independent     
        print('Beamforming correlations in: '+fitsName[-25:])
        xData,yData = bf.getSpectralArray(fitsName)
        
    
        
    

#run main function
if __name__=="__main__":
    main()
