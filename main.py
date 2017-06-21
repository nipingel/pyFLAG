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

## command line inputs
##TODO: exception handling
projectPath = sys.argv[1] ## of the form /home/gbtdata/AGBT16B_400_01
## split project path to get project string
projectPathSplit = projectPath.split('/')
projectStr = projectPathSplit[-1] + '/'
dataPath = '/lustre/projects/flag/' +  projectStr + 'BF/'

##global variables
numGPU = 10
numBanks = numGPU*2
##paths to ancillary FITS files
goFitsPath = projectPath + '/GO'

##TODO:Extend sorting from 4 FITS files to 20 FITS files.

## function to sort individual BANK data
## to the full bandpass based on XID. 
##TODO: add FRB mode functionality
def bandpassSort(xID, dataBuff, bankData, ints):
    ## which correlation mode are we in?
    ## determine this based on number of channels 
    ## in bandpasss 
    if len(dataBuff) == 500:
        ## coarse channel mode
        ## position in bandpass dictated by xID
        bandpassStartChan = xID*5
        bandpassEndChan = bandpassStartChan+5
        ## get each chunk of five contigious channels from BANK data,
        ## and place it in the proper spot in full bandpass
        for i in range(0, 5):
            bankStartChan = i * 5
            bankEndChan = bankStartChan + 5
            dataBuff[ints, bandpassStartChan:bandpassEndChan] = bankData[ints, startBankChan:endBankChan]
            
            ## increment bandpassStartChan/bandpassEndChan by 100 for proper position in full bandpass
            bandpassStartChan += 100
            bandpassEndChan = bandpassStartChan+5

##function to determine number of objects observed in session
def numObjs():        
    goFits = glob.glob(goFitsPath+'/*.fits')
    ## sort to get correct time stamp order
    goFitsSorted = np.sort(goFits)
    objList = []
    fitsList = []
    itr = 0.
    for goName in goFitsSorted:
       goHDU = fits.open(goName)
       if itr == 0:
           objList.append(goHDU[0].header['OBJECT'])
           fitsList.append([])
           fitsList[0].append(goName[-24:])
           itr+=1
       else:
           obj = goHDU[0].header['OBJECT']
           if obj not in objList:
               objList.append(goHDU[0].header['OBJECT']) 
               fitsList.append([])    
           fitsList[-1].append(goName[-24:])
    return objList, fitsList

##function to get number of integrations, integration length, and number of channels
def getScanInfo(fileName):
    fitsLst = glob.glob(dataPath + fileName[:-6] + '*.fits')
    hdu = fits.open(fitsLst[0])
    intLen = np.float(hdu[0].header['REQSTI'])
    numInts = np.int(hdu[1].header['NAXIS2'])
    form = np.float(hdu[1].header['TFORM3'][:-1])
    numChans = form/2112
    return numInts, intLen, numChans, np.sort(fitsLst)

def main():
    bf = BeamformingModule(dataPath)
    bankDict = {"A" : 0,
             "B" : 1,
             "C" : 2,
             "D" : 3,
	     "E" : 4, 
	     "F" : 5, 
             "G" : 6, 
             "H" : 7,
             "I" : 8, 
             "J" : 9,
             "K" : 10,
             "L" : 11,
             "M" : 12,
             "N" : 13,
             "O" : 14,
             "P" : 15,
             "Q" : 16,
             "R" : 17,
             "S" : 18,
             "T" : 19,
             }
##TODO: check current directory permissions -- must have writing access
## the only commandline argument should be path to project/session ancillary FITS files    
    pwd = os.getcwd()
    print('Project directory: ' + np.str(sys.argv[1]))
    print('Building Primary HDU...')  
    #os.chdir('/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/RawData/') 
    objList,fitsList = numObjs()    
    ##TODO: put in logic to sort list of fits files if observer went back to the same source... 
    for objs in range(1,2):  
        fileList = fitsList[objs]
        allBanksList = [] ## master list of all BANKS for all FITS files associated with object
        numBanksList = [] ## number of BANKS associated with FITS file
        ## loop over FITS files for one object to construct a single SINGLE DISH binary FITS table
        for dataFITSFile in fileList:
            ## above file list does not contain fits files with BANK info
            ## get number of ints, int length, # chans, and proper bankList
            numInts, intLen, numChans, bankList = getScanInfo(fileList[0])
            ## append master BANK list
            allBanksList.extend(bankList)
            ## append number of BANKS
            numBanksList.append(len(bankList))
            ## process data per beam
            for beam in range(0,7): 
                ## structure of global buffer is:
                ## dim1: scan
                ## dim2: integrations
                ## dim3: bandpass
                globalDataBuff_X = np.zeros([int(len(fileList)), numInts, numChans * numBanks]) ##TODO: make general based on CovMode
                globalDataBuff_Y = np.zeros([int(len(fileList)), numInts, numChans * numBanks]) ##TODO: make general based on CovMode
                fileIdx = -1
                cnt = 0
                for fileName in bankList:
                    print(bankList)
                    print(fileList)        
                    print('\n')                
                    print('Beamforming correlations in: '+fileName[-25:]+', Beam: '+np.str(beam)) 
                    if cnt % len(bankList) == 0:
                         dataBuff_X = np.zeros([numInts, numChans * numBanks])   ##TODO: make mode independent  
                         dataBuff_Y = np.zeros([numInts, numChans * numBanks])
                         fileIdx+=1 
                    ## bank name is ALWAYS sixth-to-last character in string
                    bank = fileName[-6] 
                    ## grab xid from dictionary
                    xID = bankDict[bank]
                
                    ## Do the beamforming; returns processed BANK data 
                    ## (cov matrices to a beam-formed bandpass) in both
                    ## XX/YY Pols; returned data are in form: 
                    ## rows: ints, columns: freqChans
                    xData,yData = bf.getSpectralArray(fileName,beam, xID, bank)       
                    ## sort based on xid number for each integration
                    for ints in range(0,numInts):
                        bandpassSort(xID, dataBuff_X, xData, ints)
		        bandpassSort(xID, dataBuff_Y, yData, ints)
                    ## fill global data bufs
                    globalDataBuff_X[fileIdx,:,:] = dataBuff_X
                    globalDataBuff_Y[fileIdx,:,:] = dataBuff_Y
                    cnt+=1

                print('\n')
                ## build metadata; inputs are FITS file for ancillary files, numInts, global data buffers, int length            
                md = MetaDataModule(projectPath, fileList, allBanksList, numBanksList, globalDataBuff_X, globalDataBuff_Y, beam, projectPath, dataPath)
                thduList = md.constuctBinTableHeader()
                dataFITSFile[:-6] = dataFITSFile
                thduList.writeto(pwd+'/' + dataFITSFile + '_Beam'+str(beam)+'.fits')
    
       
    
#run main function
if __name__=="__main__":
    main()
