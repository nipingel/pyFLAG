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
import pickle
from modules.metaDataModule import MetaDataModule
from modules.beamformingModule import BeamformingModule
import matplotlib.pyplot as pyplot

## command line inputs
##TODO: exception handling
projectPath = sys.argv[1] ## of the form /home/gbtdata/AGBT16B_400_01
badScanList = sys.argv[2] ## list of bad timestamps included as a single string (e.g. 
                          ## '2017_08_01_01:01:01 2017_08_01_02:02:02'
badBankList = sys.argv[3] ## same format as above, but listing banks which stalled
## split up by space delimiter
badScanList = badScanList.split()
badBankList = badBankList.split()

## loop through and append '.fits' to each element for later ID
badScanList = [s + '.fits' for s in badScanList]
badBankList = [s + '.fits' for s in badBankList]

##global variables
numGPU = 10
numBanks = numGPU*2
##paths to ancillary FITS files
goFitsPath = projectPath + '/GO'

##TODO:Extend sorting from 4 FITS files to 20 FITS files.


## function to sort individual BANK data
## to the full bandpass based on XID. 
##TODO: add FRB mode functionality
def bandpassSort(xID, dataBuff, bankData):
    ## which correlation mode are we in?
    ## determine this based on number of channels 
    ## in data buff (second element in shape)
    ## number of integrations is the first element of dataBuff shape
    numInts = bankData.shape[0]
    numChans = bankData.shape[1]
    for ints in range(0, numInts):
        if numChans == 25:
            ## coarse channel mode
            ## position in bandpass dictated by xID
            bandpassStartChan = xID*5
            bandpassEndChan = bandpassStartChan+5
            ## get each chunk of five contigious channels from BANK data,
            ## and place it in the proper spot in full bandpass
            for i in range(0, 5):
                bankStartChan = i * 5
                bankEndChan = bankStartChan + 5
                dataBuff[ints, bandpassStartChan:bandpassEndChan] = bankData[ints, bankStartChan:bankEndChan]
            
                ## increment bandpassStartChan/bandpassEndChan by 100 for proper position in full bandpass
                bandpassStartChan += 100
                bandpassEndChan = bandpassStartChan+5
        elif numChans == 160:
            bandpassStartChan = xID*160
            bandpassEndChan = bandpassStartChan + 160
            ## ADDED TO AVOID ERROR THAT BANKA DATA LENGTH IS ZERO
            if len(bankData[ints,:]) == 0:
                dataBuff[ints, bandpassStartChan:bandpassEndChan] = np.zeros([160], dtype='float32')
            else:
                dataBuff[ints, bandpassStartChan:bandpassEndChan] = bankData[ints, :]
    return dataBuff


def weightSort(xID, weightBuff, weightData):
    numChans = weightData.shape[0]
    if numChans == 25:
        weightPassStart = xID*5
        weightPassEnd = weightPassStart+5
        for i in range(0,5):
            weightDataStart = i * 5
            weightDataEnd = weightDataStart + 5
            weightBuff[weightPassStart:weightPassEnd] = weightData[weightDataStart:weightDataEnd]        
            weightPassStart += 100
            weightPassEnd = weightPassStart + 5
    elif numChans == 160:
        weightPassStart = xID*160
        weightPassEnd = weightPassStart + 160
        weightBuff[weightPassStart:weightPassEnd] = weightData
    return weightBuff
## DEBUG       

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
def getScanInfo(fileName, dataPath):
    fitsLst = glob.glob(dataPath + fileName[:-5] + '*.fits')
    hdu = fits.open(fitsLst[0])
    intLen = np.float(hdu[0].header['REQSTI'])
    numInts = np.int(hdu[1].header['NAXIS2'])
    form = np.float(hdu[1].header['TFORM3'][:-1])
    data = hdu[1].data['DATA']
    numChans = data.shape[1]/2112
    #numChans = np.int(form/2112)
    return numInts, intLen, numChans, sorted(fitsLst)

def main():
    ## command line inputs
    ##TODO: exception handling
    projectPath = sys.argv[1] ## of the form /home/gbtdata/AGBT16B_400_01
    ## split project path to get project string
    projectPathSplit = projectPath.split('/')
    projectStr = projectPathSplit[-1]
    dataPath = '/lustre/projects/flag/' +  projectStr + '/BF/'
    bf = BeamformingModule(dataPath)

##TODO: check current directory permissions -- must have writing access
## the only commandline argument should be path to project/session ancillary FITS files    
    pwd = os.getcwd()
    pfb = False
    print('Project directory: ' + np.str(sys.argv[1]))
    print('Building Primary HDU...')
    objList,fitsList = numObjs()
    print(objList)
    ##TODO: put in logic to sort list of fits files if observer went back to the same source... 
    for objs in range(1,2):
        ## get source
        source = objList[objs]
        fileList = fitsList[objs]
        fileList = fileList[0:2]
        #fileList = fileList[0:41]
        ## remove bad scans from file list
        for s in badScanList:
            try:
                fileList.remove(s)
            except ValueError:
                pass
        print(fileList)
        ## loop over FITS files for one object to construct a single SINGLE DISH binary FITS table
        for beam in range(0,7):
            allBanksList = [] ## list of paths to all good BANK FITS files associated with object
            numBanksList = [] ## number of good BANKS associated with a scan
            ## above file list does not contain fits files with BANK info
            ## process data per beam
            ## open first FITS file to get relevant parameters
            
            ## global lists to hold data. Each element of the list contains a numpy array in the shape of (ints, freq chans)
            globalDataBuff_X_List = []
            globalDataBuff_Y_List = []
                
            fileCnt = 0
            for dataFITSFile in fileList:
                numInts, intLen, numSpecChans, bankList = getScanInfo(dataFITSFile, dataPath)
                ## initialize bank data buffers
                dataBuff_X = np.zeros([numInts, numSpecChans * numBanks])
                dataBuff_Y = np.zeros([numInts, numSpecChans * numBanks])
                ## remove bad BANKS
                for s in badBankList:
                    try:
                        bankList.remove(dataPath + s)
                    except ValueError:
                        pass 
                ## DEBUG
                xWeightBuff =  np.zeros([numSpecChans * numBanks], dtype = 'complex64')
                yWeightBuff =  np.zeros([numSpecChans * numBanks], dtype = 'complex64')
                ## DEBUG
                ## set bank counter
                bankCnt = 0
                for fileName in bankList:
                    print('\n')                
                    print('Beamforming correlations in: '+fileName[-25:]+', Beam: '+np.str(beam)) 
                    ## bank name is ALWAYS sixth-to-last character in string
                    bank = fileName[-6] 
                    ## DEBUG
                    if fileName[-25:] == '2017_07_28_06:22:39Q.fits':
                        continue
                    ## grab xid from dictionary
                    corrHDU = fits.open(fileName)
                    nRows = corrHDU[1].header['NAXIS2']
                    data = corrHDU[1].data['DATA']
                    ## DEBUG
                     
                    xID = np.int(corrHDU[0].header['XID'])

                    if nRows != 0:
                        ## Do the beamforming; returns processed BANK data 
                        ## (cov matrices to a beam-formed bandpass) in both
                        ## XX/YY Pols; returned data are in form: 
                        ## rows: ints, columns: freqChans
                        ## DEBUG
                        ##xData,yData = bf.getSpectralArray(fileName, data, beam, xID, bank)
                        xData,yData,xWeightBP, yWeightBP = bf.getSpectralArray(fileName, data, beam, xID, bank)       
                        ## sort based on xid number for each integration
                        dataBuff_X = bandpassSort(xID, dataBuff_X, xData)
                        dataBuff_Y = bandpassSort(xID, dataBuff_Y, yData)
                        ## append good BANK file to master list
                        allBanksList.append(fileName)
                        ## increment bank counter
                        bankCnt += 1  
                        ## fill global data bufs
                        #globalDataBuff_X[fileIdx,:,:] = dataBuff_X
                        #globalDataBuff_Y[fileIdx,:,:] = dataBuff_Y
                ## append to global buffer lists
                globalDataBuff_X_List.append(dataBuff_X)
                globalDataBuff_Y_List.append(dataBuff_Y)
                ## append bankCnt to list 
                numBanksList.append(bankCnt)
            print('\n')

            ## loop through global data buffer lists to find max ints for subsequent numpy array
            numChans = globalDataBuff_X_List[0].shape[1]
            if numChans == 3200:
                pfb = True
            ## save out important variables TEST
            #if beam == 0:
            #    with open('/users/npingel/FLAG/2017Reduction/globalDataBuffs_' + projectStr + '_Beam' + np.str(beam) + '.pickle', 'wb') as f:
            #        pickle.dump([globalDataBuff_X_List, globalDataBuff_Y_List],  f)
            
            ## build metadata; inputs are FITS file for ancillary files, numInts, global data buffers, int length            
            md = MetaDataModule(projectPath, dataPath, fileList, allBanksList, numBanksList, globalDataBuff_X_List, globalDataBuff_Y_List, beam, pfb)
            thduList = md.constuctBinTableHeader()
            dataFITSFile = dataFITSFile[:-6]
            thduList.writeto(pwd+'/' + projectStr + '_' + source + '_Beam'+str(beam) + '.fits')
    
       
    
#run main function
if __name__=="__main__":
    main()
