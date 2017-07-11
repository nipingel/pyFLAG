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
    ## in bandpasss (second element in shape)
    ## Number of integrations are the first element
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
            dataBuff[ints, bandpassStartChan:bandpassEndChan] = bankData[ints, :]
    return dataBuff

## DEBUG
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
    fitsLst = glob.glob(dataPath + fileName[:-6] + '*.fits')
    hdu = fits.open(fitsLst[0])
    intLen = np.float(hdu[0].header['REQSTI'])
    numInts = np.int(hdu[1].header['NAXIS2'])
    form = np.float(hdu[1].header['TFORM3'][:-1])
    data = hdu[1].data['DATA']
    numChans = data.shape[1]/2112
    #numChans = np.int(form/2112)
    return numInts, intLen, numChans, np.sort(fitsLst)

def main():
    ## command line inputs
    ##TODO: exception handling
    projectPath = sys.argv[1] ## of the form /home/gbtdata/AGBT16B_400_01
    ## split project path to get project string
    projectPathSplit = projectPath.split('/')
    projectStr = projectPathSplit[-1]
    dataPath = '/lustre/projects/flag/' +  projectStr + '/BF/'

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
    pfb = False
    print('Project directory: ' + np.str(sys.argv[1]))
    print('Building Primary HDU...')  
    #os.chdir('/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/RawData/') 
    objList,fitsList = numObjs()
    print(objList)
    ##TODO: put in logic to sort list of fits files if observer went back to the same source... 
    for objs in range(0,1):
        fileList = fitsList[objs]
        fileList = fileList[10:12]
        print(fileList)
        allBanksList = [] ## master list of all BANKS for all FITS files associated with object
        numBanksList = [] ## number of BANKS associated with FITS file
        ## loop over FITS files for one object to construct a single SINGLE DISH binary FITS table
        for beam in range(0,7):
            ## above file list does not contain fits files with BANK info
            ## process data per beam
            ## open first FITS file to get relevant parameters
            numInts, intLen, numChans, bankList = getScanInfo(fileList[0], dataPath)
            ## structure of global buffer is:
            ## dim1: scan
            ## dim2: integrations
            ## dim3: bandpass
            ## TODO: WHAT HAPPENS WHEN WE HAVE DIFFERENT NUMBER OF INTEGRATIONS!!
            globalDataBuff_X = np.zeros([int(len(fileList)), numInts, numChans * numBanks])
            globalDataBuff_Y = np.zeros([int(len(fileList)), numInts, numChans * numBanks])
            fileIdx = 0
            for dataFITSFile in fileList: 
                numInts, intLen, numSpecChans, bankList = getScanInfo(dataFITSFile, dataPath)
                ## append master BANK list
                allBanksList.extend(bankList)
                ## append number of BANKS
                numBanksList.append(len(bankList))
                if fileIdx == 0:
                    numBanksList[0] = numBanksList[0] - 1          
                ## initialize bank data buffers
                dataBuff_X = np.zeros([numInts, numSpecChans * numBanks])
                dataBuff_Y = np.zeros([numInts, numSpecChans * numBanks])
                
                ## DEBUG
                xWeightBuff =  np.zeros([numSpecChans * numBanks], dtype = 'complex64')
                yWeightBuff =  np.zeros([numSpecChans * numBanks], dtype = 'complex64')
                ## DEBUG
                
                
                for fileName in bankList:
                    print('\n')                
                    print('Beamforming correlations in: '+fileName[-25:]+', Beam: '+np.str(beam)) 
                    ## bank name is ALWAYS sixth-to-last character in string
                    bank = fileName[-6] 
                    ## grab xid from dictionary
                    corrHDU = fits.open(fileName)
                    nRows = corrHDU[1].header['NAXIS2']
                    data = corrHDU[1].data['DATA']
                    xID = np.int(corrHDU[0].header['XID'])
                    if nRows != 0:
                        ## Do the beamforming; returns processed BANK data 
                        ## (cov matrices to a beam-formed bandpass) in both
                        ## XX/YY Pols; returned data are in form: 
                        ## rows: ints, columns: freqChans
                        ## DEBUG
                        ##xData,yData = bf.getSpectralArray(fileName, data, beam, xID, bank)
                        xData,yData,xWeightBP, yWeightBP = bf.getSpectralArray(fileName, data, beam, xID, bank)       
                        ## DEBUG
                          
                        ## DEBUG
                        ## plot weight bandpasses
                        #pyplot.figure()
                        #pyplot.title('Beamformed Bank Bandpass; XID: ' + np.str(xID))
                        #pyplot.title('Weight Bandpass')
                        #pyplot.plot(np.abs(xWeightBP), label = 'XX-Pol')
                        #pyplot.plot(np.mean(yData, axis=0), label = 'YY-Pol')
                        #pyplot.xlabel('Frequency [MHz]')
                        #pyplot.ylabel('Magnitude')
                        #pyplot.legend(loc=0)
                        #pyplot.show()
                        #pyplot.savefig('/users/npingel/FLAG/2017Reduction/PFBTests/Plots/BankBandPass_' + 'XID_' + np.str(xID) + '_' + np.str(objList[objs]) + '.pdf')
                        #pyplot.savefig('/users/npingel/FLAG/2017Reduction/PFBTests/Plots/WeightBandpass_M101_Fine.pdf')
                        ## DEBUG
                        
                        ## sort based on xid number for each integration
                        dataBuff_X = bandpassSort(xID, dataBuff_X, xData)
    		        dataBuff_Y = bandpassSort(xID, dataBuff_Y, yData)
                        
                        ## DEBUG 
                        xWeightBuff = weightSort(xID,xWeightBuff, xWeightBP)
                        yWeightBuff = weightSort(xID,yWeightBuff, yWeightBP)
                        ## DEBUG
                         
                        ## fill global data bufs
                        globalDataBuff_X[fileIdx,:,:] = dataBuff_X
                        globalDataBuff_Y[fileIdx,:,:] = dataBuff_Y
                ## increment fileIdx for global data buffers
                fileIdx += 1 
            print('\n')
            if globalDataBuff_X.shape[2] == 3200:
                pfb = True
            
            ## DEBUG 
            ## magnitude   
            pyplot.figure()
            pyplot.title('Weight Magnitude Bandpass (' + projectStr + ')')
            pyplot.plot(np.abs(xWeightBuff)**2, label = 'XX-Pol')
            pyplot.plot(np.abs(yWeightBuff)**2, label = 'YY-Pol')
            pyplot.xlabel('Frequency Chanel')
            pyplot.ylabel('Magnitude')
            pyplot.legend(loc=0)
            pyplot.savefig('/users/npingel/FLAG/2017Reduction/Plots/Weights/' + projectStr + '_weightBandPass_Beam' + np.str(beam) + '_mag.pdf')
            pyplot.show()
            
            ## amplitdue
            pyplot.figure()
            pyplot.title('Weight Amplitude Bandpass (' + projectStr + ')')
            pyplot.plot(np.abs(xWeightBuff), label = 'XX-Pol')
            pyplot.plot(np.abs(yWeightBuff), label = 'YY-Pol')
            pyplot.xlabel('Frequency Chanel')
            pyplot.ylabel('Amplitude')
            pyplot.legend(loc=0)
            pyplot.savefig('/users/npingel/FLAG/2017Reduction/Plots/Weights/' + projectStr + '_weightBandPass_Beam' + np.str(beam) + '_amp.pdf')
            pyplot.show()
            ## phase
            pyplot.figure()
            pyplot.title('Weight Phase Bandpass (' + projectStr + ')')
            pyplot.plot(np.angle(xWeightBuff), label = 'XX-Pol')
            pyplot.plot(np.angle(yWeightBuff), label = 'YY-Pol')
            pyplot.xlabel('Frequency Chanel')
            pyplot.ylabel('Phase [rad]')
            pyplot.legend(loc=0)
            pyplot.savefig('/users/npingel/FLAG/2017Reduction/Plots/Weights/' + projectStr + '_weightBandPass_Beam' + np.str(beam) + '_phase.pdf')
            pyplot.show()

            ## save out important variables TEST
            with open('/users/npingel/FLAG/2017Reduction/WeightBandpass_' + projectStr + '_Beam' + np.str(beam) + '.pickle', 'wb') as f:
                pickle.dump([xWeightBuff, yWeightBuff],  f)
            ## DEBUG
            ## build metadata; inputs are FITS file for ancillary files, numInts, global data buffers, int length            
            md = MetaDataModule(projectPath, dataPath, fileList, allBanksList, numBanksList, globalDataBuff_X, globalDataBuff_Y, beam, pfb)
            thduList = md.constuctBinTableHeader()
            dataFITSFile = dataFITSFile[:-6]
            thduList.writeto(pwd+'/' + projectStr + '_Beam'+str(beam)+'.fits')
    
       
    
#run main function
if __name__=="__main__":
    main()
