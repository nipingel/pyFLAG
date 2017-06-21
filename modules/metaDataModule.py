# -*- coding: utf-8 -*-
"""
This module (i.e. class) will collate the required metadata to construc the primary HDU and SDFITS binary table. 

@author: npingel
"""
from astropy.io import fits
from astropy.time import Time
import os 
import numpy as np
import datetime
import collections
import sys

##globalPaths
rawDataPath = 
goDataPath = '/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/GO'
antDataPath = '/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/Antenna'
lo1APath = '/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/LO1A/'
weightDataPath = '/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/Weights/'
fitsPath = '/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/'
keywordPath = '/users/npingel/FLAG/SpectralFiller/misc/sdKeywords.txt'

class MetaDataModule:
    
        
    ## initialize function. Opens raw data file, numInts, data buffers (banpasses), beam, freqChans, intLength
    ## total BANK files, and creates column object. 
    def __init__(self, projectPath, fitsList, bankFitsList, numBanksList, dataBuff_X,dataBuff_Y,beamNum, ancPath, rawDataPath):        
        self.fitsList = projectPath +      
        self.bankFitsList = bankFitsList
        self.numBanksList = numBanksList
        os.chdir(rawDataPath) ##TODO: run on flag3
        hdu = fits.open(bankFitsList[0])
        dmjd = hdu[1].data['DMJD']
        self.dmjd = dmjd
        self.dataBuff_X = dataBuff_X
        self.dataBuff_Y = dataBuff_Y
        self.beamNum = beamNum 
        self.ancPath = ancPath
        self.dataPath = rawDataPath
        self.valueArr = None
        self.newArr = None
        self.numFreqChans = len(dataBuff_X[0,0,:])
        self.Column = collections.namedtuple('Column',['param','valueArr','comment'])
        self.cols = None
        self.numPhases = 2

        ## function to initialize self.valueArr/self.newArr based on datatype and fileNum iteration
        def arrayInit(self, fileNum, numInts, dtype, valStr):
            ## if first iteration, initialize the global array to hold all ints info; 
            ## only scan phases are the two pols so self.numPhases = 2. Makes a character array if
            ## data types are not numeric
            if fileNum == 0: 
                if not any(dtype == 'float32', dtype == 'int32', dtype == 'int16']):
                    self.valueArr = np.chararray([self.numPhases * NumInts], itemsize=len(valStr))
                else   
                    self.valueArr = np.empty([self.numPhases * numInts], dtype=dtype)
                if fill == True:
                    self.valueArr
            ## update self.newArr if we are past the first iteration of the loop over scan FITS files
            else:
                if not any(dtype == 'float32', dtype == 'int32', dtype == 'int16']):
                    self.newArr = np.chararray([self.numPhases * NumInts], itemsize=len(valStr))
                else
                    self.newArr = np.empty([self.numPhases * numInts], dtype=dtype)
                    
                      

        def getSMKey(self,param,singleVal = False):
            ## iterate flag defaults to false. Is set True based on parameter
            iterate = False
            for fileNum in range(0,len(self.fitsList)):
                corrHDU = fits.open(self.bankFitsList[bandIdx])
                scanNumInts = corrHdu[1].header['NAXIS2']
                 ## determine time spent in state (e.g. calOn, calOff); 
                 ## since we do no freq sw or cal sw, this is equal to scan time
                 if param == 'DURATION':
                     self.initArr(fileNum, numScanInts, 'float32', None)
                     value = corrHDU[0].header['SCANLEN']
                 ## retrieve actual integration time
                 elif param == 'EXPOSURE':
                     paramLook = 'ACTSTI'
                     self.initArr(fileNum, numScanInts, 'float32', None)
                     value = corrHDU[0].header['ACTSTI']
                 ## retrieve correlation mode
                 elif param == 'OBSMODE':
                     self.initArr(fileNum, numScanInts, 'float32', None)
                     value = corrHDU[0].header['COVMODE']
                 ## function that is called in getModeDepParams; returns a list of correlation modes for FITS files associated with object
                 elif any(param == 'BANDWID', param == 'CRPIX1', param == 'CDELT1', param == 'FREQRES']) and singleVal == True:
                     paramLook = 'COVMODE'
                     covModeList = []
                     ## loop through first BANK files to compile list of COVMODES
                     for bankIdx in range(0, len(self.dataFitsList), self.numBanks):
                         corrHDU = fits.open(dataFitsList[bankIdx])
                         covModeList.append(corrHDU[0].header[paramLook])
                     return covModeList 
                 ## determine start time of each integration
                 elif param == 'DATE-OBS':
                     scanDMJD = corrHDU[1].data['DMJD']
                     ## cast DMJD to correct string format
                     timeObj = Time(scanDMJD, format = 'mjd', scale='utc') 
                     ## initialize array to store parameter values
                     self.initArr(fileNum, numScanInts, 'str', timeObj)
                     TimeVal == True
                 bankIdx += self.numBanksList[fileNum]
            
                 ## fill in value array for similar values across integrations/polarizations
                 if timeVal == True and fileNum == 0:
                     self.valueArr[0::2] = timeObj.isot
                     self.valueArr[1::2] = timeObj.isot
                 elif timeVal == True and FileNum > 0:
                     self.newArr[0::2] = timeObj.isot
                     self.newArr[1::2] = timeObj.isot
                     self.valueArr = np.concatenate([self.valueArr, self.newArr])
                 elif timeVal == False and FileNum == 0:
                     self.valueArr[:] = value
                 else:
                     self.newArr[:] = value
                     ## concatenate new array to global value array
                     self.valueArr = np.concatenate([self.valueArr, self.newArr])
                 ## update BANK index
                 bankIdx += self.numBanksList[fileNum]
            ## retrieve comment and construct column    
            comment = self.commentDict[param]            
            self.Column.param = param            
            self.Column.comment = comment
            self.Column.valueArr = self.valueArr
        
        ## function to get values for GO FITS file parameters    
        ##TODO: error handling (e.g. keyword not found)
        def getGOFITSParam(self,param):
            idx = 0          
            if param == 'TRGTLONG':
                paramLook = 'RA'
            elif param == 'TRGTLAT':
                paramLook = 'DEC'           
            elif param == 'CTYPE2' or param == 'CTYPE3':
                paramLook = 'COORDSYS'
            else:
                paramLook = param
            bankIdx = 0
            for fileNum in range(0,len(self.fitsList)):            
                goHDU = fits.open(self.projectPath + '/' + self.fitsList[fileNum])                
                corrHDU = fits.open(self.bankFitsList[bandIdx])
                scanNumInts = corrHDU[1].header['NAXIS2']
                if paramLook != 'TIMESTAMP':           
                    value = goHDU[0].header[paramLook]
                else:
                    value = None
                if value == 'GALACTIC' and param == 'CTYPE2':
                    valStr = 'GLON'
                    ## initialize array to hold values if first iteration
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif value == 'GALACTIC' and param == 'CTYPE3':
                    valStr = 'GLAT'
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif value == 'RADEC' and param == 'CTYPE2':
                    valStr = 'RA'
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif value == 'RADEC' and param == 'CTYPE3':
                    valStr = 'DEC'
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif value == 'HADEC' and param == 'CTYPE2':
                    valStr = 'HA'
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif value == 'HADEC' and param == 'CTYPE3':
                    valStr = 'DEC'
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif value == 'AZEL' and param == 'CTYPE2':
                    valStr = 'AZ'
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif value == 'AZEL' and param == 'CTYPE3':
                    valStr = 'EL'
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif param =='TRGTLONG':
                    valStr = value
                    self.initArr(fileNum, numScanInts, 'float32', None)
                elif param =='TRGTLAT':
                    valStr = value
                    self.initArr(fileNum, numScanInts, 'float32', None)
                elif param  == 'TIMESTAMP':
                    valStr = self.fitsList[fileNum]
                    valStr = valStr[:-6]
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif param =='VELOCITY':
                    valStr = value
                    self.initArr(fileNum, numScanInts, 'float32', None)
                elif param =='OBJECT':
                    valStr = value
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif param == 'OBSERVER':
                    valStr = value
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif param == 'OBSID':
                    valStr = 'unknown'
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif param == 'RESTFRQ':
                    valStr = 1450.00*1e6##TODO:remove for production
                    self.initArr(fileNum, numScanInts, 'float32', None)
                elif param == 'SCAN' or param == 'PROCSEQN' or param == 'PROCSIZE':
                    valStr = value
                    self.initArr(fileNum, numScanInts, 'int32', None)
                elif param == 'PROCSCAN' or param == 'PROCTYPE':
                    valStr = value
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                elif param =='EQUINOX':
                    valStr = value
                    self.initArr(fileNum, numScanInts, 'float32', None)
                elif param == 'LASTON':
                    valStr = value
                    self.initArr(fileNum, numScanInts, 'int16', None)
                elif param == 'LASTOFF':
                    valStr = value
                    self.initArr(fileNum, numScanInts, 'int16', None)
                elif param == 'RADESYS':
                    coordSys = goHDU[0].header['COORDSYS']
                    if coordSys == 'RADEC':
                        valStr = value
                    else:
                        valStr=''
                    self.initArr(fileNum, numScanInts, 'str', valStr)
                ## fill initial array
                if fileNum == 0:
                    self.valueArr[:] = valStr
                ## concatenate subsequent arrays
                else:
                    self.newArr[:] = valStr
                    self.valueArr = np.concatenate([valueArr, newArr])  
                ## update bankIdx
                bankIdx += self.numBanksList[fileNum]  
            ## retrieve comment and construct column 
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = self.valueArr
            self.Column.comment = comment

        ## function to get values for GO FITS file parameters
        def getAntFITSParam(self,param):                      
            iterate = True
            bankIdx = 0
            for fileNum in range(0,len(self.fitsList)):            
                antHDU = fits.open(self.fitsList[fileNum])
                antDMJD = antHDU[2].data['DMJD']
                ## get integrations this scan
                corrHDU = fits.open(self.bankFitsList[bankIdx])
                scanNumInts = corrHDU[1].header['NAXIS2']
                if param == 'CRVAL2':
                    ## interpolate data time samples to data time samples
                    maj = antHDU[2].data['MAJOR']
                    value = np.interp(self.dmjd,antDMJD,maj)
                    ## initialize array to hold parameter values if first iteration
                    self.initArr(fileNum, numScanInts, 'float32', None)
                elif param == 'CRVAL3':
                    minor = antHDU[2].data['MINOR']
                    value = np.interp(self.dmjd,antDMJD,minor)
                    self.initArr(fileNum, numScanInts, 'float32', None)
                elif param == 'AZIMUTH':
                    az = antHDU[2].data['MNT_AZ']
                    value = np.interp(self.dmjd,antDMJD,az)
                    self.initArr(fileNum, numScanInts, 'float32', None)
                elif param == 'ELEVATIO':
                    el = antHDU[2].data['MNT_EL']
                    value = np.interp(self.dmjd,antDMJD,el)
                    self.initArr(fileNum, numScanInts, 'float32', None)
                elif param == 'LST':
                    param1 = 'LSTSTART'
                    lstStart = antHDU[0].header[param1]
                    param2 = 'ACTSTI'
                    intLen = np.float(corHDU[0].header[param2])
                    self.initArr(fileNum, numScanInts, 'float32', None)
                    idx = 0
                    value = np.zeros(scanNumInts, dtype='float32')
                    for i in range(0,scanNumInts):
                        value[i] = lstStart + intLen/2+(i*intLen)
                elif param == 'TAMBIENT':
                    value = antHDU[0].header['AMBTEMP']
                    value+=273.0
                    if fileNum == 0:
                        valueArr = np.full([2 * scanNuInts], dtype='float32', fill_value = value)
                    else:
                        newArr = np.full([2 * scanNuInts], dtype='float32', fill_value = value)
                elif param == 'HUMIDITY':
                    value = antHDU[0].header['AMBHUMID']
                    if fileNum == 0:
                        valueArr = np.full([2 * scanNuInts], dtype='float32', fill_value = value)
                    else:
                        newArr = np.full([2 * scanNuInts], dtype='float32', fill_value = value)
                elif param == 'PRESSURE':
                    value = antHDU[0].header['AMBPRESS']
                    value=value*0.75006375541921    
                    if fileNum == 0:
                        valueArr = np.full([2 * scanNuInts], dtype='float32', fill_value = value)
                    else:
                        newArr = np.full([2 * scanNuInts], dtype='float32', fill_value = value)
                if fileNum == 0:
                    self.valueArr[0::2] = value
                    self.valueArr[0::1] = value
                else:    
                    self.newArr[0::2] = value
                    self.newArr[1::2] = value
                    self.valueArr = np.concatenate(self.valueArr, self.newArr
                bankIdx += self.numBanksList[fileNum]
            ## retrieve comment and create column 
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
        
            ## delete variables: bankIdx, valueArr, newArr
            del, bankIdx, valueArr, newArr
    
        def getRcvrFITSParam(self,param):
            if param == 'FRONTEND':            
                valStr = 'PAF' ##TODO: search for value in PAF manager for production
            valueArr = np.chararray([2*self.numInts*len(self.fitsList)],itemsize=len(value))    
            bankIdx = 0
            for fileNum in range(0,len(self.fitsList)):            
                corrHDU = fits.open(self.bankFitsList[bankIdx])
                scanNumInts = corrHDU[1].header['NAXIS2']
                self.initArr(fileNum, numScanInts, 'str', valStr)
                if fileNum == 0:
                    self.valueArr[:] = valStr
                else: 
                    self.newArr[:] = valStr
                    self.valueArr = np.concatenate(self.valueArr, self.newArr
                 
                bankIdx += self.numBanksList[fileNum]
                    
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment    
       
        def getArbParam(self,param):           
            if param == 'TSYS' or param == 'TCAL' or param == 'DOPFREQ':    
                valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='float32')                
                valueArr.fill(1.0)
            elif param == 'BEAM':    
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='float32')                
                valueArr.fill(self.beamNum)
            elif param == 'FEED' or param == 'SUBREF_STATE' or param == 'QD_BAD':    
                valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='int16')                
                valueArr.fill(1)
            elif param == 'CRVAL4':
                valueArr = np.zeros([2*self.numInts*len(self.fitsList)], dtype='int16')
                valueArr[0::2] = -5
                valueArr[1::2] = -6
            elif param == 'ZEROCHAN' or param =='TWARM' or param =='TCOLD' or param == 'QD_XEL' or param == 'QD_EL' or param == 'CALPOSITION':
                valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='float32')                
                valueArr.fill(np.nan)
            elif param == 'FDNUM' or param == 'IFNUM':
                valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='int16')                 
                valueArr.fill(0.0)
            elif param == 'SIG' or param == 'CAL' or param == 'CALTYPE':
                valStr = param
                valueArr = np.chararray([2*self.numInts*len(self.fitsList)],itemsize=1)  
                valueArr.fill('T')
            elif param == 'QD_METHOD':
                valueArr = np.chararray([2*self.numInts*len(self.fitsList)],itemsize=1)              
                valueArr.fill('C')
            elif param == 'TUNIT7':
                valStr = 'counts'
                valueArr = np.chararray([2*self.numInts*len(self.fitsList)],itemsize=len(valStr))                
                valueArr[:] = valStr 
            elif param == 'CTYPE1':
                valStr = 'FREQ-OBS'
                valueArr = np.chararray([2*self.numInts*len(self.fitsList)],itemsize=len(valStr))                
                valueArr[:] = valStr                 
            elif param == 'SRFEED':
                valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='int16')
                valueArr.fill(0)
            elif param == 'SIDEBAND':
                valueArr = np.chararray([2*self.numInts*len(self.fitsList)],itemsize=1)
                valueArr.fill('L')
            elif param =='SAMPLER':
                valStr1 = 'A1_0'    
                valStr2 = 'A2_0'
                valueArr = np.chararray([2*self.numInts*len(self.fitsList)],itemsize=len(valStr1))  
                valueArr[0::2] = valStr1
                valueArr[1::2] = valStr2
            elif param == 'DOPFREQ':
                valStr = 1432.729*1e6 ##TODO: fix for production. 
                valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='float32')
                valueArr[:] = valStr
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
        
        def getDataParam(self,param):
            if param == 'DATA':
                valueArr = np.empty([2*self.numInts*len(self.fitsList),self.numFreqChans],dtype='float32')
                valueArr[::2,:] = self.dataBuff_Y[0,:,:]
                valueArr[1::2,:] = self.dataBuff_X[0,:,:]            
            elif param == 'TDIM7':
                valStr = '['+np.str(self.numFreqChans)+',1,1,1]'
                valueArr = np.chararray([2*self.numInts*len(self.fitsList)],itemsize=len(valStr))                
                valueArr[:] = valStr             
            elif param == 'PLNUM':
                valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='int16')
                valueArr[::2] = 0
                valueArr[1::2] = 1
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
        
        def getLOFITSParam(self,param, repeat = True ): 
            os.chdir(lo1APath) ##TODO: run on flag3/GB
            if param == 'VELDEF':
                loHDU = fits.open(self.fitsList[0])   
                tblNum = 2
                value = loHDU[tblNum].header[param]
                valueArr = np.chararray([2*self.numInts*len(self.fitsList)],itemsize=len(value)) 
                valueArr[:] = value
                repeat = False
            elif param == 'VFRAME' or param == 'RVSYS':
                valueArr = np.zeros([2*self.numInts*len(self.fitsList)],dtype='float32')            
                tblNum = 3             
                repeat = True
                paramLook = param
            elif param == 'CRVAL1' or 'OBSFREQ':
                valueArr = np.zeros([2*self.numInts*len(self.fitsList)],dtype='float32')
                paramLook = 'LO1FREQ'
                tblNum = 3
            idx=0            
            for file in range(0,len(self.fitsList)):            
                loHDU = fits.open(self.fitsList[file])       
                if repeat == True:        
                    value = loHDU[tblNum].data[paramLook]
                    valueArr[(idx*2*self.numInts):(idx*2*self.numInts)+2*self.numInts:2] = value
                    valueArr[(idx*2*self.numInts)+1:(idx*2*self.numInts)+2*self.numInts:2] = value
                idx+=1
               
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
        
        def getBeamOffsets(self,param):
            os.chdir(weightDataPath)            
            hdu = fits.open('2016_07_25_04:32:35_xid0_weights.fits') ##TODO: work on flag3
            if param == 'FEEDXOFF':
                beamOff_Az = hdu[1].data['BeamOff_AZ']
                value = beamOff_Az[self.beamNum]
            else:
                beamOff_El = hdu[1].data['BeamOff_EL']
                value = beamOff_El[self.beamNum]
            valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='float32')
            valueArr.fill(value) 
            
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
                
        
        def getModeDepParams(self,param):
            if param == 'CDELT1' or param == 'FREQRES':
                modeName = getSMKey(self,param,singleVal = True)
                if modeName == 'PAF_CAL':
                    value = 303.18*1000.
                    valueArr = np.zeros([2*self.numInts*len(self.fitsList)],dtype='float32')    
            elif param == 'CRPIX1':
                modeName = getSMKey(self,param,singleVal = True)
                if modeName == 'PAF_CAL':
                    value = 500/2.
                    valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='float32')
            elif param == 'BANDWID':
                modeName = getSMKey(self,param,singleVal = True)
                if modeName == 'PAF_CAL':
                    value = 303.18*500*1000.
                    valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='float32')
            idx = 0
            for file in range(0,len(self.fitsList)):
                valueArr[(idx*2*self.numInts):(idx*2*self.numInts)+2*self.numInts:2] = value
                valueArr[(idx*2*self.numInts)+1:(idx*2*self.numInts)+2*self.numInts:2] = value
                idx+=1
                    
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
        
                
        self.commentDict = {'OBJECT':'name of source observed',
                            'OBSERVER':'name of observer(s)',                            
                            'BANDWID':'bandwidth',
                            'DATE-OBS':'date and time of observation start',
                            'DURATION':'total integration duration in seconds',
                            'EXPOSURE':'effective int time (excludes blanking) in secs',
                            'TSYS':'system temperature in Kelvin',                            
                            'OBSID':'observation description',
                            'SCAN':'scan number',
                            'RESTFRQ':'rest frequency at band center',
                            'EQUINOX':'equinox of selected coordinate reference frame',
                            'RADESYS':'Equitorial coordinate system name',
                            'TRGTLONG':'target longitude in coord. ref. frame',
                            'TRGTLAT':'target latitude in coord. ref. frame',
                            'PROCSEQN':'scan sequence number',
                            'PROCSIZE':'number of scans in procedure',
                            'PROCSCAN':'scan\'s role in the procedure',
                            'PROCTYPE':'type of the procedure',
                            'LASTON':'last \'on\' for position switching',
                            'LASTOFF':'last \'off\' for position switching',
                            'VELOCITY':'line velocity in rest frame',
                            'FRONTEND':'frontend device',
                            'TAMBIENT':'ambient temperature',
                            'PRESSURE':'ambient pressure',
                            'HUMIDITY':'relative humidity',
                            'TIMESTAMP':'date and time of scan start',
                            'AZIMUTH':'azimuth',
                            'ELEVATIO':'elevation',
                            'ZEROCHAN':'zero channel',
                            'SIG':'signal is true; reference is false',
                            'TWARM':'4mm RX ambient load temp (K)',
                            'TCOLD':'4mm RX ambient cold temp (K)',
                            'FDNUM':'feed number',
                            'IFNUM':'Spectral window (IF) number',
                            'TUNIT7':'',
                            'TDIM7':'data dimensions of the array',
                            'CTYPE1':'first data axis is frequency-like',
                            'DATA':'actual data',
                            'PLNUM':'Polarization number',
                            'OBSMODE':'observing mode',
                            'TCAL':'always 1.0',
                            'VELDEF':'velocity definition and frame',
                            'VFRAME':'radial velocity of the reference frame',
                            'RVSYS':'radial velocity, Vsource - Vtelescope',
                            'CAL':'No noise diode; always \'T\'',
                            'CALTYPE':'No noise diode; always \'T\'',
                            'CALPOSITION':'No noise diode; always \'T\'',
                            'FEED':'(signal) feed number; always 1 for PAF',
                            'SAMPLER':'Sampler ID; always A1_0 for FLAG',
                            'SRFEED':'reference feed number',
                            'FEEDXOFF':'Beam Offset in crossEl from boresight [arcmin]',
                            'FEEDEOFF':'Beam Offset in El from boresight [arcmin]',
                            'SUBREF_STATE':'Subreflector position when nodding',
                            'SIDEBAND':'Resulting sideband; Always L for PAF',
                            'QD_XEL':'QuadrantDetector cross-elevation offset ',
                            'QD_EL':'QuadrantDetector elevation offset',
                            'QD_BAD':'QuadrantDetector flag: 0=good,1=bad',
                            'QD_METHOD':'Quad. Det. method A,B,C. Blank indicates none.',
                            'CRVAL1':'',            
                            'OBSFREQ':'observed center frequency', 
                            'DOPFREQ':'Doppler tracked frequency', 
                            'FREQRES':'frequency resolution', 
                            'LST': 'local sidereal time (LST)', 
                            'CTYPE2':'second axis is longitude-like',
                            'CTYPE3': 'third axis is latitude-like',
                            'CRVAL1':'',
                            'CRPIX1':'',
                            'CDELT1':'',
                            'CRVAL2':'',
                            'CRVAL3':'',
                            'CRVAL4':'',
                            'TRGTLONG':'target longitude in coord. ref. frame',
                            'TRGTLAT':'target latitude in coord. ref. frame'
                            
              }
        self.funcDict = {'OBJECT':getGOFITSParam,
                         'CTYPE1':getArbParam,
                         'CRVAL1':getLOFITSParam,
                         'DATE-OBS':getSMKey,
                         'DURATION':getSMKey,
                         'EXPOSURE':getSMKey,
                         'TSYS':getArbParam,
                         'OBSERVER':getGOFITSParam,
                         'OBSID':getGOFITSParam,
                         'SCAN':getGOFITSParam,
                         'RESTFRQ':getGOFITSParam,
                         'EQUINOX':getGOFITSParam,
                         'RADESYS':getGOFITSParam,
                         'TRGTLONG':getGOFITSParam,
                         'TRGTLAT':getGOFITSParam,
                         'PROCSEQN':getGOFITSParam,
                         'PROCSIZE':getGOFITSParam,
                         'PROCSCAN':getGOFITSParam,
                         'PROCTYPE':getGOFITSParam,
                         'LASTON':getGOFITSParam,
                         'LASTOFF':getGOFITSParam,
                         'VELOCITY':getGOFITSParam,
                         'FRONTEND':getRcvrFITSParam,
                         'TAMBIENT':getAntFITSParam,
                         'PRESSURE':getAntFITSParam,
                         'HUMIDITY':getAntFITSParam,
                         'TIMESTAMP':getGOFITSParam,
                         'AZIMUTH':getAntFITSParam,
                         'ELEVATIO':getAntFITSParam,
                         'ZEROCHAN':getArbParam,
                         'SIG':getArbParam,
                         'TWARM':getArbParam,
                         'TCOLD':getArbParam,
                         'FDNUM':getArbParam,
                         'IFNUM':getArbParam,
                         'TUNIT7':getArbParam,
                         'TDIM7':getDataParam,
                         'DATA':getDataParam,
                         'PLNUM':getDataParam,
                         'OBSMODE':getSMKey,
                         'TCAL':getArbParam,
                         'VELDEF':getLOFITSParam,
                         'VFRAME':getLOFITSParam,
                         'RVSYS':getLOFITSParam,
                         'OBSFREQ':getLOFITSParam,
                         'VRSYS':getLOFITSParam,
                         'CAL':getArbParam,
                         'CALTYPE':getArbParam,
                         'CALPOSITION':getArbParam,
                         'FEED':getArbParam,
                         'SAMPLER':getArbParam,
                         'SRFEED':getArbParam,
                         'FEEDXOFF':getBeamOffsets,
                         'FEEDEOFF':getBeamOffsets,
                         'SUBREF_STATE':getArbParam,
                         'SIDEBAND':getArbParam, 
                         'QD_XEL':getArbParam,
                         'QD_EL':getArbParam,
                         'QD_BAD':getArbParam,
                         'QD_METHOD':getArbParam,
                         'DOPFREQ':getArbParam, 
                         'FREQRES':getModeDepParams, 
                         'LST': getAntFITSParam,
                         'CTYPE2':getGOFITSParam,
                         'CTYPE3':getGOFITSParam,
                         'CRPIX1':getModeDepParams,
                         'CDELT1':getModeDepParams,
                         'BANDWID':getModeDepParams,
                         'CRVAL2':getAntFITSParam,
                         'CRVAL3':getAntFITSParam,
                         'CRVAL4':getArbParam, 
                         'TRGTLONG':getGOFITSParam,
                         'TRGTLAT':getGOFITSParam
                         }
        ##Parameter Dictionary
        self.keyToParamDict = {'XTENSION':'BINTABLE',
              'TTYPE1':'OBJECT',
              'TFORM1':'32A',
              'TUNIT1':'',
              'TTYPE2': 'BANDWID', ##TODO: finish with other modes
              'TFORM2':'1D',
              'TUNIT2': 'Hz',
              'TTYPE3': 'DATE-OBS',
              'TFORM3': '22A',
              'TUNIT3': '',
              'TTYPE4': 'DURATION',
              'TFORM4':'1D',
              'TUNIT4': 's',
              'TTYPE5': 'EXPOSURE',
              'TFORM5': '1D',
              'TUNIT5': 's', 
              'TTYPE6': 'TSYS', 
              'TFORM6': '1D',
              'TUNIT6': 'K',
              'TTYPE7': 'DATA', 
              'TFORM7':'',
              'TUNIT7':'', 
              'TTYPE8':'TDIM7',
              'TFORM8': '16A', 
              'TUNIT8': '', 
              'TTYPE9':'TUNIT7', 
              'TFORM9':'6A',
              'TUNIT9':'',
              'TTYPE10':'CTYPE1',
              'TFORM10':'8A', 
              'TUNIT10':'Hz',
              'TTYPE11':'CRVAL1',
              'TFORM11': '1D',
              'TUNIT11': 'Hz',
              'TTYPE12': 'CRPIX1',             
              'TFORM12': '1D',
              'TUNIT12':'',
              'TTYPE13':' CDELT1',
              'TFORM13': '1D',
              'TUNIT13': 'Hz',
              'TTYPE14':'CTYPE2', 
              'TFORM14': '4A', 
              'TUNIT14': '',
              'TTYPE15':'CRVAL2',
              'TFORM15': '1D',
              'TUNIT15': 'deg',
              'TTYPE16': 'CTYPE3',
              'TFORM16':'4A',
              'TUNIT16':'',
              'TTYPE17':'CRVAL3', 
              'TFORM17':'1D',
              'TUNIT17':'deg',
              'TTYPE18':'CRVAL4', 
              'TFORM18':'1I',
              'TUNIT18':'',
              'TTYPE19':'OBSERVER',
              'TFORM19':'32A',
              'TUNIT19':'',
              'TTYPE20':'OBSID',
              'TFORM20':'32A', 
              'TUNIT20':'',
              'TTYPE21':'SCAN',
              'TFORM21':'1J',
              'TUNIT21':'',
              'TTYPE22':'OBSMODE',
              'TFORM22':'32A', 
              'TUNIT22':'',
              'TTYPE23':'FRONTEND',
              'TFORM23':'16A',
              'TUNIT23':'',
              'TTYPE24':'TCAL',
              'TFORM24':'1E',
              'TUNIT24':'K', 
              'TTYPE25':'VELDEF',
              'TFORM25':'8A', 
              'TUNIT25':'', 
              'TTYPE26':'VFRAME',
              'TFORM26':'1D',
              'TUNIT26':'m/s',
              'TTYPE27':'RVSYS',
              'TFORM27':'1D',
              'TUNIT27':'m/s', 
              'TTYPE28':'OBSFREQ',
              'TFORM28':'1D',
              'TUNIT28':'Hz',
              'TTYPE29':'LST', ##TODO: FIX SO PER INTEGRATION AND ACTUALLY CONVERT TO LST!!
              'TFORM29':'1D',
              'TUNIT29':'s',
              'TTYPE30':'AZIMUTH',
              'TFORM30':'1D',
              'TUNIT30':'deg',
              'TTYPE31':'ELEVATIO', 
              'TFORM31':'1D',
              'TUNIT31':'deg',
              'TTYPE32':'TAMBIENT',
              'TFORM32':'1D',
              'TUNIT32':'K', 
              'TTYPE33':'PRESSURE',
              'TFORM33':'1D',
              'TUNIT33':'mmHg', 
              'TTYPE34':'HUMIDITY',
              'TFORM34':'1D',
              'TUNIT34':'',
              'TTYPE35':'RESTFRQ',
              'TFORM35':'1D',
              'TUNIT35':'Hz',
              'TTYPE36':'FREQRES', ##TODO finish other COVMODES (Currently only PAF_CAL)
              'TFORM36':'1D',
              'TUNIT36':'Hz',
              'TTYPE37':'EQUINOX',
              'TFORM37':'1D',
              'TUNIT37':'',
              'TTYPE38':'RADESYS', 
              'TFORM38':'8A',
              'TUNIT38':'',
              'TTYPE39':'TRGTLONG',
              'TFORM39':'1D',
              'TUNIT39':'deg', 
              'TTYPE40':'TRGTLAT',
              'TFORM40':'1D',
              'TUNIT40':'deg',
              'TTYPE41':'SAMPLER',
              'TFORM41':'12A',
              'TUNIT41':'',  
              'TTYPE42':'FEED',
              'TFORM42':'1I',
              'TUNIT42':'',
              'TTYPE43':'SRFEED',
              'TFORM43':'1I',
              'TUNIT43':'',
              'TTYPE44':'FEEDXOFF', 
              'TFORM44':'1D',
              'TUNIT44':'deg',
              'TTYPE45':'FEEDEOFF', 
              'TFORM45': '1D',
              'TUNIT45':'deg',
              'TTYPE46':'SUBREF_STATE', 
              'TFORM46':'1I',
              'TUNIT46':'',
              'TTYPE47':'SIDEBAND',
              'TFORM47':'1A',
              'TUNIT47':'',
              'TTYPE48':'PROCSEQN',
              'TFORM48':'1I',
              'TUNIT48':'',
              'TTYPE49':'PROCSIZE',
              'TFORM49':'1I',
              'TUNIT49':'',
              'TTYPE50':'PROCSCAN',
              'TFORM50':'16A',
              'TUNIT50':'',
              'TTYPE51':'PROCTYPE',
              'TFORM51':'16A',
              'TUNIT51':'',
              'TTYPE52':'LASTON',
              'TFORM52':'1J',
              'TUNIT52':'',
              'TTYPE53':'LASTOFF',
              'TFORM53':'1J',
              'TUNIT53':'',
              'TTYPE54':'TIMESTAMP',
              'TFORM54':'22A',
              'TUNIT54':'UTC', 
              'TTYPE55':'QD_XEL',
              'TFORM55':'1D',
              'TUNIT55':'deg',
              'TTYPE56':'QD_EL', 
              'TFORM56':'1D',
              'TUNIT56':'deg',
              'TTYPE57':'QD_BAD',
              'TFORM57':'1I',
              'TUNIT57':'',
              'TTYPE58':'QD_METHOD', 
              'TFORM58':'1A', 
              'TUNIT58':'',
              'TTYPE59':'VELOCITY',
              'TFORM59':'1D', 
              'TUNIT59':'m/s',
              'TTYPE60':'ZEROCHAN',
              'TFORM60':'1E',
              'TUNIT60':'',
              'TTYPE61':'DOPFREQ', ## TODO only set to 0.0 for now; needs method for final production code. 
              'TFORM61':'1D', 
              'TUNIT61':'Hz', 
              'TTYPE62':'SIG',
              'TFORM62':'1A', 
              'TUNIT62':'',
              'TTYPE63':'CAL',
              'TFORM63':'1A',
              'TUNIT63':'',
              'TTYPE64':'CALTYPE',
              'TFORM64':'8A',
              'TUNIT64':'K',
              'TTYPE65':'TWARM',
              'TFORM65':'1E', 
              'TUNIT65':'K',
              'TTYPE66':'TCOLD',
              'TFORM66':'1E',
              'TUNIT66':'K',
              'TTYPE67':'CALPOSITION',
              'TFORM67':'16A',
              'TUNIT67':'',
              'TTYPE68':'IFNUM',
              'TFORM68':'1I',
              'TUNIT68':'',
              'TTYPE69':'PLNUM',
              'TFORM69':'1I',
              'TUNIT69':'',
              'TTYPE70':'FDNUM',
              'TFORM70':'1I',
              'TUNIT70':''
              }   
    
    def getCurrentUTC(self):
        time = datetime.datetime.utcnow()
        dateStr=str(time.year)
        dateStr+='-'+str(time.strftime('%m'))
        dateStr+='-'+str(time.strftime('%d'))
        dateStr+='T'
        dateStr+=str(time.strftime('%H'))
        dateStr+=':'+str(time.strftime('%M'))
        dateStr+=':'+str(time.strftime('%S'))
        return dateStr
        
    def readScanLog_Header(self):
        os.chdir(fitsPath)
        hdu = fits.open('ScanLog.fits')
        return hdu[0].header      
        
    def constructPriHDUHeader(self):
        ## Collect header metadata from ScanLog
        priHeader=fits.Header()
        currentUTC = self.getCurrentUTC()
        priHeader.set('DATE',currentUTC, 'date and time this HDU was created, UTC')
        scanLogHeader = self.readScanLog_Header()
        origin = scanLogHeader['ORIGIN']
        priHeader.set('ORIGIN',origin,'origin of observation')
        telescope = scanLogHeader['TELESCOP']
        priHeader.set('TELESCOP',telescope,'the telescope used')
        ##TODO:
        ##Decide on SDFITSVER and FITSVER keywords
        
        ## Set other metadata
        priHeader.set('INSTRUME','FLAGBF','backend')
        priHeader.set('SDFITVER','sdfits-bf','SDFITS format for BF')
        priHeader.set('FITSVER','fits-bf','FITS format for BF')
        prihdu = fits.PrimaryHDU(header=priHeader)
        return prihdu 

    def progressBar(self,value, endvalue, beam,bar_length=20):

        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))
        sys.stdout.write("\rPercent of FITS table for  beam "+np.str(beam)+": [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()     
    
    def getProjId(self):
        os.chdir(fitsPath)
        hdu = fits.open('ScanLog.fits')
        return hdu[0].header['PROJID']
    
    def constuctBinTableHeader(self):    
        ## construct primary HDU
        prihdu = self.constructPriHDUHeader()   
        binHeader = fits.Header()
        ## load list of required SDFITS keywords
        keywordList = np.loadtxt(keywordPath,dtype='bytes')
        keyWordArr = keywordList.astype(str)
        ##TODO ??
        commentList = []
        paramList = []
        
        ##TODO:SDFITS CORE KEYWORDS
        ## loop through each keyword
        for keyIdx in range(0,len(keyWordArr),3):
            keyword = keyWordArr[keyIdx]         
            ## utilize the keyword dictionary to get the associated parameter
            param = self.keyToParamDict[keyword]
            ## should this be here????
            self.progressBar(keyIdx, len(keyWordArr),self.beamNum)
            ## knowing the required parameter, use a dictionary of associated functions
            ## which provide required metadata information and place data
            self.funcDict[param.strip()](self,param.strip())  
            formKey = keyWordArr[keyIdx+1]
            unitKey = keyWordArr[keyIdx+2]
            form = self.keyToParamDict[formKey]
            if formKey == 'TFORM7':
                form = np.str(len(self.Column.valueArr[0,:]))+'E'
            unit = self.keyToParamDict[unitKey]
            col = fits.Column(name=self.Column.param,format=form,array=self.Column.valueArr,unit=unit)
            if keyIdx == 0:
                self.cols = fits.ColDefs([col])
            else:
                self.cols.add_col(col)
            commentList.append(self.Column.comment)
            paramList.append(self.Column.param)
        ##make preliminary table HDU     
        tblHdu = fits.BinTableHDU.from_columns(self.cols, header = binHeader)
        
        ##Now, update comments beginnning of table
        tblHdu.header.set('NAXIS',2, '2-dimensional binary table')
        rowBytes = tblHdu.size/(2*self.numInts)
        tblHdu.header.set('NAXIS1',int(rowBytes),'width of table in bytes')
        tblHdu.header.set('NAXIS2',self.numInts*2,'number of rows in table')
        tblHdu.header.set('PCOUNT',0,'size of special data area')
        tblHdu.header.set('GCOUNT',1,'one data group (required keyword)')
        tblHdu.header.set('TFIELDS',70,'number of fields in each row')
        tblHdu.header.set('EXTNAME','SINGLE DISH', 'name of this binary table extension')       
        idx = 0        
        for i in range(0,len(keyWordArr),3):
            keyword = keyWordArr[i]
            tblHdu.header.set(keyword,paramList[idx],commentList[idx])
            idx+=1
        ##insert final additional comments
        tblHdu.header.insert('TTYPE1',('COMMENT', 'Start of SDFITS CORE keywords/columns.'))
        tblHdu.header.insert('TTYPE2',('TELESCOP','NRAO_GBT','the telescope used'))
        tblHdu.header.insert('TTYPE7',('COMMENT', 'End of SDFITS CORE keywords/columns.'))
        tblHdu.header.insert('TTYPE7',('COMMENT', 'Start of SDFITS DATA column and descriptive axes.'))
        tblHdu.header.insert('TTYPE18',('CTYPE4','STOKES', 'fourth axis is Stokes'))
        tblHdu.header.insert('TTYPE19',('COMMENT', 'End of SDFITS DATA column and descriptive axes.'))
        tblHdu.header.insert('TTYPE19',('COMMENT', 'Start of SDFITS SHARED keywords/columns.'))
        projId = self.getProjId()    
        tblHdu.header.insert('TTYPE21',('PROJID',projId,'project identifier'))        
##        tblHdu.header.insert('TTYPE24',('BACKEND','FLAG','backend device'))
        tblHdu.header.insert('TTYPE35',('SITELONG', -7.983983E+01,'E. longitude of intersection of the az/el axes'))
        tblHdu.header.insert('TTYPE35',('SITELAT', 3.843312E+01,'N. latitude of intersection of the az/el axes'))
        tblHdu.header.insert('TTYPE35',('SITEELEV', 8.245950E+02,'height of the intersection of az/el axes'))
        tblHdu.header.insert('TTYPE37',('COMMENT', 'End of SDFITS SHARED keywords/columns.'))
        tblHdu.header.insert('TTYPE37',('COMMENT', 'Start of GBT-specific keywords/columns.'))
        tblHdu.header.insert('TTYPE44',('COMMENT', 'Feed offsets ARE included in the CRVAL2 and CRVAL3 columns'))
        tblHdu.header.insert('TFORM70',('COMMENT', 'End of GBT-specific keywords/columns.'))    
        
        thduList = fits.HDUList([prihdu, tblHdu])
        os.chdir(rawDataPath)   ##TODO: run on flag3; 
      
        return thduList
    
