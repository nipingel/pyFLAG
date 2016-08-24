# -*- coding: utf-8 -*-
"""
This module (i.e. class) will collate the required metadata to construc the primary HDU and SDFITS binary table. 

@author: npingel
"""
from astropy.io import fits
from astropy.time import Time
import glob
import os 
import numpy as np
import datetime
import collections


class MetaDataModule:
    
        

    def __init__(self,fitsList,numInts,dataBuff_X,dataBuff_Y,beamNum,intLen,numThreads):        
        self.fitsList = [fitsList]     
        for i in range(0,len(self.fitsList)):
            self.fitsList[i] = self.fitsList[i][:-7] + "4" + self.fitsList[i][-5:]           
        print(self.fitsList)
        ##TODO: remove for production
        self.numInts = numInts
        self.dataBuff_X = dataBuff_X
        self.dataBuff_Y = dataBuff_Y
        self.beamNum = beamNum 
        self.intLen = intLen
        self.numThreads = numThreads
        self.Column = collections.namedtuple('Column',['param','valueArr','comment'])
        self.cols = None

##TODO: check to make sure len(fitsList) / numThreads = 1in all methods
        def getSMKey(self,param,singleVal = False): ##TODO: get actual shared memory values

            valueArr = np.empty([self.numInts*len(self.fitsList)/self.numThreads],dtype='object')         
            if param == 'DURATION':
                param = 'REQSTI'
            elif param == 'EXPOSURE':
                param = 'ACTSTI'
            elif param == 'OBSMODE' and singleVal == True:
                param = 'COVMODE'
                corHDU = fits.open(fitsList[0])
                value = corHDU[0].header[param]
                return value
            elif param == 'OBSMODE':
                param = 'COVMODE'                    
            elif param == 'LST' or param == 'DATE-OBS':
                param1 = 'DATE-OBS'
                param2 = 'REQSTI'
                for file in range(0,len(self.fitsList)/self.numThreads):            
                    corHDU = fits.open(fitsList[file])
                    dateTime = corHDU[0].header[param1]
                    intLen = corHDU[0].header[param2]
                    dateTimeObj = Time(dateTime,format='isot',scale='utc')   
                    intLen = intLen/(24*3600.)
                    idx = 0
                    for i in range(0,self.numInts):
                        newTime = ((idx*self.numInts)+i)*intLen+dateTimeObj.jd
                        newTime = Time(newTime,format='isot', scale='utc')
                        value = newTime.isot
                        valueArr[(idx*self.numInts)+i] = value
                    idx+=1
            else:
                for file in range(0,len(self.fitsList)/self.numThreads):            
                    corHDU = fits.open(fitsList[file])
                    value = corHDU[0].header[param]
                    idx = 0
                    for i in range(0,self.numInts):
                        valueArr[(idx*self.numInts)+i] = value
                    idx+=1
            comment = self.commentDict[param]            
            self.Column.param = param            
            self.Column.comment = comment
            self.Column.valueArr = valueArr
            
        ##TODO: error handling (e.g. keyword not found)
        def getGOFITSParam(self,param):
            os.chdir('/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/GO') ##TODO: run on flag03          
            idx = 0  
            valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='object') 
            repeat = True            
            if param == 'TRGTLONG':
                paramLook = 'RA'
            elif param == 'TRGTLAT':
                paramLook = 'DEC'           
            elif param == 'CTYPE2' or param == 'CTYPE3':
                paramLook = 'COORDSYS'
            else:
                paramLook = param
            for file in range(0,len(self.fitsList)):            
                goHDU = fits.open(self.fitsList[file])                
                value = goHDU[0].header[paramLook]
                if value == 'GALACTIC' and param == 'CTYPE2':
                    value = 'GLON'      
                    valueArr.fill(value)
                    repeat = False
                elif value == 'GALACTIC' and param == 'CTYPE3':
                    value = 'GLAT'
                    valueArr.fill(value)
                    repeat = False
                elif value == 'RADEC' and param == 'CTYPE2':
                    value = 'RA'
                elif value == 'RADEC' and param == 'CTYPE3':
                    value = 'DEC'
                elif value == 'HADEC' and param == 'CTYPE2':
                    value = 'HA'
                elif value == 'HADEC' and param == 'CTYPE3':
                    value = 'DEC'
                elif value == 'AZEL' and param == 'CTYPE2':
                    value = 'AZ'
                elif value == 'AZEL' and param == 'CTYPE3':
                    value = 'EL'    
                elif param  == 'TIMESTAMP':
                    valStr = fitsList[file]
                    value = valStr[:-6] 
                elif repeat == True:
                    for i in range(0,self.numInts):
                        valueArr[(idx*self.numInts)+i] = value
                idx+=1
                    
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
        
        def getAntFITSParam(self,param):
            os.chdir('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/Antenna/') ##TODO: run on flag03 
            valueArr = np.empty([self.numInts*len(self.fitsList)/self.numThreads],dtype='object')           
            idx = 0  
            valueIdx = 0            
            repeat = True
            for file in range(0,len(self.fitsList)/self.numThreads):            
                antHDU = fits.open(fitsList[file])
                if param == 'DMJD':
                    valueArr[(valueIdx*self.numInts):(valueIdx*self.numInts)+self.numInts] = antHDU[2].data['DMJD']
                    repeat = False
                elif param == 'CRVAL2':
                    maj = antHDU[2].data['MAJOR']
                    samplesPerInt = len(maj)/self.numInts
                    remainder, integer = np.modf(samplesPerInt)
                    skip = int(integer)
                    intCnt = 0
                    for i in range(0,len(maj),skip):
                        firstSampBin = maj[i:(i+skip)]
                        nextSampBin = maj[(i+skip):(i+skip)+skip]
                        interpSamp = firstSampBin[-1]+(nextSampBin[0] - firstSampBin[-1])*remainder
                        finalBin = np.hstack((firstSampBin,interpSamp))
                        meanVal = np.mean(finalBin)
                        valueArr[intCnt] = meanVal
                        intCnt+=1
                        repeat = False
                elif param == 'CRVAL3':
                    minor = antHDU[2].data['MINOR']
                    samplesPerInt = len(minor)/self.numInts
                    remainder, integer = np.modf(samplesPerInt)
                    skip = int(integer)
                    intCnt = 0
                    for i in range(0,len(minor),skip):
                        firstSampBin = minor[i:(i+skip)]
                        nextSampBin = minor[(i+skip):(i+skip)+skip]
                        interpSamp = firstSampBin[-1]+(nextSampBin[0] - firstSampBin[-1])*remainder
                        finalBin = np.hstack((firstSampBin,interpSamp))
                        meanVal = np.mean(finalBin)
                        valueArr[intCnt] = meanVal
                        intCnt+=1
                        repeat = False
                
                elif param == 'AZIMUTH':
                    az = antHDU[2].data['MNT_AZ']
                    samplesPerInt = len(az)/self.numInts
                    remainder, integer = np.modf(samplesPerInt)
                    skip = int(integer)
                    intCnt = 0
                    for i in range(0,len(az),skip):
                        firstSampBin = az[i:(i+skip)]
                        nextSampBin = az[(i+skip):(i+skip)+skip]
                        interpSamp = firstSampBin[-1]+(nextSampBin[0] - firstSampBin[-1])*remainder
                        finalBin = np.hstack((firstSampBin,interpSamp))
                        meanVal = np.mean(finalBin)
                        valueArr[intCnt] = meanVal
                        intCnt+=1
                        repeat = False
                elif param == 'ELEVATIO':
                    el = antHDU[2].data['MNT_EL']
                    samplesPerInt = len(el)/self.numInts
                    remainder, integer = np.modf(samplesPerInt)
                    skip = int(integer)
                    intCnt = 0
                    for i in range(0,len(el),skip):
                        firstSampBin = el[i:(i+skip)]
                        nextSampBin = el[(i+skip):(i+skip)+skip]
                        interpSamp = firstSampBin[-1]+(nextSampBin[0] - firstSampBin[-1])*remainder
                        finalBin = np.hstack((firstSampBin,interpSamp))
                        meanVal = np.mean(finalBin)
                        valueArr[intCnt] = meanVal
                        intCnt+=1
                        repeat = False
                elif param == 'TAMBIENT':
                    param = 'AMBTEMP'
                elif param == 'HUMIDITY':
                    param = 'AMBHUMID'
                elif param == 'PRESSURE':
                    value = antHDU[0].header['AMBPRESS']
                if param == 'PRESSURE':
                    value=value*0.75006375541921                 
                elif param == 'TAMBIENT':
                    value+=273.0
                if repeat == 'True':
                    for i in range(0,self.numInts):
                        valueArr[(idx*self.numInts)+i] = value
                        idx+=1
                valueIdx+=1   
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
            
        def getRcvrFITSParam(self,param):
            os.chdir('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/Rcvr1_2/') ##TODO: run on flag03 
            valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='object')           
            idx = 0 
            RcvrHDU = fits.open('2005_05_27_00:00:00.fits')
            for file in range(0,len(self.fitsList)):            
                value = RcvrHDU[0].header[param]               
                for i in range(0,self.numInts):
                    valueArr[(idx*self.numInts)+i] = value
                idx+=1
                    
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment    
       
        def getArbParam(self,param):           
            if param == 'TSYS' or param == 'TCAL' or param == 'DOPFREQ':    
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='float32')                
                valueArr.fill(1.0)
            if param == 'BEAM':    
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='float32')                
                valueArr.fill(self.beamNum)
            elif param == 'FEED' or 'SUBREF_STATE' or 'QD_BAD':    
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='int16')                
                valueArr.fill(1)
            elif param == 'ZEROCHAN' or param =='TWARM' or param =='TCOLD' or param == 'QD_XEL' or param == 'QD_EL':
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='float32')                
                valueArr.fill(np.nan)
            elif param == 'FDNUM' or param == 'IFNUM':
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='int16')                 
                valueArr.fill(0.0)
            elif param == 'SIG' or param == 'CAL' or param == 'CALTYPE':
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='s1')                
                valueArr.fill('T')
            elif param == 'QD_METHOD':
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='s1')                
                valueArr.fill('C')
            elif param == 'TUNIT7':
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='s6')                
                valueArr.fill('counts')
            elif param == 'CTYPE1':
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='s8')                
                valueArr.fill('FREQ-OBS')
            elif param == 'SRFEED':
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='int16')
                valueArr.fill(0)
            elif param == 'SIDEBAND':
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='s1')
                valueArr.fill('L')
            
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
        
        def getDataParam(self,param):
            if param == 'DATA':
                valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='float32')
                valueArr[::2] = self.dataBuff_Y
                valueArr[1::2] = self.dataBuff_X            
            elif param == 'TDIM7':
                valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='float32')
                numChans = np.float(len(self.dataBuff_X[0,:,0]))
                valueArr.fill(numChans)
            elif param == 'PLNUM':
                valueArr = np.empty([2*self.numInts*len(self.fitsList)],dtype='int16')
                valueArr[::2] = 0
                valueArr[1::2] = 1
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
        
        def getLOFITSParam(self,param):
            repeat = True            
            if param == 'VELDEF':
                valueArr = np.zeros([self.numInts*len(self.fitsList)],dtype='float32')            
                tblNum = 2
            elif param == 'VFRAME' or param == 'RVSYS':
                valueArr = np.zeros([self.numInts*len(self.fitsList)],dtype='float32')            
                tblNum = 3
                repeat == False
            elif param == 'CRVAL1' or 'OBSFREQ':
                valueArr = np.zeros([self.numInts*len(self.fitsList)],dtype='float32')
                param = 'LO1FREQ'
                tblNum = 3
            os.chdir('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/LO1A/') ##TODO: run on flag3/GB
            idx=0            
            for file in range(0,len(self.fitsList)):            
                loHDU = fits.open(fitsList[file])       
                if repeat == True:        
                    value = loHDU[tblNum].header[param]
                    valueArr.fill(value)                    
                else:
                    for i in range(0,self.numInts):
                        valueCol = loHDU[tblNum].data[param]                       
                        value = valueCol[len(valueCol)/2] ##get val corresponding to middle DMJD
                        valueArr[(idx*self.numInts)+i] = value
                idx+=1
               
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
        
        def getBeamOffsets(self,param):
            hdu = fits.open('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/Sim2Data/weight.fits')
            if param == 'FEEDXOFF':
                beamOff_Az = hdu[1].data['BeamOff_AZ']
                value = beamOff_Az[self.beamName]
            else:
                beamOff_El = hdu[1].data['BeamOff_EL']
                value = beamOff_El[self.beamName]
            valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='float32')
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
                    valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='float32')
                    valueArr.fill(value)      
            elif param == 'CRPIX1':
                modeName = getSMKey(self,param,singleVal = True)
                if modeName == 'PAF_CAL':
                    value = 500/2.
                    valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='float32')
                    valueArr.fill(value)  
            elif param == 'BANDWID':
                modeName = getSMKey(self,param,singleVal = True)
                if modeName == 'PAF_CAL':
                    value = 303.18*500*1000.
                    valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='float32')
                    valueArr.fill(value) 
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
                            'BANDWID':'bandwidth'
                            
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
                         'TIMESTAMP':getAntFITSParam,
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
                         'LST': getSMKey,
                         'CTYPE2':getGOFITSParam,
                         'CTYPE3':getGOFITSParam,
                         'CRPIX1':getModeDepParams,
                         'BANDWID':getModeDepParams
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
              'TTYPE65':'TWARM',
              'TFORM65':'1E', 
              'TUNIT65':'K',
              'TTYPE66':'TCOLD',
              'TFORM66':'1E',
              'TUNIT66':'K',
              'TTYPE67':'CALPOSITION',
              'TFORM67':'16A' ,
              'TTYPE68':'IFNUM',
              'TFORM68':'1I',
              'TTYPE69':'PLNUM',
              'TFORM69':'1I',
              'TTYPE70':'FDNUM',
              'TFORM70':'1I',
              }   
    ## returns scanLog binary table
    def readScanLog_Data(self):
        scanLogHduList = fits.open('ScanLog.fits')
        return scanLogHduList[1].data            
                    
    ## returns scanLog header
    def readScanLog_Header(self):
        scanLogHduList = fits.open('ScanLog.fits')
        return scanLogHduList[0].header
    
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
        
    def contstructPriHDUHeader(self):
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
        return priHeader   

    def constuctBinTableHeader(self):    
        
        

        binHeader = fits.Header()
        keywordList = np.loadtxt('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/sdKeywords.txt',dtype='bytes')
        keyWordArr = keywordList.astype(str)

        binHeader.set('BITPIX','BINTABLE', 'binary table extension')
        binHeader.set('NAXIS',2, '2-dimensional binary table')
        binHeader.set('NAXIS1',33362,'width of table in bytes') ##TODO: update
        binHeader.set('PCOUNT',0,'size of special data area')
        binHeader.set('NAXIS2',1664, 'number of rows in table') ##TODO: update
        binHeader.set('GCOUNT',1,'one data group (required keyword)')
        binHeader.set('TFIELDS',70,'number of fields in each row') ##TODO: update
        binHeader['COMMENT'] = 'Start of SDFITS CORE keywords/columns.'
        ##TODO:SDFITS CORE KEYWORDS
        for keyIdx in range(0,len(keyWordArr)):
            keyword = keyWordArr[keyIdx]         
            if keyword[0:5] != 'TFORM' and keyword[0:5] != 'TUNIT':
                param = self.keyToParamDict[keyword]
                self.funcDict[param](self,param)  
                formKey = keyWordArr[keyIdx+1]
                unitKey = keyWordArr[keyIdx+2]
                form = self.keyToParamDict[formKey]
                if formKey == 'TFORM7':
                    form = np.str(len(self.valueArr))+'E'
                unit = self.keyToParamDict[unitKey]
                col = fits.Column(name=self.Column.param,format=form,array=self.Column.valueArr,unit=unit)    
                if keyIdx == 0:
                    self.cols = fits.ColDefs([col])
                else:
                    self.cols.add_col(col)
        
            
        binHeader['COMMENT'] = 'End of SDFITS CORE keywords/columns.'
        binHeader['COMMENT'] = 'Start of SDFITS DATA column and descriptive axes.'
        ##TODO: SDFITS DATA KEYWORDS (including Beamformer specific)
        binHeader['COMMENT'] = 'End of SDFITS DATA column and descriptive axes.'
        binHeader['COMMENT'] = 'Start of SDFITS SHARED keywords/columns.'
        ##TODO: SDFITS SHARED KEYWORDS
        binHeader['COMMENT'] = 'End of SDFITS SHARED keywords/columns.'
        binHeader['COMMENT'] = 'Start of GBT-specific keywords/columns.'
        ##TODO: GBT-SPECIFIC KEYWORDS
        binHeader['COMMENT'] = 'Feed offsets ARE included in the CRVAL2 and CRVAL3 columns.'
        ##TODO: MORE GBT-SPECIFIC KEYWORDS
        ## TODO: comment; find this 'PROJID'
        ##TODO: 'CTYPE4':'STOKES', ##TODO: comment
        ##TODO: 'BACKEND'
        binHeader['COMMENT'] = 'End of GBT-specific keywords/columns.'
        binHeader.set('EXTNAME','SINGLE DISH', 'name of this binary table extension')
        
        ##TODO:
        #'SITELONG':getAntFits, ## find and define this
        #'SITELAT':getAntFits,## find and define this
        #'SITEELEV':getAntFits,## find and define this        
        return binHeader
    def constructBinTableData(self):
        return
