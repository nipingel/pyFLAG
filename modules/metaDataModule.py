# -*- coding: utf-8 -*-
"""
This module (i.e. class) will collate the required metadata to construc the primary HDU and SDFITS binary table. 

@author: npingel
"""
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as pyplot
import pyslalib as pysla
import os 
import numpy as np
import datetime
import collections
import glob
from . import RadVelCorr
import sys


##globalPaths
keywordPath = '/users/npingel/FLAG/SpectralFiller/misc/sdKeywords.txt'

class MetaDataModule:
    
        
    ## initialize function. Opens raw data file, numInts, data buffers (banpasses), beam, freqChans, intLength
    ## total BANK files, and creates column object. 
    def __init__(self, projectPath, rawDataPath, fitsList, bankFitsList, numBanksList, dataBuff_X,dataBuff_Y,beamNum, pfb):        
        self.projectPath = projectPath
        ## list of scan time stamps 
        self.fitsList = fitsList  
        ## list containing all relevant BANK files
        self.bankFitsList = bankFitsList
        self.numBanksList = numBanksList
        self.dataBuff_X = dataBuff_X
        self.dataBuff_Y = dataBuff_Y
        self.beamNum = beamNum 
        self.dataPath = rawDataPath
        self.pfb = pfb
        self.valueArr = None
        self.newArr = None
        ## get numFreq Chans
        firstRow = dataBuff_X[0]
        self.numFreqChans = len(firstRow[0,:])
        self.Column = collections.namedtuple('Column',['param','valueArr','comment'])
        self.cols = None
        self.numPhases = 2
        ## arrays to hold coordinate transformation
        self.beamOff_CrossEl = None
        self.beamOff_El = None
        self.oldRaArr = None
        self.oldDecArr = None
        self.oldAzArr = None
        self.oldElArr = None
        self.lstArr = None
        self.dataDMJDList = []
        self.refractList = []
        self.GBTLAT  = np.deg2rad(38.4331294)
        self.GBTLONG = np.deg2rad(79.8398397)
        self.GBTHGT  = 824.36                     # meters above the ellipsoid
        ## initialize radial velocity correction module
        self.radvelcorrObj = RadVelCorr.RadVelCorr()
        ## initialize dictionaries
        self.commentDict = {'OBJECT':'name of source observed',
                    'OBSERVER':'name of observer(s)',                            
                    'BANDWID':'bandwidth',
                    'DATE-OBS':'date and time of observation start',
                    'DURATION':'total integration duration in seconds',
                    'EXPOSURE':'effective int time (excludes blanking) in secs',
                    'TSYS':'system temperature in Kelvin',                            
                    'OBSID':'observation description',
                    'SCAN':'scan number',	
                    'RESTFREQ':'rest frequency at band center',
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
        self.funcDict = {'OBJECT':self.getGOFITSParam,
                 'CTYPE1':self.getArbParam,
                 'CRVAL1':self.getLOFITSParam,
                 'DATE-OBS':self.getSMKey,
                 'DURATION':self.getSMKey,
                 'EXPOSURE':self.getSMKey,
                 'TSYS':self.getArbParam,
                 'OBSERVER':self.getGOFITSParam,
                 'OBSID':self.getGOFITSParam,
                 'SCAN':self.getGOFITSParam,
                 'RESTFREQ':self.getGOFITSParam, ## TODO: make this general. For now simply set to HI in TOPO frame
                 'EQUINOX':self.getGOFITSParam,
                 'RADESYS':self.getGOFITSParam,
                 'TRGTLONG':self.getGOFITSParam,
                 'TRGTLAT':self.getGOFITSParam,
                 'PROCSEQN':self.getGOFITSParam,
                 'PROCSIZE':self.getGOFITSParam,
                 'PROCSCAN':self.getGOFITSParam,
                 'PROCTYPE':self.getGOFITSParam,
                 'LASTON':self.getGOFITSParam,
                 'LASTOFF':self.getGOFITSParam,
                 'VELOCITY':self.getGOFITSParam,
                 'FRONTEND':self.getRcvrFITSParam,
                 'TAMBIENT':self.getAntFITSParam,
                 'PRESSURE':self.getAntFITSParam,
                 'HUMIDITY':self.getAntFITSParam,
                 'TIMESTAMP':self.getGOFITSParam,
                 'AZIMUTH':self.getAntFITSParam,
                 'ELEVATIO':self.getAntFITSParam,
                 'ZEROCHAN':self.getArbParam,
                 'SIG':self.getArbParam,
                 'TWARM':self.getArbParam,
                 'TCOLD':self.getArbParam,
                 'FDNUM':self.getArbParam,
                 'IFNUM':self.getArbParam,
                 'TUNIT7':self.getArbParam,
                 'TDIM7':self.getDataParam,
                 'DATA':self.getDataParam,
                 'PLNUM':self.getDataParam,
                 'OBSMODE':self.getGOFITSParam,
                 'TCAL':self.getArbParam,
                 'VELDEF':self.getLOFITSParam,
                 'VFRAME':self.getLOFITSParam,
                 'RVSYS':self.getLOFITSParam,
                 'OBSFREQ':self.getLOFITSParam,
                 'VRSYS':self.getLOFITSParam,
                 'CAL':self.getArbParam,
                 'CALTYPE':self.getArbParam,
                 'CALPOSITION':self.getArbParam,
                 'FEED':self.getArbParam,
                 'SAMPLER':self.getArbParam,
                 'SRFEED':self.getArbParam,
                 'FEEDXOFF':self.getBeamOffsets,
                 'FEEDEOFF':self.getBeamOffsets,
                 'SUBREF_STATE':self.getArbParam,
                 'SIDEBAND':self.getArbParam, 
                 'QD_XEL':self.getArbParam,
                 'QD_EL':self.getArbParam,
                 'QD_BAD':self.getArbParam,
                 'QD_METHOD':self.getArbParam,
                 'DOPFREQ':self.getArbParam, 
                 'FREQRES':self.getModeDepParams, 
                 'LST':self.getAntFITSParam,
                 'CTYPE2':self.getGOFITSParam,
                 'CTYPE3':self.getGOFITSParam,
                 'CRPIX1':self.getModeDepParams,
                 'CDELT1':self.getModeDepParams,
                 'BANDWID':self.getModeDepParams,
                 'CRVAL2':self.getAntFITSParam,
                 'CRVAL3':self.getAntFITSParam,
                 'CRVAL4':self.getArbParam, 
                 'TRGTLONG':self.getGOFITSParam,
                 'TRGTLAT':self.getGOFITSParam
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
              'TTYPE35':'RESTFREQ',
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

    ## function to initialize self.valueArr/self.newArr based on datatype and fileNum iteration
    def initArr(self, fileNum, numInts, dtype, valStr):
        ## if first iteration, initialize the global array to hold all ints info; 
        ## only scan phases are the two pols so self.numPhases = 2. Makes a character array if
        ## data types are not numeric
        if valStr == 'DATA':
            if fileNum == 0:
                self.valueArr = np.empty([self.numPhases * numInts, self.numFreqChans], dtype=dtype)
            else:
                self.newArr = np.empty([self.numPhases * numInts, self.numFreqChans], dtype=dtype)
        else:
            if fileNum == 0: 
                if not any([dtype == 'float64', dtype == 'int32', dtype == 'int16']):
                    self.valueArr = np.chararray([self.numPhases * numInts], itemsize=len(valStr))
                else:  
                    self.valueArr = np.empty([self.numPhases * numInts], dtype=dtype)
            ## update self.newArr if we are past the first iteration of the loop over scan FITS files
            else:
                if not any([dtype == 'float64', dtype == 'int32', dtype == 'int16']):
                    self.newArr = np.chararray([self.numPhases * numInts], itemsize=len(valStr))
                else:
                    self.newArr = np.empty([self.numPhases * numInts], dtype=dtype)

    ## function to get number of integrations in a scan 
    def getNumInts(self, idx):
            ## get associated row for DATA table. numInts 
            ## is then simply the length of the first dimension. 
            dataRow = self.dataBuff_Y[idx]
            numScanInts = len(dataRow[:,0])
            return numScanInts

    

    def getSMKey(self,param):
        ## iterate flag defaults to false. Is set True based on parameter
        iterate = False
        bankIdx = 0
        timeVal = False
        for fileNum in range(0,len(self.fitsList)):
            corrHDU = fits.open(self.bankFitsList[bankIdx])
            ## get number of ints knowing which scan we're working with
            numScanInts = self.getNumInts(fileNum)
            ## if numScanIts is 0, we have a bad file and should skip filling
            ## metadata. 
            if numScanInts == 0:
                ## update BANK index
                bankIdx += self.numBanksList[fileNum]
                continue
            ## determine time spent in state (e.g. calOn, calOff); 
            ## since we do no freq sw or cal sw, this is equal to scan time
            if param == 'DURATION':
                self.initArr(fileNum, numScanInts, 'float64', None)
                value = corrHDU[0].header['ACTSTI']
            ## retrieve actual integration time
            elif param == 'EXPOSURE':
                paramLook = 'ACTSTI'
                self.initArr(fileNum, numScanInts, 'float64', None)
                value = corrHDU[0].header['ACTSTI']
            elif param == 'DATE-OBS':
                scanDMJD = corrHDU[1].data['DMJD']
                ## cast DMJD to correct string format
                timeObj = Time(scanDMJD, format = 'mjd', scale='utc', precision=2)
                ## initialize array to store parameter values
                self.initArr(fileNum, numScanInts, 'str', timeObj[0].isot)
                timeVal = True
        
            ## fill in value array for similar values across integrations/polarizations
            if timeVal == True and fileNum == 0:
                for ind in range(0, numScanInts):
                    self.valueArr[ind * self.numPhases] = timeObj[ind].isot
                    self.valueArr[ind * self.numPhases + 1] = timeObj[ind].isot
            elif timeVal == True and fileNum > 0:
                for ind in range(0, numScanInts):
                    self.newArr[ind * self.numPhases] = timeObj[ind].isot
                    self.newArr[ind * self.numPhases + 1] = timeObj[ind].isot
                self.valueArr = np.concatenate([self.valueArr, self.newArr])
            elif timeVal == False and fileNum == 0:
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
        elif param == 'RESTFREQ':
            paramLook = 'RESTFRQ'
        else:
            paramLook = param
        bankIdx = 0
        for fileNum in range(0,len(self.fitsList)):            
            goHDU = fits.open(self.projectPath + '/GO/' + self.fitsList[fileNum])                
            corrHDU = fits.open(self.bankFitsList[bankIdx])
            ## get chansel
            chanSel = np.int(corrHDU[0].header['CHANSEL'])
            numScanInts = self.getNumInts(fileNum)
            if numScanInts == 0:
                ## update BANK index
                bankIdx += self.numBanksList[fileNum]
                continue
            if paramLook != 'TIMESTAMP' and paramLook != 'OBSMODE':           
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
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param =='TRGTLAT':
                valStr = value
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param  == 'TIMESTAMP':
                valStr = self.fitsList[fileNum]
                valStr = valStr[:-6]
                self.initArr(fileNum, numScanInts, 'str', valStr)
            elif param == 'OBSMODE':
                valStr = goHDU[0].header['PROCNAME']
                valStr = valStr + ':' + goHDU[0].header['SWSTATE']
                valStr= valStr + ':' + goHDU[0].header['SWTCHSIG']
                self.initArr(fileNum, numScanInts, 'str', valStr)
            elif param =='VELOCITY':
                valStr = value
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param =='OBJECT':
                valStr = value
                self.initArr(fileNum, numScanInts, 'str', valStr)
            elif param == 'OBSERVER':
                valStr = value
                self.initArr(fileNum, numScanInts, 'str', valStr)
            elif param == 'OBSID':
                valStr = 'unknown'
                self.initArr(fileNum, numScanInts, 'str', valStr)
            elif param == 'RESTFREQ':
                valStr = 1420.4057517667*1e6##TODO:remove for production
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param == 'SCAN' or param == 'PROCSEQN' or param == 'PROCSIZE':
                valStr = value
                self.initArr(fileNum, numScanInts, 'int32', None)
            elif param == 'PROCSCAN' or param == 'PROCTYPE':
                valStr = value
                self.initArr(fileNum, numScanInts, 'str', valStr)
            elif param =='EQUINOX':
                valStr = value
                self.initArr(fileNum, numScanInts, 'float64', None)
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
                self.valueArr = np.concatenate([self.valueArr, self.newArr])  
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
            antHDU = fits.open(self.projectPath + '/Antenna/' + self.fitsList[fileNum])
            ## open weight file to get cross-el/el offset
            weightFiles = glob.glob(self.dataPath + 'weight_files/*FullGrid.FITS')
            wHDU = fits.open(weightFiles[0])
            beamOff_Az = wHDU[1].data['BeamOFF_AZ'][self.beamNum]
            beamOff_El = wHDU[1].data['BeamOFF_EL'][self.beamNum]
            ## get pointing model information 
            smntAz = antHDU[0].header['SMNTC_AZ']
            sobscAz = antHDU[0].header['SOBSC_AZ']
            sobscEl = antHDU[0].header['SOBSC_EL']
            smntEl = antHDU[0].header['SMNTC_EL']
            
            ## calculate pointing model
            azPt = smntAz - sobscAz
            elPt = smntEl - sobscEl
            ## get integrations this scan
            corrHDU = fits.open(self.bankFitsList[bankIdx])
            corrDMJD = corrHDU[1].data['DMJD']
            ## set offsets if first iteration
            if fileNum == 0:
              self.beamOff_CrossEl = beamOff_Az
              self.beamOff_El = beamOff_El
            ## get Antenna DMJD values
            antDMJD = antHDU[2].data['DMJD']
            ## get integrations this scan
            numScanInts = self.getNumInts(fileNum)
            if numScanInts == 0:
                ## update BANK index
                bankIdx += self.numBanksList[fileNum]
                continue
            if param == 'CRVAL2':
                ## interpolate antenna position samples to data time samples
                maj = antHDU[2].data['MAJOR']
                value = np.interp(corrDMJD, antDMJD,maj)
                ## initialize array to hold parameter values if first iteration
                self.initArr(fileNum, numScanInts, 'float64', None)
                
                ## interpolate refraction samples to data time samples
                interpRefract = np.interp(corrDMJD, antDMJD, antHDU[2].data['REFRACT'])
                ## extend data DMJD list for later use when converting from 
                ## geocentric apparent to mean J2000.0 Ra/Dec (after beam offsets are applied)
                self.dataDMJDList.extend(corrDMJD)
                ## extend list containing refraction corrections for elevation. These are applied
                ## before converting from geocentric apparent to mean J2000.0
                self.refractList.extend(interpRefract)
            elif param == 'CRVAL3':
                minor = antHDU[2].data['MINOR']
                value = np.interp(corrDMJD,antDMJD,minor)
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param == 'AZIMUTH': 
                az = antHDU[2].data['MNT_AZ'] - azPt ## subtract pointing model contribution
                value = np.interp(corrDMJD,antDMJD,az)
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param == 'ELEVATIO':
                el = antHDU[2].data['MNT_EL'] - elPt ## subtract pointing model contribution
                value = np.interp(corrDMJD,antDMJD,el)
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param == 'LST':
                param1 = 'LSTSTART'
                lstStart = antHDU[0].header[param1]
                param2 = 'ACTSTI'
                intLen = np.float(corrHDU[0].header[param2])
                self.initArr(fileNum, numScanInts, 'float64', None)
                value = np.zeros(numScanInts, dtype='float64')
                for idx in range(0,numScanInts):
                    #value[i] = lstStart + intLen/2+(i*intLen)
                     val = lstStart + (idx*intLen)
                     value[idx] = val
            elif param == 'TAMBIENT':
                value = antHDU[0].header['AMBTEMP']
                value+=273.0
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param == 'HUMIDITY':
                value = antHDU[0].header['AMBHUMID']
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param == 'PRESSURE':
                value = antHDU[0].header['AMBPRESS']
                value=value*0.75006375541921    
                self.initArr(fileNum, numScanInts, 'float64', None)
            if fileNum == 0:
                self.valueArr[::2] = value
                self.valueArr[1::2] = value
            else: 
                self.newArr[::2] = value
                self.newArr[1::2] = value
                self.valueArr = np.concatenate([self.valueArr, self.newArr])
            bankIdx += self.numBanksList[fileNum]
        ## determine if we need to save valueArr for later correction 
        if param == 'CRVAL2':
          self.oldRaArr = self.valueArr
        elif param == 'CRVAL3':
          self.oldDecArr = self.valueArr
        elif param == 'AZIMUTH':
          self.oldAzArr = self.valueArr
        elif param == 'ELEVATIO':
          self.oldElArr = self.valueArr
        elif param == 'LST':
          self.lstArr = self.valueArr/3600*15

        ## retrieve comment and create column 
        comment = self.commentDict[param]
        self.Column.param = param
        self.Column.valueArr = self.valueArr
        self.Column.comment = comment
    
    ## function to write backend to data; will need to get from Manager FITS file at production    
    def getRcvrFITSParam(self,param):
        if param == 'FRONTEND':            
            valStr = 'PAF' ##TODO: search for value in PAF manager for production
        bankIdx = 0
        for fileNum in range(0,len(self.fitsList)):            
            corrHDU = fits.open(self.bankFitsList[bankIdx])
            numScanInts = self.getNumInts(fileNum)
            if numScanInts == 0:
                ## update BANK index
                bankIdx += self.numBanksList[fileNum]
                continue
            self.initArr(fileNum, numScanInts, 'str', valStr)
            if fileNum == 0:
                self.valueArr[:] = valStr
            else: 
                self.newArr[:] = valStr
                self.valueArr = np.concatenate([self.valueArr, self.newArr])
             
            bankIdx += self.numBanksList[fileNum]
                
        comment = self.commentDict[param]
        self.Column.param = param
        self.Column.valueArr = self.valueArr
        self.Column.comment = comment    
    ## function that fills column with a single value
    def getArbParam(self,param):
        bankIdx = 0
        ## flag for multiple values
        multiVal = False
        for fileNum in range(0,len(self.fitsList)):
            corrHDU = fits.open(self.bankFitsList[bankIdx])
            numScanInts = self.getNumInts(fileNum)
            if numScanInts == 0:
                ## update BANK index
                bankIdx += self.numBanksList[fileNum]
                continue          
            if param == 'TSYS' or param == 'TCAL' or param == 'DOPFREQ':    
                self.initArr(fileNum, numScanInts, 'float64', None)
                value = 1.0
            elif param == 'BEAM':    
                self.initArr(fileNum, numScanInts, 'float64', None)
                value = self.beamNum
            elif param == 'FEED' or param == 'SUBREF_STATE' or param == 'QD_BAD':
                self.initArr(fileNum, numScanInts, 'int16', None)
                value = 1
            elif param == 'CRVAL4':
                self.initArr(fileNum, numScanInts, 'int16', None)
                value1 = -6
                value2 = -5
                multiVal = True
            elif param == 'ZEROCHAN' or param =='TWARM' or param =='TCOLD' or param == 'QD_XEL' or param == 'QD_EL' or param == 'CALPOSITION':
                self.initArr(fileNum, numScanInts, 'float64', None)
                value = np.nan
            elif param == 'FDNUM' or param == 'IFNUM':
                self.initArr(fileNum, numScanInts, 'int16', None)
                value = 0
            elif param == 'SIG' or param == 'CAL' or param == 'CALTYPE':
                value = 'T'
                self.initArr(fileNum, numScanInts, 'str', value)
            elif param == 'QD_METHOD':
                value = 'C' 
                self.initArr(fileNum, numScanInts, 'str', value)
            elif param == 'TUNIT7':
                value = 'counts'
                self.initArr(fileNum, numScanInts, 'str', value)
            elif param == 'CTYPE1':
                value = 'FREQ-OBS'
                self.initArr(fileNum, numScanInts, 'str', value)
            elif param == 'SRFEED':
                value = 0
                self.initArr(fileNum, numScanInts, 'int16', None)
            elif param == 'SIDEBAND':
                value = 'L'
                self.initArr(fileNum, numScanInts, 'str', value)
            elif param =='SAMPLER':
                value = 'A1_0'    
                self.initArr(fileNum, numScanInts, 'str', value)
            elif param == 'DOPFREQ':
                value = 1450*1e6 ##TODO: fix for production. 
                self.initArr(fileNum, numScanInts, 'float64', None)

            if multiVal == True and fileNum == 0:
                self.valueArr[0::2] = value1
                self.valueArr[1::2] = value2
            elif multiVal == True and fileNum > 0:
                self.newArr[0::2] = value1
                self.newArr[1::2] = value2
                self.valueArr = np.concatenate([self.valueArr, self.newArr])
            elif multiVal == False and fileNum == 0:
                self.valueArr[:] = value
            else:
                self.newArr[:] = value
                ## concatenate new array to global value array
                self.valueArr = np.concatenate([self.valueArr, self.newArr])

            bankIdx += self.numBanksList[fileNum]
   
        comment = self.commentDict[param]
        self.Column.param = param
        self.Column.valueArr = self.valueArr
        self.Column.comment = comment
    ## function to process keywords related to the actual data
    def getDataParam(self,param):
        bankIdx = 0
        multiVal = False
        for fileNum in range(0,len(self.fitsList)):
            corrHDU = fits.open(self.bankFitsList[bankIdx])
            numScanInts = self.getNumInts(fileNum)
            if numScanInts == 0:
                ## update BANK index
                bankIdx += self.numBanksList[fileNum]
                continue
            if param == 'DATA':
                self.initArr(fileNum, numScanInts, 'float64', 'DATA')
                value1 = self.dataBuff_Y[fileNum]
                value2 = self.dataBuff_X[fileNum]
                self.initArr(fileNum, numScanInts, 'float64', 'DATA')
                multiVal = True
            elif param == 'TDIM7':
                value = '['+np.str(self.numFreqChans)+',1,1,1]'
                self.initArr(fileNum, numScanInts, 'str', value)
            elif param == 'PLNUM':
                value1 = 0 
                value2 = 1
                multiVal = True
                self.initArr(fileNum, numScanInts, 'int16', None) 
            
            if multiVal == True and fileNum == 0:
                 if param == 'DATA':
                    self.valueArr[::2,:] = value1
                    self.valueArr[1::2,:] = value2
                 else:
                    self.valueArr[::2] = value1
                    self.valueArr[1::2] = value2
            elif multiVal == True and fileNum > 0:
                if param == 'DATA':
                    self.newArr[::2,:] = value1
                    self.newArr[1::2,:] = value2
                else:
                    self.newArr[::2] = value1
                    self.newArr[1::2] = value2
                self.valueArr = np.concatenate([self.valueArr, self.newArr])
            elif multiVal == False and fileNum == 0:
                self.valueArr[:] = value
            else:
                self.newArr[:] = value
                self.valueArr = np.concatenate([self.valueArr, self.newArr])

            bankIdx += self.numBanksList[fileNum] 
         
        comment = self.commentDict[param]
        self.Column.param = param
        self.Column.valueArr = self.valueArr
        self.Column.comment = comment
    
    def getLOFITSParam(self,param): 
        bankIdx = 0
        for fileNum in range(0,len(self.fitsList)):
            corrHDU = fits.open(self.bankFitsList[bankIdx])
            ## get chansel
            chanSel = np.int(corrHDU[0].header['CHANSEL'])
            numScanInts = self.getNumInts(fileNum)
            if numScanInts == 0:
                ## update BANK index
                bankIdx += self.numBanksList[fileNum]
                continue
            if param == 'VELDEF':
                value = 'OPTI-OBS'
                self.initArr(fileNum, numScanInts, 'str', value)
            elif param == 'VFRAME' or param == 'RVSYS':
                value = 0
                self.initArr(fileNum, numScanInts, 'int16', None)
            elif param == 'CRVAL1' or 'OBSFREQ':
                value = 1450.00*1e6##TODO:remove for production
                if self.pfb == True:
                    lowEnd = value - (250-(chanSel*100))*.30318*1e6
                    highEnd = value-(250-(chanSel*100+100))*.30318*1e6
                    value = highEnd - (highEnd-lowEnd)/2
                self.initArr(fileNum, numScanInts, 'float64', None)
            
            if fileNum == 0:
                 self.valueArr[:] = value
            else:
                self.newArr[:] = value
                self.valueArr = np.concatenate([self.valueArr, self.newArr])
            bankIdx += self.numBanksList[fileNum]              
        
        comment = self.commentDict[param]
        self.Column.param = param
        self.Column.valueArr = self.valueArr
        self.Column.comment = comment
    
    def getBeamOffsets(self,param):
        bankIdx = 0
        for fileNum in range(0,len(self.fitsList)):
            corrHDU = fits.open(self.bankFitsList[bankIdx])
            numScanInts = self.getNumInts(fileNum)
            if numScanInts == 0:
                ## update BANK index
                bankIdx += self.numBanksList[fileNum]
                continue
            weightFiles = glob.glob(self.dataPath + 'weight_files/*FullGrid.FITS') 
            wHDU = fits.open(weightFiles[0])
            if param == 'FEEDXOFF':
                beamOff_Az = wHDU[1].data['BeamOff_AZ']
                value = beamOff_Az[self.beamNum]
                self.initArr(fileNum, numScanInts, 'float64', None)
            else:
                beamOff_El = wHDU[1].data['BeamOff_EL']
                value = beamOff_El[self.beamNum]
                self.initArr(fileNum, numScanInts, 'float64', None)
            if fileNum == 0:
                self.valueArr[:] = value
            else:
                self.newArr[:] = value
                self.valueArr = np.concatenate([self.valueArr, self.newArr])
            bankIdx += self.numBanksList[fileNum]

        comment = self.commentDict[param]
        self.Column.param = param
        self.Column.valueArr = self.valueArr
        self.Column.comment = comment
    
    def getModeDepParams(self,param):
        bankIdx = 0
        for fileNum in range(0,len(self.fitsList)):
            corrHDU = fits.open(self.bankFitsList[bankIdx])
            numScanInts = self.getNumInts(fileNum)
            if numScanInts == 0:
                ## update BANK index
                bankIdx += self.numBanksList[fileNum]
                continue
            modeName = corrHDU[0].header['MODENAME']
            if param == 'CDELT1' or param == 'FREQRES':
                if modeName == 'FLAG_CALCORR_MODE':
                    value = 303.18*1000.
                elif modeName == 'FLAG_PFBCORR_MODE':
                    value = 303.18*5/160.*1000.
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param == 'CRPIX1':
                if modeName == 'FLAG_CALCORR_MODE':
                    value = value = 500/2   
                elif modeName == 'FLAG_PFBCORR_MODE':
                    value = 3200/2
                self.initArr(fileNum, numScanInts, 'float64', None)
            elif param == 'BANDWID':
                if modeName == 'FLAG_CALCORR_MODE':
                    value = 303.18*500*1000.
                elif modeName == 'FLAG_PFBCORR_MODE':
                    value = 100*303.18*1000
                self.initArr(fileNum, numScanInts, 'float64', None)
            if fileNum == 0:
                self.valueArr[:] = value
            else:
                self.newArr[:] = value
                self.valueArr = np.concatenate([self.valueArr, self.newArr])
            bankIdx += self.numBanksList[fileNum]
        
        comment = self.commentDict[param]
        self.Column.param = param
        self.Column.valueArr = self.valueArr
        self.Column.comment = comment
               
 
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
        hdu = fits.open(self.projectPath + '/ScanLog.fits')
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
        hdu = fits.open(self.projectPath + '/ScanLog.fits')
        return hdu[0].header['PROJID']

    def az2ra(self,LST, Az, El, dec):
        ## convert all angles to radians
        azRad = np.deg2rad(360-Az) ## Azimuth of NCP is 0 for GBT. This code assumes 180. 
        elRad = np.deg2rad(El)
        decRad = np.deg2rad(dec)
        ha = -1*np.arccos(1/np.cos(decRad)*(np.sin(elRad)*np.cos(self.GBTLAT) - np.cos(elRad)*np.cos(azRad)*np.sin(self.GBTLAT)))
        return LST + np.rad2deg(ha)

    def el2dec(self, LST, Az, El):
        ## convert all angles to radians
        azRad = np.deg2rad(360 - Az) ## Azimuth of NCP is 0 for GBT. This code assumes 180. 
        elRad = np.deg2rad(El)
        return np.rad2deg(np.arcsin(np.sin(elRad)*np.sin(self.GBTLAT) + np.cos(elRad)*np.cos(azRad)*np.cos(self.GBTLAT)))

    def offsetCorrection(self, hdu):
        newAzArr = np.zeros(len(self.oldAzArr))
        newElArr = np.zeros([len(self.oldElArr)])
        newRaArr = np.zeros([len(self.oldElArr)])
        newDecArr = np.zeros([len(self.oldElArr)])
        """
        for coordIdx in range(0, len(newAzArr)):
          ## make conversion from cross-el to azimuth
          azOffVal = self.beamOff_CrossEl/np.cos(np.deg2rad(self.oldElArr[coordIdx]))
          elOffVal = self.beamOff_El
          newAzArr[coordIdx] = self.oldAzArr[coordIdx] + azOffVal
          newElArr[coordIdx] = self.oldElArr[coordIdx] + elOffVal

          ## now, compute the new ra/dec with az/el/lst in hand...
          newDecArr[coordIdx] = self.el2dec(self.lstArr[coordIdx], newAzArr[coordIdx], newElArr[coordIdx])
          newRaArr[coordIdx] =  self.az2ra(self.lstArr[coordIdx], newAzArr[coordIdx], newElArr[coordIdx], newDecArr[coordIdx])
          ## add systemic offset
          newRaArr = newRaArr + sysRaOff
          newDecArr = newDecArr + sysDecOff
          if newRaArr[coordIdx] < 0:
            newRaArr[coordIdx] = newRaArr[coordIdx] + 360
        """
        ## we have refract value for a single pol
        ## extend by 2x to describe both XX and YY pol
        extRefractArr = np.zeros(len(self.refractList)*2)
        extRefractArr[0::2] = self.refractList
        extRefractArr[1::2] = self.refractList 
        azOffVal = self.beamOff_CrossEl/np.cos(np.deg2rad(self.oldElArr)) #- extRefractArr))
        elOffVal = self.beamOff_El
        newAzArr = self.oldAzArr + azOffVal
        newElArr = self.oldElArr - extRefractArr + elOffVal
        newDecArr = self.el2dec(self.lstArr, newAzArr, newElArr)
        newRaArr = self.az2ra(self.lstArr, newAzArr, newElArr, newDecArr)
        if any(t < 0 for t in newRaArr):
          newRaArr = newRaArr + 360
        ## we have data DMJD for a single pol 
        ## extend by 2x to describe both XX and YY pol 
        dataDMJDArr = np.zeros(len(self.dataDMJDList)*2)
        dataDMJDArr[0::2] = self.dataDMJDList
        dataDMJDArr[1::2] = self.dataDMJDList 
        for coordIdx in range(0, len(newRaArr)):
           raVal = newRaArr[coordIdx]
           decVal = newDecArr[coordIdx]
           dmjdVal = dataDMJDArr[coordIdx]
           ## convert from geocentric apparant to mean place (J2000.0) 
           newRa, newDec = pysla.slalib.sla_amp(np.deg2rad(raVal), np.deg2rad(decVal), dmjdVal, 2000.0)
           ## convert radians to degrees
           newRaArr[coordIdx] = np.rad2deg(newRa)
           newDecArr[coordIdx] = np.rad2deg(newDec)
        #print(np.mean(newDecArr - self.oldDecArr))
        #print(np.mean(newRaArr - self.oldRaArr))
        #pyplot.plot(self.oldDecArr, label='Measured')
        #pyplot.plot(newDecArr, label='Calculated')
        #pyplot.xlabel('Coordinate Element')
        #pyplot.ylabel('RA [deg]')
        #pyplot.legend(loc=0)
        #pyplot.savefig('MeasuredVCalculated_Ra_Dec.pdf')
        #pyplot.show()
        
        ## update header values
        hdu.data['CRVAL2'] = newRaArr
        hdu.data['CRVAL3'] = newDecArr
        hdu.data['TRGTLONG'] = newRaArr
        hdu.data['TRGTLAT'] = newDecArr
        hdu.data['AZIMUTH'] = newAzArr
        hdu.data['ELEVATIO'] = newElArr
        return hdu

    def radVelCorrection(self, hdu):
        raArr = hdu.data['CRVAL2']
        decArr = hdu.data['CRVAL3']
        utDate_Time = hdu.data['DATE-OBS']
        cenFreqsArr = hdu.data['CRVAL1']
        restFreqArr = hdu.data['RESTFREQ']
        c = 299792458.0 ## m/s
        for velIter in range(0,len(raArr)):
            raVal = raArr[velIter]
            decVal = decArr[velIter]
            utDateTimeVal = utDate_Time[velIter]
            cenFreqVal = cenFreqsArr[velIter]
            restFreqVal = restFreqArr[velIter]
            t = Time(utDateTimeVal, format='isot')
            utdate = t.iso[0:10]
            uttime = t.iso[11:]
            radVelCorrection = self.radvelcorrObj.correctVel(utdate,uttime,raVal,decVal)
            ## compute optical velocity of ref freq
            vOpt = c*(1-cenFreqVal/restFreqVal)
            ## add radial correction
            newVOpt = vOpt + radVelCorrection
            ## now convert back to frequency
            newCenFreq = (1-newVOpt/c)*restFreqVal 
            ## update reference frequency value
            cenFreqsArr[velIter] = newCenFreq
        ## update columns
        hdu.data['CRVAL1'] = cenFreqsArr
        hdu.data['OBSFREQ'] = cenFreqsArr
        hdu.data['VELDEF'] = 'OPTI-BAR'
        return hdu    

    def constuctBinTableHeader(self):    
        ## construct primary HDU
        prihdu = self.constructPriHDUHeader()   
        binHeader = fits.Header()
        ## load list of required SDFITS keywords
        keywordList = np.loadtxt(keywordPath,dtype='str')
        keyWordArr = keywordList.astype(str)
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
            self.funcDict[param.strip()](param.strip())  
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
        corrHDU = fits.open(self.bankFitsList[0])
        ## correct spatial offsets
        tblHdu = self.offsetCorrection(tblHdu)
        ## with updated RA/Dec, perform doppler correction
        tblHdu = self.radVelCorrection(tblHdu)       
        ## loop through globalBuffer list to get total number of scans
        totalInts = 0
        for idx in range(0,len(self.fitsList)):
            dataRow = self.dataBuff_Y[idx]
            totalInts += len(dataRow[:,0])    

        ##Now, update comments beginnning of table
        tblHdu.header.set('NAXIS',2, '2-dimensional binary table')
        rowBytes = tblHdu.size/(2 * totalInts)
        tblHdu.header.set('NAXIS1',int(rowBytes),'width of table in bytes')
        tblHdu.header.set('NAXIS2',totalInts * 2,'number of rows in table')
        tblHdu.header.set('PCOUNT',0,'size of special data area')
        tblHdu.header.set('GCOUNT',1,'one data group (required keyword)')
        tblHdu.header.set('TFIELDS',70,'number of fields in each row')
        tblHdu.header.set('EXTNAME','SINGLE DISH', 'name of this binary table extension')       
        tblHdu.header.set('XTENSION', 'BINTABLE', 'binary table extension')
        idx = 0        
        for i in range(0,len(keyWordArr),3):
            keyword = keyWordArr[i]
            tblHdu.header.set(keyword,paramList[idx],commentList[idx])
            idx+=1
        ## insert final additional comments
        tblHdu.header.insert('TTYPE1',('COMMENT', 'Start of SDFITS CORE keywords/columns.'))
        tblHdu.header.insert('TTYPE2',('TELESCOP','NRAO_GBT','the telescope used'))
        tblHdu.header.insert('TTYPE7',('COMMENT', 'End of SDFITS CORE keywords/columns.'))
        tblHdu.header.insert('TTYPE7',('COMMENT', 'Start of SDFITS DATA column and descriptive axes.'))
        tblHdu.header.insert('TTYPE18',('CTYPE4','STOKES', 'fourth axis is Stokes'))
        tblHdu.header.insert('TTYPE19',('COMMENT', 'End of SDFITS DATA column and descriptive axes.'))
        tblHdu.header.insert('TTYPE19',('COMMENT', 'Start of SDFITS SHARED keywords/columns.'))
        projId = self.getProjId()    
        tblHdu.header.insert('TTYPE21',('PROJID',projId,'project identifier'))        
        ## tblHdu.header.insert('TTYPE24',('BACKEND','FLAG','backend device'))
        tblHdu.header.insert('TTYPE35',('SITELONG', -7.983983E+01,'E. longitude of intersection of the az/el axes'))
        tblHdu.header.insert('TTYPE35',('SITELAT', 3.843312E+01,'N. latitude of intersection of the az/el axes'))
        tblHdu.header.insert('TTYPE35',('SITEELEV', 8.245950E+02,'height of the intersection of az/el axes'))
        tblHdu.header.insert('TTYPE37',('COMMENT', 'End of SDFITS SHARED keywords/columns.'))
        tblHdu.header.insert('TTYPE37',('COMMENT', 'Start of GBT-specific keywords/columns.'))
        tblHdu.header.insert('TTYPE44',('COMMENT', 'Feed offsets ARE included in the CRVAL2 and CRVAL3 columns'))
        tblHdu.header.insert('TFORM70',('COMMENT', 'End of GBT-specific keywords/columns.'))    
        
        

        thduList = fits.HDUList([prihdu, tblHdu])
      
        return thduList
    
