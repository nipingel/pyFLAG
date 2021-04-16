# -*- coding: utf-8 -*-
"""
01/09/2019
Module containing methods/functions for collating metadata associated with a FLAG scan. Several key
aspects about the SDFITS file under creation are given upon instantiation. These include:
Inputs from PAF_Filler.py:
projectPath - path to ancillary FITS files
rawDataPath - path to raw BANK FITS files
weightPath - path to weight FITS files
fitsList - list of fits associated with the observed object and beam
bankFitsList - list of BANK FITS files used in the processing thus far
numBanksList - list of the number of BANK FITS files 
dataBuff_X - global data buffer for the sorted beamformed XX polarization data
dataBuff_Y - global data buffer for the sorted beamformed YY polarization data
beamNum - beam number that is being processed
pfb - boolean flag that is True when we are processing data taken in FLAG_PFB_MODE 
firstIter - boolean flag that is True when processing the first binary FITS table, so a Primary HDU is necessary

PAF_Filler.py will call the method constructBinHeader(), which will drive the creation of the SDFITS files that
contains all of the necessary SDFITS keywords for use in the GBTIDL computing environment. The method, constructBinHeader(), 
ultimately returns an HDUList that PAF_Filler.py writes out as an SDFITS file to disk. 

Eample:
__email__ = ""Nickolas.Pingel@anu.edu.au"
__status__ = "Production"
"""
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
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
global keywordPath
keywordPath = '/users/npingel/FLAG/SpectralFiller/misc/sdKeywords.txt'
## make beam dictionary to map from BYU to WVU convention (e.g. BYU 0 -> WVU 1)
wvuBeamDict = {'0':'1', '1':'2', '2':'6', '3':'0', '4':'3', '5':'5','6':'4'}
byuBeamDict = {'1':'0', '2':'1', '6':'2', '0':'3', '3':'4', '5':'5', '4':'6'}

class MetaDataModule:
    
  ## initialize function. Opens raw data file, numInts, data buffers (banpasses), beam, freqChans, intLength
  ## total BANK files, and creates column object. 
  def __init__(self, projectPath, rawDataPath, weightPath, fitsList, restfreq, centralfreq, bankFitsList, numBanksList, dataBuff_X,dataBuff_Y,beamNum, pfb, firstIter):        
    self.projectPath = projectPath
    self.dataPath = rawDataPath
    self.weightPath = weightPath
    ## list of scan time stamps 
    self.fitsList = fitsList  
    self.restfreq = restfreq
    self.centralfreq = centralfreq
    ## list containing all relevant BANK files
    self.bankFitsList = bankFitsList
    self.numBanksList = numBanksList
    self.dataBuff_X = dataBuff_X
    self.dataBuff_Y = dataBuff_Y
    self.beamNum = np.int(beamNum) 
    self.pfb = pfb
    self.firstIter = firstIter
    self.valueArr = None
    self.newArr = None
    ## get numFreq Chans
    firstRow = dataBuff_X[0]
    self.numFreqChans = len(firstRow[0,:])
    self.Column = collections.namedtuple('Column',['param','valueArr','comment'])
    self.cols = None
    self.numPhases = 2
    ## arrays to hold coordinate transformation
    self.coordSysList = []
    self.offCoordSysList = []
    self.beamOff_CrossEl = None
    self.beamOff_El = None
    self.newRaArr = None
    self.newDecArr = None
    self.raList  = []
    self.decList = []
    self.mntAzList = []
    self.mntElList = []
    self.obscAzList = []
    self.obscElList = []
    self.smntAzList = []
    self.smntElList = []
    self.sobscAzList = []
    self.sobscElList = []
    self.lstArr = None
    self.dataDMJDList = []
    self.refractList = []
    self.polMotionXList = []
    self.polMotionYList = []
    self.utCorrList = []
    self.tempList = []
    self.pressList = []
    self.humList = []
    self.waveLengthList = []
    self.lapseRateList = []
    self.GBTLAT  = 38.4331294 ## deg
    self.GBTLONG = 79.8398397 ## deg
    self.GBTHGT  = 824.36                     # meters above the ellipsoid
    self.c = 299792458.0 ## m/s
    ## initialize radial velocity correction module
    self.radvelcorrObj = RadVelCorr.RadVelCorr()
    
    ## initialize dictionaries

    ## dictionary to associate parameter with header comments
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
    ## dictionary to associate parameter to collection function
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
      'CTYPE2':self.getAntFITSParam,
      'CTYPE3':self.getAntFITSParam,
      'CRPIX1':self.getModeDepParams,
      'CDELT1':self.getModeDepParams,
      'BANDWID':self.getModeDepParams,
      'CRVAL2':self.getAntFITSParam,
      'CRVAL3':self.getAntFITSParam,
      'CRVAL4':self.getArbParam, 
      'TRGTLONG':self.getGOFITSParam,
      'TRGTLAT':self.getGOFITSParam
      }

    ## dictionary to associate SDFITS keywords to specific parameter values
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

  ## method to initialize self.valueArr/self.newArr attributes based on datatype and fileNum iteration
  def initArr(self, fileNum, numInts, dtype, valStr):
      """
      if first iteration, initialize the global array to hold all ints info; 
      only scan phases are the two pols so self.numPhases = 2. Makes a character array if
      data types are not numeric.
      """
      if valStr == 'DATA': ## if we are sorting the actual beamformed spectra, we need a 2D array (ints*phase by freq channels)
        if fileNum == 0:
          self.valueArr = np.empty([self.numPhases * numInts, self.numFreqChans], dtype=dtype)
        else:
          self.newArr = np.empty([self.numPhases * numInts, self.numFreqChans], dtype=dtype)
      else: ## otherwise, a 1D vector will do. 
        if fileNum == 0: ## if first iteration, initialze the valueArr to store metadata values
          if not any([dtype == 'float64', dtype == 'int32', dtype == 'int16']):
            self.valueArr = np.chararray([self.numPhases * numInts], itemsize=len(valStr))
          else:  
            self.valueArr = np.empty([self.numPhases * numInts], dtype=dtype)
        else: ## update self.newArr if we are past the first iteration of the loop over scan FITS files. self.newArr will be appended to self.valueArr
              ## as subsequent scans are looped over
          if not any([dtype == 'float64', dtype == 'int32', dtype == 'int16']):
            self.newArr = np.chararray([self.numPhases * numInts], itemsize=len(valStr))
          else:
            self.newArr = np.empty([self.numPhases * numInts], dtype=dtype)

  ## utility method to get number of integrations in a scan 
  def getNumInts(self, idx):
    ## get associated row for DATA table. numInts 
    ## is then simply the length of the first dimension. 
    dataRow = self.dataBuff_Y[idx]
    numScanInts = len(dataRow[:,0])
    return numScanInts

  """  
  method to populate a FITS column based on the given parameter that is queried from the shared memory written out
  from the data acquistion stage
  """
  def getSMKey(self, param):
    ## iterate flag defaults to false. Is set True based on parameter
    iterate = False
    bankIdx = 0
    timeVal = False
    for fileNum in range(0,len(self.fitsList)):    ## loop through scan FITS files
      corrHDU = fits.open(self.bankFitsList[bankIdx]) ## open up raw data bank FITS files to get number of ints per scan
      numScanInts = self.getNumInts(fileNum) ## get number of ints knowing which scan we're working with
      
      """
      if numScanIts is 0, we have a bad file and should skip filling
      metadata.
      """ 
      if numScanInts == 0:
        ## update BANK index (skip to next scan)
        bankIdx += self.numBanksList[fileNum]
        continue
      
      """
      determine time spent in state (e.g. calOn, calOff); 
      since we do no freq sw or cal sw, this is equal to scan time
      """
      if param == 'DURATION':
        self.initArr(fileNum, numScanInts, 'float64', None) ## initialize array to store parameter values
        value = corrHDU[0].header['ACTSTI']
      elif param == 'EXPOSURE': ## retrieve actual integration time
        paramLook = 'ACTSTI'
        self.initArr(fileNum, numScanInts, 'float64', None)
        value = corrHDU[0].header['ACTSTI']
      elif param == 'DATE-OBS': ## These values are NOT constant, we will have to calculate them... set timeVal flag to True
        scanDMJD = corrHDU[1].data['DMJD'] ## retrieve DMJD values for scan and cast to correct string format
        timeObj = Time(scanDMJD, format = 'mjd', scale='utc', precision=2)
        self.initArr(fileNum, numScanInts, 'str', timeObj[0].isot) ## initialize array to store parameter values
        timeVal = True
      
      ## if first iteration, place modified DMJD into the array that will become the FITS column 
      if timeVal == True and fileNum == 0:
        for ind in range(0, numScanInts):
          self.valueArr[ind * self.numPhases] = timeObj[ind].isot
          self.valueArr[ind * self.numPhases + 1] = timeObj[ind].isot
      ## if greater than first iteration, place modified DMJD into a local data buffer to be added to array that will become FITS column
      elif timeVal == True and fileNum > 0:
        for ind in range(0, numScanInts):
          self.newArr[ind * self.numPhases] = timeObj[ind].isot
          self.newArr[ind * self.numPhases + 1] = timeObj[ind].isot
        self.valueArr = np.concatenate([self.valueArr, self.newArr]) ## add to array that will become the FITS column
      elif timeVal == False and fileNum == 0:
        self.valueArr[:] = value ## value is constant, so set every element equal to one another; place in FITS column array
      else:
        self.newArr[:] = value ## add to local data buffer to be added to FITS column array
        self.valueArr = np.concatenate([self.valueArr, self.newArr]) ## concatenate new array to global value array

      bankIdx += self.numBanksList[fileNum] ## update BANK index

     
    ## retrieve comment and construct column   
    comment = self.commentDict[param]      
    self.Column.param = param            
    self.Column.comment = comment
    self.Column.valueArr = self.valueArr
  
  """  
  method to populate a FITS column based on the given parameter that is queried from the GO FITS file written out
  by Astrid
  """  
  def getGOFITSParam(self,param):

    ## interpret the input parameter and associate it with a FITS keyword
    if param == 'TRGTLONG':
      paramLook = 'RA'
    elif param == 'TRGTLAT':
      paramLook = 'DEC'           
    elif param == 'CTYPE2' or param == 'CTYPE3':
      paramLook = 'COORDSYS'
    elif param == 'RESTFREQ':
      paramLook = 'RESTFRQ'
    else:
      paramLook = param ## parameter exists in current form
    
    bankIdx = 0 ## index to keep track of which BANK file we are on
    for fileNum in range(0,len(self.fitsList)):  ## loop through scan FITS files          
      goHDU = fits.open(self.projectPath + '/GO/' + self.fitsList[fileNum] + '.fits')             
      corrHDU = fits.open(self.bankFitsList[bankIdx]) ## open bank FITS file to get CHANSEL
      chanSel = np.int(corrHDU[0].header['CHANSEL']) ## get CHANSEL
      numScanInts = self.getNumInts(fileNum) ## determine the number of integrations
      """
      if numScanIts is 0, we have a bad file and should skip filling
      metadata.
      """ 
      if numScanInts == 0:
        ## update BANK index
        bankIdx += self.numBanksList[fileNum]
        continue ## move on to next scan
      """
      Based on the parameter, the data will be grabbed from the header.
      Some parameters are constant, while some require calculation. 
      """
      if paramLook != 'TIMESTAMP' and paramLook != 'OBSMODE':    
        ## try-except to catch an empty EQUINOX value (when mapping in other than J2000)
        try:       
          value = goHDU[0].header[paramLook] ## get parameter value from GO FITS file
        except KeyError:
          value = 2000.0
      else:
        value = None
      if param =='TRGTLONG': ## source longitude coordinate
        valStr = value
        self.initArr(fileNum, numScanInts, 'float64', None)
      elif param =='TRGTLAT': # source latitude coordinate
        valStr = value
        self.initArr(fileNum, numScanInts, 'float64', None)
      elif param  == 'TIMESTAMP': ## time stamp of scan 
        valStr = self.fitsList[fileNum]
        valStr = valStr
        self.initArr(fileNum, numScanInts, 'str', valStr)
      elif param == 'OBSMODE': ## observing mode
        valStr = goHDU[0].header['PROCNAME']
        valStr = valStr + ':' + goHDU[0].header['SWSTATE']
        valStr= valStr + ':' + goHDU[0].header['SWTCHSIG']
        self.initArr(fileNum, numScanInts, 'str', valStr)
      elif param =='VELOCITY': ## source velocity 
        valStr = value
        self.initArr(fileNum, numScanInts, 'float64', None)
      elif param =='OBJECT': ## observed Object
        valStr = value
        self.initArr(fileNum, numScanInts, 'str', valStr)
      elif param == 'OBSERVER': ## name of Observer
        valStr = value
        self.initArr(fileNum, numScanInts, 'str', valStr)
      elif param == 'OBSID': ## Project name
        valStr = 'unknown'
        self.initArr(fileNum, numScanInts, 'str', valStr)
      elif param == 'RESTFREQ': ## Rest frequency 
        valStr = self.restfreq
        self.initArr(fileNum, numScanInts, 'float64', None)
      elif param == 'SCAN' or param == 'PROCSEQN' or param == 'PROCSIZE': ## Scan number, sequence in procedures, or total scans in procedure
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
          self.coordSys = coordSys
          self.initArr(fileNum, numScanInts, 'str', valStr)
        elif coordSys == 'GALACTIC':
          self.coordSys = coordSys
          valStr = value
          self.initArr(fileNum, numScanInts, 'float64', valStr)
        else:
          valStr=np.str(value)
          self.initArr(fileNum, numScanInts, 'str', valStr)
      ## fill initial array
      if fileNum == 0:
        self.valueArr[:] = valStr
      else:
        self.newArr[:] = valStr
        self.valueArr = np.concatenate([self.valueArr, self.newArr])  ## concatenate subsequent arrays
      ## update bankIdx
      bankIdx += self.numBanksList[fileNum]  
    
    ## retrieve comment and construct column 
    comment = self.commentDict[param]
    self.Column.param = param
    self.Column.valueArr = self.valueArr
    self.Column.comment = comment

  """  
  method to populate a FITS column based on the given parameter that is queried from the Antenna FITS file written out
  by the servos
  """  
  def getAntFITSParam(self,param):                      
    
    bankIdx = 0 ## index to keep track of which BANK file we are on
    for fileNum in range(0,len(self.fitsList)): ## loop through scan FITS files               
      antHDU = fits.open(self.projectPath + '/Antenna/' + self.fitsList[fileNum] + '.fits')
      numScanInts = self.getNumInts(fileNum) ## get integrations this scan

      ## open weight file to get cross-el/el offset
      weightFiles = glob.glob(self.weightPath)
      wHDU = fits.open(weightFiles[0])
      
      ## get beam offsets
      beamOff_XEl = wHDU[1].data['BeamOFF_XEL'][self.beamNum]
      beamOff_El = wHDU[1].data['BeamOFF_EL'][self.beamNum]
          
      ## get integrations this scan
      corrHDU = fits.open(self.bankFitsList[bankIdx])
      corrDMJD = corrHDU[1].data['DMJD']
      
      ## set offsets if first iteration
      if fileNum == 0:
        self.beamOff_CrossEl = beamOff_XEl
        self.beamOff_El = beamOff_El
      
      antDMJD = antHDU[2].data['DMJD'] ## get Antenna DMJD values
      """
      if numScanIts is 0, we have a bad file and should skip filling
      metadata.
      """ 
      if numScanInts == 0:
        ## update BANK index
        bankIdx += self.numBanksList[fileNum]
        continue ## move on to next scan 
      
      """
      Since Antenna samples at 10 Hz, we need to interpolate antenna samples to data samples
      """
      if param == 'CRVAL2':
        maj = antHDU[2].data['MAJOR'] ## get antenna values for MAJOR axis coordinate
        value = np.interp(corrDMJD, antDMJD, maj) ## interpolate Antenna data to raw correlation data DMJD using Antenna DMJD
        
        ## initialize array to hold parameter values if first iteration
        self.initArr(fileNum, numScanInts, 'float64', None)
        interpRefract = np.interp(corrDMJD, antDMJD, antHDU[2].data['REFRACT']) ## interpolate refraction samples to data time samples
        """
        extend data DMJD list for later use when converting from 
        geocentric apparent to mean J2000.0 Ra/Dec (after beam offsets are applied)
        """
        self.dataDMJDList.extend(corrDMJD)
        """
        extend list containing refraction corrections for elevation. These are applied
        before converting from geocentric apparent to mean J2000.0
        """
        self.refractList.extend(interpRefract)

        """
        repeat interpolate for all other coordinates
        """
      elif param == 'CRVAL3':
        minor = antHDU[2].data['MINOR']
        value = np.interp(corrDMJD,antDMJD,minor)
        self.initArr(fileNum, numScanInts, 'float64', None)
      elif param == 'AZIMUTH': 
        mntAz = antHDU[2].data['MNT_AZ']
        obscAz = antHDU[2].data['OBSC_AZ']
        ## place pointing model information in list
        smntAz = antHDU[0].header['SMNTC_AZ']
        self.smntAzList.extend([smntAz] * numScanInts)
        sobscAz = antHDU[0].header['SOBSC_AZ']
        self.sobscAzList.extend([sobscAz] * numScanInts)

        value = np.interp(corrDMJD,antDMJD, mntAz)
        self.mntAzList.extend(value)
        self.obscAzList.extend(np.interp(corrDMJD, antDMJD, obscAz))
        self.initArr(fileNum, numScanInts, 'float64', None)
      elif param == 'ELEVATIO':
        el = antHDU[2].data['MNT_EL']
        obscEl = antHDU[2].data['OBSC_EL']

        ## place pointing model information in list
        smntEl = antHDU[0].header['SMNTC_EL']
        self.smntElList.extend([smntEl] * numScanInts)
        sobscEl = antHDU[0].header['SOBSC_EL']
        self.sobscElList.extend([sobscEl] * numScanInts)

        value = np.interp(corrDMJD,antDMJD,el)
        self.mntElList.extend(value)
        self.obscElList.extend(np.interp(corrDMJD, antDMJD, obscEl))
        self.initArr(fileNum, numScanInts, 'float64', None)
        ## now, add in the J2000 coordinates for beam offset applications later on
        antRa = antHDU[2].data['RAJ2000']
        antDec = antHDU[2].data['DECJ2000']
        interpRa = np.interp(corrDMJD, antDMJD, antRa)
        interpDec = np.interp(corrDMJD, antDMJD, antDec)
        self.raList.extend(interpRa)
        self.decList.extend(interpDec)
      elif param == 'CTYPE2' or param == 'CTYPE3':
        value = antHDU[0].header['INDICSYS']
        ## create temporary list to hold the scan's indicated coordinates system and extend the object list
        if param == 'CTYPE2':
          self.coordSysList.extend([value] * numScanInts)
          self.offCoordSysList.extend([antHDU[0].header['OFFSYS']] * numScanInts)
        if value == 'GALACTIC' and param == 'CTYPE2':
          valStr = 'GLON'
          ## initialize array to hold values if first iteration
          self.initArr(fileNum, numScanInts, 'str', valStr) 
        if value == 'GALACTIC' and param == 'CTYPE3':
          valStr = 'GLAT'
          self.initArr(fileNum, numScanInts, 'str', valStr)
        elif value == 'RADEC' and param == 'CTYPE2': ## Specify that we're in J2000 Coordinates
          valStr = 'RA'
          self.initArr(fileNum, numScanInts, 'str', valStr)
        elif value == 'RADEC' and param == 'CTYPE3':
          valStr = 'DEC'
          self.initArr(fileNum, numScanInts, 'str', valStr)
        elif value == 'HADEC' and param == 'CTYPE2': ## Specify that we're in Hour Angle/Dec Coordinates
          valStr = 'HA'
          self.initArr(fileNum, numScanInts, 'str', valStr)
        elif value == 'HADEC' and param == 'CTYPE3':
          valStr = 'DEC'
          self.initArr(fileNum, numScanInts, 'str', valStr) ## Specify that we're in Horizontal Coordinates
        elif value == 'AZEL' and param == 'CTYPE2':
          valStr = 'AZ'
          self.initArr(fileNum, numScanInts, 'str', valStr)
        elif value == 'AZEL' and param == 'CTYPE3':
          valStr = 'EL'
          self.initArr(fileNum, numScanInts, 'str', valStr)
        elif value == 'AZEL_ENCODER' and param == 'CTYPE2':
          valStr = 'AZ'
          self.initArr(fileNum, numScanInts, 'str', valStr)
        elif value == 'AZEL_ENCODER' and param == 'CTYPE3':
          valStr = 'EL'
          self.initArr(fileNum, numScanInts, 'str', valStr)
        value = valStr
        
      ## we will need to calculate this based on LST start time, DMJD value, and scan DATE-OBS (start time) -> all stored in Antenna FITS file
      elif param == 'LST': ## we will need to calculate this based on LST start and integration time
        param1 = 'LSTSTART' ## get LSTSTART value
        lstStart = antHDU[0].header[param1]
        scanStartStr = antHDU[0].header['DATE-OBS']
        #startObj = Time(scanStartStr, format = 'isot', scale = 'utc')
        startObj = Time(scanStartStr, format = 'fits')
        startMJD = startObj.mjd
        dmjdArr = np.array(self.dataDMJDList)
        lstArr = (1.00273790935 * (corrDMJD - startMJD)*86400.0) + lstStart
        self.initArr(fileNum, numScanInts, 'float64', None) ## intialize array
        #value = np.zeros(numScanInts, dtype='float64')
        value = np.copy(lstArr)
      elif param == 'TAMBIENT': ## air temperature 
        value = antHDU[0].header['AMBTEMP']
        value += 273.15 ## set to units of Kelvin
        self.initArr(fileNum, numScanInts, 'float64', None)
        ## because we're dealing with the Antenna FITS file, fill out lists needed for beam offsets later on...
        polMotionX = antHDU[0].header['IERSPMX'] ## radians
        polMotionY = antHDU[0].header['IERSPMY'] ## radians
        utCorr = antHDU[0].header['DELTAUTC'] ## sec

        ## fill lists with the value and set length based on number of integrations
        tempPolXList = [polMotionX] * len(self.valueArr)
        tempPolYList =  [polMotionY] * len(self.valueArr)
        tempUtCorrList = [utCorr] * len(self.valueArr)

        self.polMotionXList.extend(tempPolXList)
        self.polMotionYList.extend(tempPolYList)
        self.utCorrList.extend(tempUtCorrList)

        ## extend object lists with temporary lists created above
      elif param == 'HUMIDITY': ## humidity 
        value = antHDU[0].header['AMBHUMID']
        self.initArr(fileNum, numScanInts, 'float64', None)
      elif param == 'PRESSURE': ## air pressure
        value = antHDU[0].header['AMBPRESS']
        value=value*0.75006375541921    
        self.initArr(fileNum, numScanInts, 'float64', None)

      if fileNum == 0: 
        self.valueArr[::2] = value ## fill FITS column array
        self.valueArr[1::2] = value
      else: 
        self.newArr[::2] = value ## fill local data buffer to be added to FITS column array
        self.newArr[1::2] = value
        self.valueArr = np.concatenate([self.valueArr, self.newArr]) ## add to FITS column array
      bankIdx += self.numBanksList[fileNum] ## upate BANK index 

    ## retrieve comment and create column 
    comment = self.commentDict[param]
    self.Column.param = param
    self.Column.valueArr = self.valueArr
    self.Column.comment = comment
    
  """  
  method to populate a FITS column based on the given parameter that is queried from the FITS file written out
  by PAF Manager. Currently, the PAF Manager is removed and basically a dummy program that Astrid can communicate with. 
  This method will likely need work once fully integrated into the GBO system. Most of these keywords do not matter beyond
  their mere presence in the SDFITS binary table. 
  """     
  def getRcvrFITSParam(self,param):
    if param == 'FRONTEND':            
          valStr = 'PAF' 
    bankIdx = 0
    for fileNum in range(0,len(self.fitsList)): ## loop through scan FITS files              
      corrHDU = fits.open(self.bankFitsList[bankIdx]) ## open raw data FITS file
      numScanInts = self.getNumInts(fileNum) ## get number of integrations per scan
      
      """
      if numScanIts is 0, we have a bad file and should skip filling
      metadata.
      """ 
      if numScanInts == 0:
        ## update BANK index
        bankIdx += self.numBanksList[fileNum]
        continue ## move on to process next scan 
      self.initArr(fileNum, numScanInts, 'str', valStr)
      if fileNum == 0:
        self.valueArr[:] = valStr ## fill FITS column array
      else: 
        self.newArr[:] = valStr ## fill local data buffer to be added to FITS column array
        self.valueArr = np.concatenate([self.valueArr, self.newArr]) ## add to FITS column array
         
      bankIdx += self.numBanksList[fileNum] ## update bank index
    
    ## retrieve comment and create column         
    comment = self.commentDict[param]
    self.Column.param = param
    self.Column.valueArr = self.valueArr
    self.Column.comment = comment    
  
  """  
  method to populate a FITS column based with a single arbitrary value. The actual data values do not matter in the data reduction stage. 
  """  
  def getArbParam(self,param):
    bankIdx = 0 ## set bank index to keep try of what bank file we process/open
    multiVal = False ## flag for multiple values
    for fileNum in range(0,len(self.fitsList)): ## loop of scan FITS files
      corrHDU = fits.open(self.bankFitsList[bankIdx]) ## open raw data FITS file
      numScanInts = self.getNumInts(fileNum) ## get number of integrations per scan

      """
      if numScanIts is 0, we have a bad file and should skip filling
      metadata.
      """ 
      if numScanInts == 0:
        ## update BANK index
        bankIdx += self.numBanksList[fileNum]
        continue ## move on to next scan          
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
        value = 'T' ## set to True since we do not have multiple phases (on polarizations)
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
        value = 1449.8943*1e6 ##TODO: fix for production. 
        self.initArr(fileNum, numScanInts, 'float64', None)
      ## if we have multiple values, they are different for XX/YY Pol
      if multiVal == True and fileNum == 0: ## fill FITS column array
          self.valueArr[0::2] = value1
          self.valueArr[1::2] = value2
      elif multiVal == True and fileNum > 0: ## fill array to be added to FITS column array
          self.newArr[0::2] = value1
          self.newArr[1::2] = value2
          self.valueArr = np.concatenate([self.valueArr, self.newArr]) ## add to FITS column array
      elif multiVal == False and fileNum == 0:
          self.valueArr[:] = value ## fill FITS column array with single value
      else:
          self.newArr[:] = value ## fill array to be added to FITS column with a single value 
          self.valueArr = np.concatenate([self.valueArr, self.newArr]) ## add to FITS column array
      bankIdx += self.numBanksList[fileNum] ## update bank index

    ## retrieve comment and create column 
    comment = self.commentDict[param]
    self.Column.param = param
    self.Column.valueArr = self.valueArr
    self.Column.comment = comment

  """  
  method to populate a FITS column with DATA specific keywords. Much these parameters
  are in the raw bank FITS files.
  """  
  def getDataParam(self,param):
    bankIdx = 0 ## set bank index to keep try of what bank file we process/open
    multiVal = False ## flag for multiple values
    for fileNum in range(0,len(self.fitsList)): ## loop through scan FITS files
      corrHDU = fits.open(self.bankFitsList[bankIdx]) ## open raw bank FITS file
      numScanInts = self.getNumInts(fileNum)
     
      """
      if numScanIts is 0, we have a bad file and should skip filling
      metadata.
      """ 
      if numScanInts == 0:
        ## update BANK index
        bankIdx += self.numBanksList[fileNum]
        continue ## move on to next scan    
      if param == 'DATA':
        self.initArr(fileNum, numScanInts, 'float64', 'DATA') ## initialize 2D array (ints*phases x freqChan)
        value1 = self.dataBuff_Y[fileNum] ## with with scan's XX beamformed data
        value2 = self.dataBuff_X[fileNum] ## with with scan's YY beamformed data
        multiVal = True
      elif param == 'TDIM7':
        value = '['+np.str(self.numFreqChans)+',1,1,1]'
        self.initArr(fileNum, numScanInts, 'str', value) ## initialize 1D FITS column 
      elif param == 'PLNUM':
        value1 = 0 ## YY Pol
        value2 = 1 ## XX Pol
        multiVal = True
        self.initArr(fileNum, numScanInts, 'int16', None) 
      """
      if data for parameter has varying values, fill FITS column array based on polarization
      """
      if multiVal == True and fileNum == 0:
        if param == 'DATA': ## if actual beamformed spectra, fill the
          self.valueArr[::2,:] = value1 ## fill the YY Pol beamformed data along freq axis
          self.valueArr[1::2,:] = value2 ## fill the XX Pol beamformed data along freq axis
        else:
          self.valueArr[::2] = value1 ## fill the YY Pol data (1D)
          self.valueArr[1::2] = value2 ## fill the XX Pol data (1D)
        """
        if data for parameter has varying values, fill column to be added to FITS array based on polarization
        """
      elif multiVal == True and fileNum > 0:
        if param == 'DATA':
          self.newArr[::2,:] = value1 ## fill the YY Pol beamformed data along freq axis
          self.newArr[1::2,:] = value2 ## fill the XX Pol beamformed data along freq axis
        else:
          self.newArr[::2] = value1 ## fill the YY Pol data (1D)
          self.newArr[1::2] = value2 ## fill the XX Pol data (1D)
        self.valueArr = np.concatenate([self.valueArr, self.newArr]) ## add to FITS column array
        """
        if data for parameter does not vary and it is the first iteration, fill FITS column array with single value
        """
      elif multiVal == False and fileNum == 0:
        self.valueArr[:] = value
        """
        if data for paramter does not vary and it is past the first iteration, fill array to be added to FITS column array
        """
      else:
        self.newArr[:] = value
        self.valueArr = np.concatenate([self.valueArr, self.newArr])

      bankIdx += self.numBanksList[fileNum] 
     ## retrieve comment and create column  
    comment = self.commentDict[param]
    self.Column.param = param
    self.Column.valueArr = self.valueArr
    self.Column.comment = comment

  """  
  method to populate a FITS column with LO specific keywords. 
  """  
  def getLOFITSParam(self,param): 
    bankIdx = 0 ## set bank index to keep try of what bank file we process/open
    for fileNum in range(0,len(self.fitsList)): ## loop through scan FITS files
      corrHDU = fits.open(self.bankFitsList[bankIdx]) ## open raw bank FITS file
      chanSel = np.float(corrHDU[0].header['CHANSEL'])
      numScanInts = self.getNumInts(fileNum)
      #LOHdu = fits.open(self.projectPath + '/LO1B/' + self.fitsList[0] + '.fits')

      """
      if numScanIts is 0, we have a bad file and should skip filling
      metadata.
      """ 
      if numScanInts == 0:
        ## update BANK index
        bankIdx += self.numBanksList[fileNum]
        continue ## move on to next scan   

      if param == 'VELDEF': ## set velocity definition to OPTICAL-TOPOCENTRIC; this will be updated later
        value = 'OPTI-OBS'
        self.initArr(fileNum, numScanInts, 'str', value) ## initialize 1D FITS column vector
      elif param == 'VFRAME' or param == 'RVSYS': ## irrelevant keywords for now
        value = 0
        self.initArr(fileNum, numScanInts, 'int16', None)
        """
        Since only the LO1B FITS file is available, the RESTFRQ in primary header 
        corresponds to the TESTTONE FREQUENCY. As of Winter 2018, this is always 10 MHz ABOVE the LO frequency setting
        (in the TOPOCENTRIC reference frame). Really, this is the central frequency
        """
      elif param == 'CRVAL1' or param == 'OBSFREQ': 
        paramLook = 'RESTFRQ'
        #value = LOHdu[0].header[paramLook] - 10e6 ## subtract 10 MHz
        value = self.centralfreq
        if self.pfb == True: ## if we are in PFB mode, the chanSel is used to determine the frequency values for fine channels
            lowEnd = value - (250-(chanSel*100))*.30318*1e6
            highEnd = value- (250-(chanSel*100+100))*.30318*1e6
            value = highEnd - (highEnd-lowEnd)/2
        self.initArr(fileNum, numScanInts, 'float64', None) ## initialize FITS column array
      if fileNum == 0:
        self.valueArr[:] = value ## fill FITS column array
      else:
        self.newArr[:] = value ## fill array to be added FITS column array
        self.valueArr = np.concatenate([self.valueArr, self.newArr]) ## add to FITS column array
      bankIdx += self.numBanksList[fileNum]              
    ## retrieve comment and create column 
    comment = self.commentDict[param]
    self.Column.param = param
    self.Column.valueArr = self.valueArr
    self.Column.comment = comment
    
  """  
  method to populate a FITS column with offsets for beam being processed.
  """  
  def getBeamOffsets(self,param):
    bankIdx = 0 ## set bank index to keep try of what bank file we process/open
    for fileNum in range(0,len(self.fitsList)): ## loop through scan FITS files
      corrHDU = fits.open(self.bankFitsList[bankIdx]) ## open raw bank FITS file
      numScanInts = self.getNumInts(fileNum)

      """
      if numScanIts is 0, we have a bad file and should skip filling
      metadata.
      """ 
      if numScanInts == 0:
        ## update BANK index
        bankIdx += self.numBanksList[fileNum]
        continue ## move on to next scan   

      weightFiles = glob.glob(self.weightPath) ## open weight file 
      wHDU = fits.open(weightFiles[0]) ## grab first weight file (all values in headers are similar)
      if param == 'FEEDXOFF': 
        beamOff_XEl = wHDU[1].data['BeamOff_XEL'] ## get Cross-El offsets (in degs)
        value = beamOff_XEl[self.beamNum]
        self.initArr(fileNum, numScanInts, 'float64', None) ## initialize FITS column array
      else:
        beamOff_El = wHDU[1].data['BeamOff_EL'] ## get El offsets (in degs)
        value = beamOff_El[self.beamNum]
        self.initArr(fileNum, numScanInts, 'float64', None)
      if fileNum == 0: ## if first iteration, fill FITS column vector
          self.valueArr[:] = value
      else:
        self.newArr[:] = value ## if not first iteration, fill array to be added to FITS column array
        self.valueArr = np.concatenate([self.valueArr, self.newArr]) ## add to FITS column array
      bankIdx += self.numBanksList[fileNum] ## update bank index

    ## retrieve comment and create column 
    comment = self.commentDict[param]
    self.Column.param = param
    self.Column.valueArr = self.valueArr
    self.Column.comment = comment
    
  """  
  method to populate a FITS column with PAF mode dependent metadata
  """  
  def getModeDepParams(self,param):
    bankIdx = 0 ## set bank index to keep try of what bank file we process/open
    for fileNum in range(0,len(self.fitsList)): ## loop through scan FITS files
      corrHDU = fits.open(self.bankFitsList[bankIdx]) ## open raw bank FITS file
      numScanInts = self.getNumInts(fileNum)

      """
      if numScanIts is 0, we have a bad file and should skip filling
      metadata.
      """ 
      if numScanInts == 0:
        ## update BANK index
        bankIdx += self.numBanksList[fileNum]
        continue ## move on to next scan   

      modeName = corrHDU[0].header['MODENAME'] ## get mode name
      if param == 'CDELT1' or param == 'FREQRES':
        if modeName == 'FLAG_CALCORR_MODE':
          value = 303.18*1000. ## calculate freq resolution for calibration mode
        elif modeName == 'FLAG_PFBCORR_MODE':
          value = 303.18*5/160.*1000. ## calculate freq resolution for fine channel mode
        self.initArr(fileNum, numScanInts, 'float64', None) ## initialize FITS column array
      elif param == 'CRPIX1':
        if modeName == 'FLAG_CALCORR_MODE':
          value = value = 500/2 ## calculate reference pixel for calibration mode
        elif modeName == 'FLAG_PFBCORR_MODE':
          value = 3200/2 ## calculate reference pixel for fine channel mode
        self.initArr(fileNum, numScanInts, 'float64', None)
      elif param == 'BANDWID':
        if modeName == 'FLAG_CALCORR_MODE':
          value = 303.18*500*1000. ## calculate bandwith for calibraiton mode
        elif modeName == 'FLAG_PFBCORR_MODE':
          value = 100*303.18*1000 ## calculate bandwith for fine channel mode
        self.initArr(fileNum, numScanInts, 'float64', None)
      if fileNum == 0: ## if first iteration, fill FITS column array
        self.valueArr[:] = value
      else: ## if not first iteration, fill array to be added to FITS column array
        self.newArr[:] = value
        self.valueArr = np.concatenate([self.valueArr, self.newArr]) ## add to FITS column array
      bankIdx += self.numBanksList[fileNum] ## update bank index
  
    ## retrieve comment and create column 
    comment = self.commentDict[param]
    self.Column.param = param
    self.Column.valueArr = self.valueArr
    self.Column.comment = comment
  
  """  
  method to return current UT time (for FITS file generation)
  """               
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
        
  """  
  method to return ScanLog FITS file header
  """  
  def readScanLog_Header(self):
    hdu = fits.open(self.projectPath + '/ScanLog.fits')
    return hdu[0].header      
  
  """  
  method to construct a primary header for the SDFITS binary Table
  """  
  def constructPriHDUHeader(self):
    priHeader=fits.Header() ## blank primary header
    currentUTC = self.getCurrentUTC() ## get current UTC time to place in SDFITS header
    priHeader.set('DATE',currentUTC, 'date and time this HDU was created, UTC')
    scanLogHeader = self.readScanLog_Header() ## Collect header metadata from ScanLog
    origin = scanLogHeader['ORIGIN'] ## GBT
    priHeader.set('ORIGIN',origin,'origin of observation')
    telescope = scanLogHeader['TELESCOP'] ## GBT
    priHeader.set('TELESCOP',telescope,'the telescope used')
    
    ## Set other metadata
    priHeader.set('INSTRUME','FLAGBF','backend')
    priHeader.set('SDFITVER','sdfits-bf','SDFITS format for BF')
    priHeader.set('FITSVER','fits-bf','FITS format for BF')
    prihdu = fits.PrimaryHDU(header=priHeader)
    return prihdu 

  """  
  method for progress bar that is shown in the terminal
  """  
  def progressBar(self,value, endvalue, beam, bar_length=20):
    beamName = np.str(beam)
    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write("\rPercent of FITS table for  beam "+ beamName +": [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()     

  """  
  method to return the project ID as a string
  """  
  def getProjId(self):
      hdu = fits.open(self.projectPath + '/ScanLog.fits')
      return hdu[0].header['PROJID']
  """  
  method to convert Azimuth to Right Acension (geoapparent)
  """ 
  def az2ra(self,LST, Az, El, dec):
    ## convert all angles to radians
    azRad = np.deg2rad(Az) ## Azimuth of NCP is 0 for GBT. This code assumes 180. 
    posInds = np.where(Az > 360)
    negInds = np.where(Az < 180) ## need to rotate multiply by negative 1 while in 
    elRad = np.deg2rad(El)
    decRad = np.deg2rad(dec)
    ha = (-1) * np.arccos(1/np.cos(decRad)*(np.sin(elRad)*np.cos(self.GBTLAT) - np.cos(elRad)*np.cos(azRad)*np.sin(self.GBTLAT)))
    ha[posInds] = ha[posInds] * (-1)
    ha[negInds] = ha[negInds] * (-1)
    return LST + np.rad2deg(ha)
  
  """  
  method to convert Elevation to Declination (geoapparent)
  """ 
  def el2dec(self, LST, Az, El):
    ## convert all angles to radians
    azRad = np.deg2rad(360 - Az) ## Azimuth of NCP is 0 for GBT. This code assumes 180. 
    elRad = np.deg2rad(El)
    return np.rad2deg(np.arcsin(np.sin(elRad)*np.sin(self.GBTLAT) + np.cos(elRad)*np.cos(azRad)*np.cos(self.GBTLAT)))

  """
  method to transform the beam offsets to Ra/Dec coordinates. The current state of the binary table HDU is passed in. 
  That same binary table HDU with the updated coordinates is returned. 
  """
  def offsetCorrection(self, hdu):
    print('\n')
    print('Applying beam offsets...')
    ## extend lists that hold Antenna FITS files by 2x to describe both XX and YY pol
    extCoordSys = [self.coordSysList[0]] * 2 * len(self.coordSysList)
    extCoordSys[0::2] = self.coordSysList
    extCoordSys[1::2] = self.coordSysList
    extOffCoordSys = [self.offCoordSysList[0]] * 2 * len(self.offCoordSysList)
    extOffCoordSys[0::2] = self.offCoordSysList
    extOffCoordSys[1::2] = self.offCoordSysList
    extObscAzArr = np.zeros(len(self.refractList)*2)
    extObscAzArr[0::2] = self.obscAzList
    extObscAzArr[1::2] = self.obscAzList
    extObscElArr = np.zeros(len(self.refractList)*2)
    extObscElArr[0::2] = self.obscElList
    extObscElArr[1::2] = self.obscElList
    extSMntAzArr = np.zeros(len(self.refractList)*2)
    extSMntAzArr[0::2] = self.smntAzList
    extSMntAzArr[1::2] = self.smntAzList   
    extSMntElArr = np.zeros(len(self.refractList)*2)
    extSMntElArr[0::2] = self.smntElList
    extSMntElArr[1::2] = self.smntElList
    extSObscAzArr = np.zeros(len(self.refractList)*2)
    extSObscAzArr[0::2] = self.sobscAzList
    extSObscAzArr[1::2] = self.sobscAzList
    extSObscElArr = np.zeros(len(self.refractList)*2)
    extSObscElArr[0::2] = self.sobscElList
    extSObscElArr[1::2] = self.sobscElList
    extRefractArr = np.zeros(len(self.refractList)*2)
    extRefractArr[0::2] = self.refractList
    extRefractArr[1::2] = self.refractList 
    extRaArr = np.zeros(len(self.refractList)*2)
    extRaArr[0::2] = self.raList
    extRaArr[1::2] = self.raList
    extDecArr = np.zeros(len(self.refractList)*2)
    extDecArr[0::2] = self.decList
    extDecArr[1::2] = self.decList
    extMntAzArr = np.zeros(len(self.refractList)*2)
    extMntAzArr[0::2] = self.mntAzList
    extMntAzArr[1::2] = self.mntAzList
    extMntElArr = np.zeros(len(self.refractList)*2)
    extMntElArr[0::2] = self.mntElList
    extMntElArr[1::2] = self.mntElList    

    ## we have data DMJD, Maj, Min, Az, El, Ra, Dec, velocity, temp, pressure, & humidity for a single pol 
    ## extend by 2x to describe both XX and YY pol
    extCoordSys.extend(extCoordSys)
    newMajArr = np.zeros(len(self.refractList)*2)
    newMinArr = np.zeros(len(self.refractList)*2)
    newAzArr = np.zeros(len(self.refractList)*2)
    newElArr = np.zeros(len(self.refractList)*2)
    newRaArr = np.zeros(len(self.refractList)*2)
    newDecArr = np.zeros(len(self.refractList)*2)
    extLstArr = np.zeros(len(self.refractList)*2)
    dataDMJDArr = np.zeros(len(self.refractList)*2)
    extvelArr = np.zeros(len(self.refractList)*2)
    dataDMJDArr = np.zeros(len(self.dataDMJDList)*2)
    dataDMJDArr[0::2] = self.dataDMJDList
    dataDMJDArr[1::2] = self.dataDMJDList 

    beamElArr = np.deg2rad(hdu.data['FEEDEOFF'])
    beamXElArr = np.deg2rad(hdu.data['FEEDXOFF'])

    extFreqArr = hdu.data['CRVAL1']
    velArr = hdu.data['VELOCITY'] / 1000.0 ## km/s

    tempArr = hdu.data['TAMBIENT']
    pressArr = hdu.data['PRESSURE']
    pressArr*=1/0.75006375541921 ## mb
    humArr = hdu.data['HUMIDITY']

    extLstArr = hdu.data['LST']

    lstStartRads = np.deg2rad(extLstArr/3600.0*15) ## radians

    ## loop through each coordinate to precess from J2000 (mean place) to geoapparent RA/Dec
    for coordIdx in range(0, len(extRaArr)):
      raVal = extRaArr[coordIdx]
      decVal = extDecArr[coordIdx]
      velVal = velArr[coordIdx]
      dmjdVal = dataDMJDArr[coordIdx]
      deltaUT = self.utCorrList[coordIdx]
      polXVal = self.polMotionXList[coordIdx]
      polYVal = self.polMotionYList[coordIdx]
      tempVal = tempArr[coordIdx]
      pressVal = pressArr[coordIdx]
      humVal = humArr[coordIdx]
      waveVal = self.c / extFreqArr[coordIdx]/1e6 ## micrometers
      lapseVal = 0.0065

      ## convert J2000 -> geoapparent (center of Earth)
      geoRa, geoDec = pysla.slalib.sla_map(np.deg2rad(raVal), np.deg2rad(decVal), 0.0, 0.0, 0.0, velVal, 2000.0, dmjdVal)

      ## convert from center-of-earth to observed at GBO
      obsAz, obsZen, obsHA, obsDec, obsRa = pysla.slalib.sla_aop(geoRa, geoDec, dmjdVal, deltaUT, np.deg2rad(self.GBTLONG), np.deg2rad(self.GBTLAT), self.GBTHGT, polXVal, polYVal, tempVal, pressVal, humVal, 299792458.0/1420.405752e6*1e6, 0.0065)
      
      ## convert from obs eq to horizontal (use local definitation of hour angle)
      #obsAz, obsEl = pysla.slalib.sla_e2h(lstStartRads[coordIdx] - geoRa, geoDec, np.deg2rad(self.GBTLAT))
      obsAz, obsEl = pysla.slalib.sla_e2h(lstStartRads[coordIdx] - obsRa, obsDec, np.deg2rad(self.GBTLAT))

      ## apply refraction correction and beam offsets
      obsEl_refract = obsEl + np.deg2rad(extRefractArr[coordIdx])
      newElVal = obsEl_refract - beamElArr[coordIdx] 
      newAzVal = obsAz - (beamXElArr[coordIdx] / np.cos(newElVal))

      ## place in arrays
      newAzArr[coordIdx] = np.rad2deg(newAzVal)
      newElArr[coordIdx] = np.rad2deg(newElVal)

      ## put the major/minor array in the correct coordinate system based on INDICSYS
      if extCoordSys[coordIdx] == 'AZEL':
        # major and minor are az and el
        newMajArr[coordIdx] = np.rad2deg(newAzVal)
        newMinArr[coordIdx] = np.rad2deg(newElVal)
      ## convert to Cross-Elevation & Elevation by subtracting off pointing model if in engineering coordinates
      elif extOffCoordSys[coordIdx] == 'AZEL_ENCODER':
        newMajArr[coordIdx] = ((extMntAzArr[coordIdx] - extObscAzArr[coordIdx]) - (extSMntAzArr[coordIdx] - extSObscAzArr[coordIdx]))*np.cos(newElVal)
        newMinArr[coordIdx] = (extMntElArr[coordIdx] - extObscElArr[coordIdx]) - (extSMntElArr[coordIdx] - extSObscElArr[coordIdx])
      else:
        ## convert from obs horiz to eq 
        newObsHA, newObsDec = pysla.slalib.sla_h2e(newAzVal, newElVal, np.deg2rad(self.GBTLAT))
        newObsRa = lstStartRads[coordIdx] - newObsHA

        ## convert from observed at GBO to geocentric (center of Earth)
        newGeoRa, newGeoDec = pysla.slalib.sla_oap('R', newObsRa, newObsDec, dmjdVal, deltaUT, np.deg2rad(self.GBTLONG), np.deg2rad(self.GBTLAT), self.GBTHGT, polXVal, polYVal, tempVal, pressVal, humVal, 299792458.0/1420.405752e6*1e6, 0.0065)
        ## convert from geocenteric to mean
        newRa, newDec = pysla.slalib.sla_amp(newGeoRa, newGeoDec, dmjdVal, 2000.0)

        ## convert radians to degrees
        newRaArr[coordIdx] = np.rad2deg(newRa)
        newDecArr[coordIdx] = np.rad2deg(newDec)
        newMajArr[coordIdx] = np.rad2deg(newRa)
        newMinArr[coordIdx] = np.rad2deg(newDec)

        if extCoordSys[coordIdx] == 'GALACTIC':
          glon, glat = pysla.slalib.sla_eqgal(newRa, newDec)
          newMajArr[coordIdx] = np.rad2deg(glon)
          newMinArr[coordIdx] = np.rad2deg(glat)

    ## update ra/dec array to be used for radial velocity correction
    self.newRaArr = newRaArr
    self.newDecArr = newDecArr

    ## update header values
    hdu.data['CRVAL2'] = newMajArr
    hdu.data['CRVAL3'] = newMinArr
    hdu.data['TRGTLONG'] = newMajArr
    hdu.data['TRGTLAT'] = newMinArr
    hdu.data['AZIMUTH'] = newAzArr
    hdu.data['ELEVATIO'] = newElArr

    ## return SDFITS binary table with updated values
    return hdu
  """
  method to calculate doppler correction and put frequency in OPTI-HEL vel def/ref frame 
  """
  def radVelCorrection(self, hdu):
    print('\n')
    print('Applying Doppler Correction to HELIOCENTRIC...')

    ## collect all relevant data for doppler calculation
    utDate_Time = hdu.data['DATE-OBS']
    cenFreqsArr = hdu.data['CRVAL1']
    restFreqArr = hdu.data['RESTFREQ']
    for velIter in range(0,len(self.newRaArr)): ## loop through to calculate each integration/polarization's correciton 
      raVal = self.newRaArr[velIter] ## RA
      decVal = self.newDecArr[velIter] ## Dec 
      utDateTimeVal = utDate_Time[velIter] ## UT Time
      cenFreqVal = cenFreqsArr[velIter] ## center frequency value
      restFreqVal = restFreqArr[velIter] ## restfreq

      ## put time in necessary format
      t = Time(utDateTimeVal, format='isot')
      utdate = t.iso[0:10]
      uttime = t.iso[11:]
      radVelCorr_HEL, radVelCorr_LSR = self.radvelcorrObj.correctVel(utdate,uttime,raVal,decVal) ## calculate correction
      ## compute optical velocity of ref freq

      vOpt = self.c*(restFreqVal/cenFreqVal - 1)
      ## add radial correction
      newVOpt = vOpt + radVelCorr_HEL
      ## now convert back to frequency
      newCenFreq = restFreqVal/(1+newVOpt/self.c)
      ## update reference frequency value
      cenFreqsArr[velIter] = newCenFreq
    ## update columns
    hdu.data['CRVAL1'] = cenFreqsArr
    hdu.data['OBSFREQ'] = cenFreqsArr
    hdu.data['VELDEF'] = 'OPTI-HEL'
    return hdu

  """
  method that drives the creation of a list of FITS columns to be placed in a single binary table. 
  The main task is to loop over SDFITS keywords and collate the associated metadata with calls to various
  internal methods. Spatial coordinate offsets and Doppler corrections are calculated/applied before making final
  edits to the SDITS binary table header. The full HDU is returned to PAF_Filler.py to be written out to disk.
  """
  def constuctBinTableHeader(self):    
        
    ## if first iteration, construct primary HDU
    if self.firstIter:
      prihdu = self.constructPriHDUHeader()   
    binHeader = fits.Header()
    
    ## load list of required SDFITS keywords
    keywordList = np.loadtxt(keywordPath,dtype='str')
    keyWordArr = keywordList.astype(str)
    
    ## define empty lists to hold comments and params
    commentList = []
    paramList = []
    
    ## loop through each keyword
    for keyIdx in range(0,len(keyWordArr),3):
      keyword = keyWordArr[keyIdx]         
      
      ## utilize the keyword dictionary to get the associated parameter
      param = self.keyToParamDict[keyword]
      """
      knowing the required parameter, use a dictionary of associated functions
      which provide required metadata information and place data
      """
      self.funcDict[param.strip()](param.strip())  
      formKey = keyWordArr[keyIdx+1] ## gets the required TFORM and TUNIT keys from the SDFITS keyword array
      unitKey = keyWordArr[keyIdx+2]
      form = self.keyToParamDict[formKey]
      if formKey == 'TFORM7':
          form = np.str(len(self.Column.valueArr[0,:]))+'E'
      unit = self.keyToParamDict[unitKey]
      col = fits.Column(name=self.Column.param, format=form, array=self.Column.valueArr, unit=unit) ## construct the FITS column with returned array of values
      if keyIdx == 0: ## if first iteration, create the object FITS column attribute. This will be added to each iteration
        self.cols = fits.ColDefs([col])
      else:
        self.cols.add_col(col)
      commentList.append(self.Column.comment) ## append comment to list
      paramList.append(self.Column.param) ## append parameter to list

      ## update user
      self.progressBar(keyIdx, len(keyWordArr),self.beamNum)

      ## reset self.Column object
      self.Column = collections.namedtuple('Column',['param','valueArr','comment'])

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
    rowBytes = tblHdu.size/(2 * totalInts) ## calculate row bytes based on the total integrations and polarizations divided by total table size
    tblHdu.header.set('NAXIS1',int(rowBytes),'width of table in bytes')
    tblHdu.header.set('NAXIS2',totalInts * 2,'number of rows in table')
    tblHdu.header.set('PCOUNT',0,'size of special data area')
    tblHdu.header.set('GCOUNT',1,'one data group (required keyword)')
    tblHdu.header.set('TFIELDS',70,'number of fields in each row')
    tblHdu.header.set('EXTNAME','SINGLE DISH', 'name of this binary table extension')       
    tblHdu.header.set('XTENSION', 'BINTABLE', 'binary table extension')
    idx = 0        
    for i in range(0,len(keyWordArr),3): ## loop through to set keywords in header. A new keyword is every third element in list
        keyword = keyWordArr[i]
        tblHdu.header.set(keyword, paramList[idx], commentList[idx]) 
        idx+=1
    ## insert final header comments
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
    
    ## finally, create HDU List with primary header if first iteration, else just a HDU list with generated table
    if self.firstIter == True:
      thduList = fits.HDUList([prihdu, tblHdu])
    else: 
      thduList = tblHdu

    ## return finalized HDU to PAF_Filler.py
    return thduList
    
