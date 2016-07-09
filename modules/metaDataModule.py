# -*- coding: utf-8 -*-
"""
This module (i.e. class) will collate the required metadata to construc the primary HDU and SDFITS binary table. 

@author: npingel
"""
from astropy.io import fits
import glob
import os 
import numpy as np
import datetime
import collections


class MetaDataModule:
    
        

    def __init__(self,fitsList,numInts):
        self.fitsList = fitsList
        self.numInts = numInts
        self.Column = collections.namedtuple('Column',['param','valueArr','comment'])
        self.cols = None

        def getSMKey(self,param): ##TODO: get actual shared memory values
            comment = self.commentDict[param]
            value = 0.
            self.Column.param = param
            valueArr = np.full([self.numInts*len(self.fitsList)],value)         
            self.Column.valueArr = valueArr
            """
            for file in range(0,len(self.fitsList)):            
                HDU = fits.open(fitsList[file])
                value = goHDU[0].header[param]
                idx = 0
                for i in range(0,self.numInts):
                    valueArr[(idx*self.numInts)+i] = value
            """
            self.Column.comment = comment
            
        ##TODO: error handling (e.g. keyword not found)
        def getGOFITSParam(self,param):
            os.chdir('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/GO/') ##TODO: run on flag03 
            valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='object')           
            idx = 0  
            if param == 'TRGTLONG':
                param = 'RA'
            if param == 'TRGTLAT':
                param = 'DEC'
            for file in range(0,len(self.fitsList)):            
                goHDU = fits.open(fitsList[file])
                value = goHDU[0].header[param]
                for i in range(0,self.numInts):
                    valueArr[(idx*self.numInts)+i] = value
                idx+=1
                    
            comment = self.commentDict[param]
            self.Column.param = param
            self.Column.valueArr = valueArr
            self.Column.comment = comment
        
        def getAntFITSParam(self,param):
            os.chdir('/Users/npingel/Desktop/Research/FLAG/pros/exampleData/Antenna/') ##TODO: run on flag03 
            valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='object')           
            idx = 0  
            for file in range(0,len(self.fitsList)):            
                goHDU = fits.open(fitsList[file])
                if param == 'PRESSURE':
                    param = 'AMBPRESS'
                elif param == 'TAMBIENT':
                    param = 'AMBTEMP'
                elif param == 'HUMIDITY':
                    param = 'AMBHUMID'
                value = antHDU[0].header[param]
                if param == 'PRESSURE':
                    value=value*0.75006375541921                 
                elif param == 'TAMBIENT':
                    value+=273.0
                for i in range(0,self.numInts):
                    valueArr[(idx*self.numInts)+i] = value
                idx+=1
                    
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
       
       def getTsys(self,param):
            valueArr = np.empty([self.numInts*len(self.fitsList)],dtype='object')           
            idx = 0
            for file in range(0,len(self.fitsList)):            
                for i in range(0,self.numInts):
                    valueArr[(idx*self.numInts)+i] = 1.0
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
                            'HUMIDITY':'relative humidity'
                            
              }
        self.funcDict = {'OBJECT':getGOFITSParam,
                         'BANDWID':getSMKey,
                         'DATE-OBS':getSMKey,
                         'DURATION':getSMKey,
                         'EXPOSURE':getSMKey,
                         'TSYS':getTsys,
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
                         }
        ##Parameter Dictionary
        self.keyToParamDict = {'XTENSION':'BINTABLE',
              'TTYPE1':'OBJECT',
              'TFORM1':'32A',
              'TUNIT1':'',
              'TTYPE2': 'BANDWID', ##TODO: get from shared memory
              'TFORM2':'1D',
              'TUNIT2': 'Hz',
              'TTYPE3': 'DATE-OBS', ##TODO: get from shared memory (MJD)
              'TFORM3': '22A',
              'TUNIT3': '',
              'TTYPE4': 'DURATION', ##TODO: get from shared memory
              'TFORM4':'1D',
              'TUNIT4': 's',
              'TTYPE5': 'EXPOSURE', ##TODO: get from shared memory
              'TFORM5': '1D',
              'TUNIT5': 's', 
              'TTYPE6': 'TSYS', 
              'TFORM6': '1D',
              'TUNIT6': 'K',
              'TTYPE7': 'DATA', ##TODO: comment
              'TFORM7':'8192E', ##TODO: comment;define this
              'TUNIT7':'', 
              'TTYPE8':'', ##TODO: comment and find dimesions of data array
              'TFORM8': '16A', 
              'TUNIT8': '', 
              'TTYPE9':'TUNIT7', ##TODO: 
              'TFORM9':'6A',
              'TUNIT9':'',
              'TTYPE10':'CTYPE1', ##TODO: comment
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
              'TTYPE14':'CTYPE2', ##TODO: comment
              'TFORM14': '4A', 
              'TUNIT14': '',
              'TTYPE15':'CRVAL2',
              'TFORM15': '1D',
              'TUNIT15': 'deg',
              'TTYPE16': 'CTYPE3', ##TODO: comment
              'TFORM16':'4A',
              'TUNIT16':'',
              'TTYPE17':'CRVAL3', 
              'TFORM17':'1D',
              'TUNIT17':'deg',
              'CTYPE4':'STOKES', ##TODO: comment
              'TTYPE18':'CRVAL4', 
              'TFORM18':'1I',
              'TUNIT18':'',
              'TTYPE19':'OBSERVER',
              'TFORM19':'32A',
              'TUNIT19':'',
              'TTYPE20':'OBSID',
              'TFORM20':'32A', 
              'TUNIT20':'',
              'PROJID':'', ## TODO: comment; find this
              'TTYPE21':'SCAN',
              'TFORM21':'1J',
              'TUNIT21':'',
              'TTYPE22':'OBSMODE', ##TODO comment; depends if there is noise injection?
              'TFORM22':'32A', 
              'TUNIT22':'',
              'TTYPE23':'FRONTEND',
              'TFORM23':'16A',
              'TUNIT23':'',
              'BACKEND':'', ##TODO comment; define this 
              'TTYPE24':'TCAL', ##TODO comment; find this
              'TFORM24':'1E',
              'TUNIT24':'K', 
              'TTYPE25':'VELDEF', ##TODO comment; find this 
              'TFORM25':'8A', 
              'TUNIT25':'', 
              'TTYPE26':'VFRAME', ## TODO coment; find this
              'TFORM26':'1D',
              'TUNIT26':'m/s',
              'TTYPE27':'RVSYS', ## TODO comment; find this 
              'TFORM27':'1D',
              'TUNIT27':'m/s', 
              'TTYPE28':'OBSFREQ', ##TODO comment; find this
              'TFORM28':'1D',
              'TUNIT28':'Hz',
              'TTYPE29':'LST', ##TODO comment; find this
              'TFORM29':'1D',
              'TUNIT29':'s',
              'TTYPE30':'AZIMUTH', ##TODO comment; find this
              'TFORM30':'1D',
              'TUNIT30':'deg',
              'TTYPE31':'ELEVATIO', ##TODO comment; find this
              'TFORM31':'1D',
              'TUNIT31':'deg',
              'TTYPE32':'TAMBIENT',
              'TFORM32':'1D',
              'TUNIT32':'K', 
              'TTYPE33':'PRESSURE', ## TODO comment; find this
              'TFORM33':'1D',
              'TUNIT33':'mmHg', 
              'TTYPE34':'HUMIDITY', ## TODO comment; find this
              'TFORM34':'1D',
              'TUNIT34':'',
              'SITELONG':'', ## find and define this
              'SITELAT':'',## find and define this
              'SITEELEV':'',## find and define this
              'TTYPE35':'RESTFRQ',
              'TFORM35':'1D',
              'TUNIT35':'Hz',
              'TTYPE36':'FREQRES', ##TODO comment; find this
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
              'TTYPE41':'SAMPLER', ##TODO comment
              'TFORM41':'12A',
              'TUNIT41':'FEED', ##TODO comment  
              'TTYPE42':'1I',
              'TFORM42':'',
              'TUNIT42':'',
              'TTYPE43':'SRFEED', ##TODO:
              'TFORM43':'1I',
              'TUNIT43':'',
              'TTYPE44':'FEEDXOFF', ##TODO:
              'TFORM44':'1D',
              'TUNIT44':'deg',
              'TTYPE45':'FEEDEOFF', ##TODO:
              'TFORM45': '1D',
              'TUNIT45':'deg',
              'TTYPE46':'SUBREF_STATE', ##TODO:
              'TFORM46':'1I',
              'TUNIT46':'',
              'TTYPE47':'SIDEBAND', ##TODO:
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
              'TTYPE54':'TIMESTAMP', ##TODO:
              'TFORM54':'22A',
              'TUNIT54':'UTC', 
              'TTYPE55':'QD_XEL', ##TODO:
              'TFORM55':'1D',
              'TUNIT55':'deg',
              'TTYPE56':'QD_EL', ##TODO:
              'TFORM56':'1D',
              'TUNIT56':'deg',
              'TTYPE57':'QD_BAD', ##TODO:
              'TFORM57':'1I',
              'TUNIT57':'',
              'TTYPE58':'QD_METHOD', ##TODO:
              'TFORM58':'1A', 
              'TUNIT58':'',
              'TTYPE59':'VELOCITY',
              'TFORM59':'1D', 
              'TUNIT59':'m/s',
              'TTYPE60':'ZEROCHAN', ##TODO:
              'TFORM60':'1E',
              'TUNIT60':'',
              'TTYPE61':'DOPFREQ', ## TODO:
              'TFORM61':'1D', 
              'TUNIT61':'Hz', 
              'TTYPE62':'SIG', ##TODO:
              'TFORM62':'1A', 
              'TUNIT62':'',
              'TTYPE63':'CAL', ##TODO;
              'TFORM63':'1A',
              'TUNIT63':'',
              'TTYPE64':'CALTYPE', ##TODO:
              'TFORM64':'8A',
              'TTYPE65':'TWARM', ##TODO:
              'TFORM65':'1E', 
              'TUNIT65':'K',
              'TTYPE66':'TCOLD',
              'TFORM66':'1E',
              'TUNIT66':'K',
              'TTYPE67':'CALPOSITION', ##TODO: 
              'TFORM67':'16A' ,
              'TTYPE68':'IFNUM', ##TODO:
              'TFORM68':'1I',
              'TTYPE69':'PLNUM', ##TODO:
              'TFORM69':'1I',
              'TTYPE70':'FDNUM', ##TODO:
              'TFORM70':'1I',
              'EXTNAME':'SINGLE DISH', ##TODO:
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
                print(self.Column.valueArr)
                formKey = keyWordArr[keyIdx+1]
                unitKey = keyWordArr[keyIdx+2]
                form = self.keyToParamDict[formKey]
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
        binHeader['COMMENT'] = 'End of GBT-specific keywords/columns.'
        binHeader.set('EXTNAME','SINGLE DISH', 'name of this binary table extension')
        return binHeader
    def constructBinTableData(self):
        return
