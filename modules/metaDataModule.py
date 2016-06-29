# -*- coding: utf-8 -*-
"""
This module (i.e. class) will collate the required metadata to construc the primary HDU and SDFITS binary table. 

@author: npingel
"""
from astropy.io import fits
import numpy as np
import datetime

class MetaDataModule:
    
    ## parameter dictionary
    stuff = {'XTENSION':'BINTABLE',
              'BITPIX':8,
              'NAXIS':2,
              'NAXIS1': '',## TODO: width
              'NAXIS2': '',## TODO: rows
              'PCOUNT':0,
              'GCOUNT':1,
              'TFIELDS':'',##TODO: num of fields
              'TTYPE1':'OBJECT', ##TODO: comment
              'TFORM1':'32A',
              'TUNIT1':'',
              'TELESCOP':'NRAO_GBT', ##TODO: comment
              'TTYPE2': 'BANDWID', ##TODO: coment
              'TFORM2':'1D',
              'TUNIT2': 'Hz',
              'TTYPE3': 'DATE-OBS', ##TODO: comment
              'TFORM3': '22A',
              'TUNIT3': '',
              'TTYPE4': 'DURATION', ##TODO: comment
              'TFORM4':'1D',
              'TUNIT4': 's',
              'TTYPE5': 'EXPOSURE', ##TODO: comment
              'TFORM5': '1D',
              'TUNIT5': 's', 
              'TTYPE6': 'TSYS', ##TODO: comment
              'TFORM6': '1D',
              'TUNIT6': 'K',
              'TTYPE7': 'DATA', ##TODO: comment
              'TFORM7':'8192E', ##TODO: comment;define this
              'TUNIT7':'', 
              'TTYPE8':'', ##TODO: comment and find dimesions of data array
              'TFORM8': '16A', 
              'TUNIT8': '', 
              'TTYPE9':'TUNIT7',
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
              'TTYPE19':'OBSERVER',##TODO: comment; find this
              'TFORM19':'32A',
              'TUNIT19':'',
              'TTYPE20':'OBSID', ##TODO comment; find this
              'TFORM20':'32A', 
              'TUNIT20':'',
              'PROJID':'', ## TODO: comment; find this
              'TTYPE21':'SCAN',##TODO: comment; find this
              'TFORM21':'1J',
              'TUNIT21':'',
              'TTYPE22':'OBSMODE', ##TODO comment; find this
              'TFORM22':'32A', 
              'TUNIT22':'',
              'TTYPE23':'FRONTEND', ##TODO comment; find this
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
              'TFORM26':'1D', '
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
              'TTYPE32':'TAMBIENT', ##TODO comment; find this
              'TFORM32':'1D,
              'TUNIT32':'K', 
              'TTYPE33':'PRESSURE', ## TODO comment; find this
              'TFORM33':'1D'
              'TUNIT33':'mmHg', 
              'TTYPE34':'HUMIDITY', ## TODO comment; find this
              'TFORM34':'1D',
              'TUNIT34':'',
              'SITELONG':'', ## find and define this
              'SITELAT':'',## find and define this
              'SITEELEV':'',## find and define this
              'TTYPE35':'RESTFREQ', ##TODO: comment; find this
              'TFORM35':'1D'
              'TUNIT35':'Hz'
              'TTYPE36':'FREQRES', ##TODO comment; find this
              'TFORM36':'1D',
              'TUNIT36':'Hz',
              'TTYPE38':'RADESYS', ##TODO comment; find this. 
              'TFORM38':'8A',
              'TTYPE39':'TRGTLONG', ##TODO comment
              'TFORM39':'1D',
              'TUNIT39':'deg', 
              'TTYPE40':'TRGTLAT', ##TODO comment
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
              'TTYPE48':'PROCSEQN', ##TODO:
              'TFORM48':'1I',
              'TUNIT48':''
              'TTYPE49':'PROCSIZE', ##TODO:
              'TFORM49':'1I',
              'TUNIT49':'',
              'TTYPE50':'PROCSCAN', ##TODO:
              'TFORM50':'16A',
              'TUNIT50':'',
              'TTYPE51':'PROCTYPE', ##TODO:
              'TFORM51':'16A',
              'TUNIT51':'',
              'TTYPE52':'LASTON', ##TODO:
              'TFORM52':'1J',
              'TUNIT52':'',
              'TTYPE53':'LASTOFF', ##TODO:
              'TFORM53':'1J'
              'TUNIT53':''
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
              'TUNIT57':''
              'TTYPE58':'QD_METHOD', ##TODO:
              'TFORM58':'1A', 
              'TUNIT58':'',
              'TTYPE59':'VELOCITY' ##TODO:
              'TFORM59':'1D', 
              'TUNIT59':'m/s'
              'TTYPE60':'ZEROCHAN', ##TODO:
              'TFORM60':'1E',
              'TUNIT60':''
              'TTYPE61':'DOPFREQ', ## TODO:
              'TFORM61':'1D', 
              'TUNIT61':'Hz';, 
              'TTYPE62':'SIG', ##TODO:
              'TFORM62':'1A', 
              'TUNIT62':'',
              'TTYPE63':'CAL', ##TODO;
              'TFORM63':'1A',
              'TUNIT63':'',
              'TTYPE64':'CALTYPE', ##TODO:
              'TFORM64':'8A'
              'TTYPE65':'TWARM', ##TODO:
              'TFORM65':'1E', 
              'TUNIT65':'K',
              'TTYPE66':'TCOLD',
              'TFORM66':'1E',
              'TUNIT66':'K',
              'TTYPE67':'CALPOSITION', ##TODO: 
              'TFORM67':'16A' 
              'TTYPE68':'IFNUM', ##TODO:
              'TFORM68':'1I'
              'TTYPE69':'PLNUM', ##TODO:
              'TFORM69':'1I'
              'TTYPE70':'FDNUM', ##TODO:
              'TFORM70':'1I'
              'EXTNAME':'SINGLE DISH', ##TODO:
              }

    def __init__(self):
        return
    
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
        
        for keyIdx in range(0,len(keyWordArr)):
            binHeader.set(keyWordArr[keyIdx],'BINTABLE', 'binary table extension')
            binHeader.set(keyWordArr[])
        binHeader.set('BITPIX','BINTABLE', 'binary table extension')
        binHeader.set('NAXIS1',2, '2-dimensional binary table')
        ##TODO:Descriptive keywords about table properties: NAXIS1, NAXIS2, PCOUNT, GCOUNT, TFIELDS
        binHeader['COMMENT'] = 'Start of SDFITS CORE keywords/columns.'
        ##TODO:SDFITS CORE KEYWORDS
        for keyIdx in range(0,len(keyWordArr)):
            binHeader.set(keyWordArr[keyIdx],'BINTABLE', 'binary table extension')
            binHeader.set(keyWordArr[])
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
