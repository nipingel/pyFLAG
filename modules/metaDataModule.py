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
    stuff = {'XTENSION':'BINTABLE'
              'BITPIX':8
              'NAXIS':2
              'NAXIS1':## TODO: width
              'NAXIS2':## TODO: rows
              'PCOUNT':0
              'GCOUNT':1
              'TFIELDS':##TODO: num of fields
              'TTYPE1':'OBJECT'
              'TFORM1':'32A'
              'TUNIT1':''
              'TELESCOP':'NRAO_GBT'
              'TTYPE2':
              'TFORM2':
              'TUNIT2':
              'TTYPE3':
              'TFORM3':
              'TUNIT3':
              'TTYPE4':
              'TFORM4':
              'TUNIT4':
              'TTYPE5':
              'TFORM5':
              'TUNIT5':
              'TTYPE6':
              'TFORM6':
              'TUNIT6':
              'COMMENT':
              'COMMENT':
              'TTYPE7':
              'TFORM7':
              'TUNIT7':
              'TTYPE8':
              'TFORM8':
              'TUNIT8':
              'TTYPE9':
              'TFORM9':
              'TUNIT9':
              'TTYPE10':
              'TFORM10':
              'TUNIT10':
              'TTYPE11':
              'TFORM11':
              'TUNIT11':
              'TTYPE12':
              'TFORM12':
              'TUNIT12':
              'TTYPE13':
              'TFORM13':
              'TUNIT13':
              'TTYPE14':
              'TFORM14':
              'TUNIT14':
              'TTYPE15':
              'TFORM15':
              'TUNIT15':
              'TTYPE16':
              'TFORM16':
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
