# -*- coding: utf-8 -*-
"""
This module (i.e. class) will collate the required metadata to construc the primary HDU and SDFITS binary table. 

@author: npingel
"""
from astropy.io import fits
import numpy as np
import datetime

class MetaDataModule:

## returns scanLog binary table
    def readScanLog_Data():
        scanLogHduList = fits.open('ScanLog.fits')
        return scanLogHduList[1].data
    
    ##returns scanLog header
    def readScanLog_Header():
        scanLogHduList = fits.open('ScanLog.fits')
        return scanLogHduList[0].header
    
    def getCurrentUTC():
        time = datetime.datetime.utcnow()
        dateStr=str(time.year)
        dateStr+='-'+str(time.strftime('%m'))
        dateStr+='-'+str(time.strftime('%d'))
        dateStr+='T'
        dateStr+=str(time.strftime('%H'))
        dateStr+=':'+str(time.strftime('%M'))
        dateStr+=':'+str(time.strftime('%S'))
        return dateStr
        
    def contstructPriHDUHeader():
        ## Collect header metadata from ScanLog
        priHeader=fits.Header()
        currentUTC = getCurrentUTC()
        priHeader.set('DATE',currentUTC, 'date and time this HDU was created, UTC')
        scanLogHeader = readScanLog_Header()
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

    def constuctBinTableHeader():    
        binHeader = fits.Header()
        binHeader.set('XTENSION','BINTABLE', 'binary table extension')
        binHeader.set('BITPIX','BINTABLE', 'binary table extension')
        binHeader.set('NAXIS1',2, '2-dimensional binary table')
        ##TODO:Descriptive keywords about table properties: NAXIS1, NAXIS2, PCOUNT, GCOUNT, TFIELDS
        binhdr['COMMENT'] = 'Start of SDFITS CORE keywords/columns.'
        ##TODO:SDFITS CORE KEYWORDS
        binhdr['COMMENT'] = 'End of SDFITS CORE keywords/columns.'
        binhdr['COMMENT'] = 'Start of SDFITS DATA column and descriptive axes.'
        ##TODO: SDFITS DATA KEYWORDS (including Beamformer specific)
        binhdr['COMMENT'] = 'End of SDFITS DATA column and descriptive axes.'
        binhdr['COMMENT'] = 'Start of SDFITS SHARED keywords/columns.'
        ##TODO: SDFITS SHARED KEYWORDS
        binhdr['COMMENT'] = 'End of SDFITS SHARED keywords/columns.'
        binhdr['COMMENT'] = 'Start of GBT-specific keywords/columns.'
        ##TODO: GBT-SPECIFIC KEYWORDS
        binhdr['COMMENT'] = 'Feed offsets ARE included in the CRVAL2 and CRVAL3 columns.'
        ##TODO: MORE GBT-SPECIFIC KEYWORDS
        binhdr['COMMENT'] = 'End of GBT-specific keywords/columns.'
        binHeader.set('EXTNAME','SINGLE DISH', 'name of this binary table extension')
        return binHeader
    def constructBinTable():
        
        
    
    
    
    
