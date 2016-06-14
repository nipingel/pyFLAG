# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 10:44:55 2016
Module containing methods/functions for beamforming. Raw data will be correlations (in FISHFITS order). 
Since the raw correlations are initially in a long 1D array, they must first be rearranged in matrices. 
The beamforming weights are then read in and applied to these covariance matrices to create spectra for 
each individual beam 
@author: npingel
"""
##imports
from astropy.io import fits
import numpy as np
import os
import glob

class BeamformingModule:
    
    def __init__ (self):
        return
    
    def getRawCorrelations(self,fitsName):
        hdu = fits.open(fitsName)
        corrData = hdu[1].data.field('DATA')
        return corrData
    ##TODO:
    ##What is the format of the weight files?        
    def getWeights():
        return    
    
    def processPol(self,corrMatrix,pol):
        R = np.copy(corrMatrix[0:19, 0:19]) ##drop empty data stream
        Reval, Revec = np.linalg.eig(R)
        
        aomega = Revec[:,Reval==max(Reval)] ##TODO: should this be here?
        
        ##TODO:
        ##where do we read weights from?
        w = np.ones([19])
        spectrum = np.dot(w.conj().T,np.dot(R,w))
        return spectrum
        
    
    def getCorrelationCube(self,dataVector):
            num_chan = 40
            # The data order in the 'cross_accum' buffer is ch1Xch1_freq1_real,
            # ch1Xch1_freq1_imag, ch1Xch2_freq1_real... ch1Xch20_freq1_imag,
            # ch2Xch1_freq1_real...ch2Xch20_freq1_imag, ch3Xch1_freq1_real...
            # ch1Xch1_freq2_real, ch1Xch1_freq2_imag...ch20Xch20_freq2_imag,
            # ch1Xch1_freq3_real...ch20Xch20_freqN_imag
          
            num_cor = 820.
            num_freq = len(dataVector) / num_cor
            ret_dat = np.zeros([num_freq, num_chan, num_chan], dtype=np.complex64)
            offset = 0
            for row in range(num_chan) :
                for col in range(row, num_chan) :
                    ret_dat[:,row,col] = dataVector[offset::num_cor]
                    if row != col :
                        ret_dat[:,col,row] = dataVector[offset::num_cor].conj()
                    offset += 1   
            return ret_dat
    
    ## Put 1D correlation array into 20x20 matrix (element 20,20 irrelevant)
    def unpackCorrelations(self,dat):
        xpol = [0,1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33,36,37]
        xdip = [1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20]
     
        ypol = [2,3,6,7,10,11,14,15,18,19,22,23,26,27,30,31,34,35,38,39]
        ydip = [1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20]
    
        xdat = np.empty_like(dat[0:20,0:20])
        ydat = np.empty_like(dat[0:20,0:20])
    
        for cnt1 in range(20) :
          for cnt2 in range(20) :
             dipind1 = xdip[cnt1] - 1
             dipind2 = xdip[cnt2] - 1
             polind1 = xpol[cnt1]
             polind2 = xpol[cnt2]
             xdat[dipind1, dipind2] = dat[polind1, polind2]
             xdat[dipind2, dipind1] = dat[polind2, polind1]
    
             dipind1 = ydip[cnt1] - 1
             dipind2 = ydip[cnt2] - 1
             polind1 = ypol[cnt1]
             polind2 = ypol[cnt2]
             ydat[dipind1, dipind2] = dat[polind1, polind2]
             ydat[dipind2, dipind1] = dat[polind2, polind1]
    
        return xdat, ydat 
    
        
    def getSpectralArray(self,fitsName):
        dataArr = self.getRawCorrelations(fitsName)
        for ints in range(0,len(dataArr[:,0])):     
            dataVector=dataArr[ints,:]
            corrCube = self.getCorrelationCube(dataVector)
            corrShape = corrCube.shape
            num_freq = corrShape[0]
            spectrumArr_X = np.zeros([num_freq], dtype='float32')    
            spectrumArr_Y = np.zeros([num_freq], dtype='float32')  
            for z in range(0,num_freq):
                dat = corrCube[z,:,:]
                xdat, ydat = self.unpackCorrelations(dat)
                spectrumArr_X[z] = self.processPol(xdat,0) ##0 is XX
                spectrumArr_Y[z] = self.processPol(ydat,1) ##1 is YY
        return spectrumArr_Sort_X,spectrumArr_Sort_Y
        
