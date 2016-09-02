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
import sys

class BeamformingModule:
    
    def __init__ (self):
        self.mapVector = np.loadtxt('/Users/npingel/Desktop/Research/FLAG/pros/SpectralFiller/misc/gpuToNativeMap.dat', dtype='int')

    def progressBar(self,value, endvalue, beam,xid,bar_length=20):

        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))

        sys.stdout.write("\rPercent of integrations filled in beam "+np.str(beam)+",xid "+str(xid)+": [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()    
        
    def getRawCorrelations(self,fitsName):
        hdu = fits.open(fitsName)
        corrData = hdu[1].data.field('DATA')
        return corrData      
    def getWeights(self,numChans,xid):
        hdu = fits.open('/Users/npingel/Desktop/Research/data/FLAG/TGBT16A_508/TGBT16A_508_03/Weights/2016_07_25_04:32:35_xid'+np.str(xid)+'_weights.fits')
        ##TODO: make above line more general        
        xWeights = np.zeros([numChans,40,7], dtype='complex64') 
        yWeights = np.zeros([numChans,40,7], dtype='complex64')  
        for i in range(0,7):
            polStr = 'X'
            data = hdu[1].data['Beam'+np.str(i)+polStr]
            totalElemInChan = 64*2
            for chan in range(0,numChans):
                weightPairs = data[(totalElemInChan*chan):(totalElemInChan*chan)+80]                 
                xWeights[chan,:,i].real = weightPairs[0:len(weightPairs):2]
                xWeights[chan,:,i].imag = weightPairs[1:len(weightPairs):2]
            polStr = 'Y'
            data = hdu[1].data['Beam'+np.str(i)+polStr]
            for chan in range(0,numChans):
                weightPairs = data[(totalElemInChan*chan):(totalElemInChan*chan)+80]                 
                yWeights[chan,:,i].real = weightPairs[0:len(weightPairs):2]
                yWeights[chan,:,i].imag = weightPairs[1:len(weightPairs):2]
        return xWeights,yWeights
    
    def processPol(self,corrMatrix,pol,w):
        R = corrMatrix
        spectrum = np.dot(w.conj().T,np.dot(R,w))
        return spectrum
        
    
    def getCorrelationCube(self,dataVector):
           ##Reorder to FISHIFITS
           newDataVector= np.zeros([820*25], dtype='complex64') ##make mode dependent 
           FITS_strt_idx = 0 
           for i in range(0,52800,2112): ##make mode dependent
               singleChanData = dataVector[i:i+2112]
               for z in range(0,len(self.mapVector)):
                   newDataVector[z+FITS_strt_idx] = singleChanData[self.mapVector[z]]
               FITS_strt_idx+=820   
            # The data order in the 'cross_accum' buffer is ch1Xch1_freq1_real,
            # ch1Xch1_freq1_imag, ch1Xch2_freq1_real... ch1Xch20_freq1_imag,
            # ch2Xch1_freq1_real...ch2Xch20_freq1_imag, ch3Xch1_freq1_real...
            # ch1Xch1_freq2_real, ch1Xch1_freq2_imag...ch20Xch20_freq2_imag,
            # ch1Xch1_freq3_real...ch20Xch20_freqN_imag
           num_chan = 40
           num_cor = 820
           num_freq = int(len(newDataVector) / num_cor)
           ret_dat = np.zeros([num_freq, num_chan, num_chan], dtype=np.complex64)
           offset = 0                
           for col in range(num_chan) :
               for row in range(col, num_chan) :
                   ret_dat[:,row,col] = newDataVector[offset::num_cor]
                   if row != col :
                       ret_dat[:,col,row] = newDataVector[offset::num_cor].conj()
                   offset += 1   
           return ret_dat
    
    ##TODO:Remove?    
    ## Put 1D correlation array into 20x20 matrix (element 20,20 irrelevant)
    def unpackCorrelations(self,dat):
        """        
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
        """
        xdat = np.empty_like(dat[0:20,0:20])
        ydat = np.empty_like(dat[0:20,0:20])
        xdat = dat[0:19]
        ydat = dat[20:39]
        return xdat, ydat 
        
    def getSpectralArray(self,fitsName,beam,xid):
        dataArr = self.getRawCorrelations(fitsName)
        dataVector=dataArr[0,:]
        corrCube = self.getCorrelationCube(dataVector)
        corrShape = corrCube.shape
        num_freq = corrShape[0] ##TODO: get from header?
        spectrumArr_X = np.zeros([len(dataArr[:,0]),num_freq], dtype='float32')    
        spectrumArr_Y = np.zeros([len(dataArr[:,0]),num_freq], dtype='float32') 
        xWeight,yWeight = self.getWeights(num_freq,xid)
        for ints in range(0,len(dataArr[:,0])):   
            self.progressBar(ints,len(dataArr[:,0]),beam,xid)            
            dataVector=dataArr[ints,:]
            corrCube = self.getCorrelationCube(dataVector)
            corrShape = corrCube.shape
            for z in range(0,num_freq):
                dat = corrCube[z,:,:]
                spectrumArr_X[ints,z] = np.real(self.processPol(dat,0,xWeight[z,:,beam])) ##0 is XX
                spectrumArr_Y[ints,z] = np.real(self.processPol(dat,1,yWeight[z,:,beam])) ##1 is YY
        return spectrumArr_X,spectrumArr_Y,ints+1 ##TODO: why return ints?
        
