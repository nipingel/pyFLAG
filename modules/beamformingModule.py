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

## class to read in BANK data from FITS; processes each integration to go from raw covariances to beam-formed spectra
class BeamformingModule:
    ## initialize function
    ## grabs the map from covariances to FISHFITS order and assigns path to data
    def __init__ (self, dataPath):
        self.mapVector = np.loadtxt('/users/npingel/FLAG/SpectralFiller/misc/gpuToNativeMap.dat', dtype='int')
        self.dataPath = dataPath
        ## get project ID
        dataPathSplit = dataPath.split('/')
        self.projectId = dataPathSplit[-2] ## always second last
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
    def getWeights(self,numChans): ## TODO: finish weights
        hdu = fits.open(dataPath+'/weight_files/' + 'w_' + self.projectID + '*.fits')
        ## weight array in the form of freqchans, element data, beams
        xWeights = np.zeros([numChans,40,7], dtype='complex64') 
        yWeights = np.zeros([numChans,40,7], dtype='complex64')  
        for i in range(0,7):
            ## process XX Pol
            polStr = 'X'
            ## open file
            data = hdu[1].data['Beam'+np.str(i)+polStr]
            ## 64 complex pairs per freq channel
            totalElemInChan = 64*2
            ## can drop last 24 correlation pairs as they are unused and zero. First 80 elements are good. 
            for chan in range(0,numChans):
                weightPairs = data[(totalElemInChan*chan):(totalElemInChan*chan)+80]                 
                xWeights[chan,:,i].real = weightPairs[0:len(weightPairs):2]
                xWeights[chan,:,i].imag = weightPairs[1:len(weightPairs):2]
            ## open and process YY Pol
            polStr = 'Y'
            data = hdu[1].data['Beam'+np.str(i)+polStr]
            for chan in range(0,numChans):
                weightPairs = data[(totalElemInChan*chan):(totalElemInChan*chan)+80]                 
                yWeights[chan,:,i].real = weightPairs[0:len(weightPairs):2]
                yWeights[chan,:,i].imag = weightPairs[1:len(weightPairs):2] 
        return xWeights,yWeights
    
    def processPol(self,R,pol,w):
        spectrum = np.dot(w.conj().T,np.dot(R,w))
        return spectrum
        
    ## function which reorders correlation bandpass to FISHFITS order;
    ## inputs is a covariance bandpass for a single integration and number of freq chans
    def getCorrelationCube(self,dataVector, numFreqs):
           ##Reorder to FISHFITS
           ## explicit order is below, but consists of 820 complex pairs* numFreqs
           numCorr = 820
           newDataVector= np.zeros([820*numFreqs], dtype='complex64') ##make mode dependent 
           ## set starting index for new data buffer
           FITS_strt_idx = 0
           ## put correlations into real/imag pairs
           for i in range(0,freqChans*2112, 2112):
               singleChanData = dataVector[i:i+2112]
               for z in range(0,len(self.mapVector)):
                   newDataVector[z+FITS_strt_idx] = singleChanData[self.mapVector[z]]
               ## update starting index
               FITS_strt_idx+=820   
           # The data order in the 'cross_accum' buffer is ch1Xch1_freq1_real,
           # ch1Xch1_freq1_imag, ch1Xch2_freq1_real... ch1Xch20_freq1_imag,
           # ch2Xch20_freq1_imag ...  ch3Xch3_freq1_real ... ch40xch40_freq1_imag
           # ch1Xch1_freq2_real ... 
           # ch1Xch1_freq3_real...ch20Xch20_freqN_imag
           ## Once sorted, we can construct a retCube with dim1 = numElements, dim2 = numElements, dim3 = freqChans
           numElem = 40
           retCube = np.zeros([numElem, numElem, numFreqs], dtype=np.complex64)
           offset = 0
           ## we wish to fill out the cube such that 0,0,0 (numElem, numElem, freqChan) is equal to ch1Xch1_freq1,
           ## Since the correlations are redundant, we can loop
           ## through columns while decreasing the number of rows looped through each column iteration
           ## by one. (filling in 0,0 ... 39,0; new column: 0,1 ... 38,1). Using the fact that 820 corr pairs 
           ## make up each freq channel, once on a specific row, col element, we select each 820 pair to fill 
           ## in the frequency axis. The 'offset', which denotes the correlation pair, must be increased by one. 
           ## The lower index of rows is increased to avoid placing a redundant correlation
           for col in range(numElem) :
               for row in range(col, numElem) :
                   retCube[row,col, :] = newDataVector[offset::num_cor]
                   if row != col :
                       retCube[col,row,:] = newDataVector[offset::num_cor].conj()
                   offset += 1   
           return retCube
    
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
        
    def getSpectralArray(self,fitsName,beam):
        ## open BANK FITS file
        dataArr = self.getRawCorrelations(fitsName)
        
        ## get number of freq channels
        num_freqs = dataArr.shape[1] / 2112 ## always 2112 complex pairs per frequency channel
        
        ## create bandpass containers. Rows are ints; colums represent integration bandpasses
        spectrumArr_X = np.zeros([len(dataArr[:,0]),num_freq], dtype='float32')    
        spectrumArr_Y = np.zeros([len(dataArr[:,0]),num_freq], dtype='float32') 
        
        ## unpack the weights from the associated binary file
        xWeight,yWeight = self.getWeights(num_freq,xid)
        ## loop through integrations and process covariance bandpass to beam-formed spectra
        for ints in range(0,len(dataArr[:,0])):   
            ## TODO: best place for this??
            self.progressBar(ints,len(dataArr[:,0]),beam)            
            ## grab a covariance bandpass for a single integration
            dataVector=dataArr[ints,:]
            ## send to get in FISHFITS order before sorting into 'cube' of shape (40,40,freqChans)
            corrCube = self.getCorrelationCube(dataVector)
            ## loop through frequency channels to apply weights. Inputs into processPol
            ## are a 40x40 covariance matrix, pol flag, and a length 40 weight vector
            for z in range(0,num_freq):
                dat = corrCube[z,:,:]
                spectrumArr_X[ints,z] = np.real(self.processPol(dat,0,xWeight[z,:,beam])) ##0 is XX
                spectrumArr_Y[ints,z] = np.real(self.processPol(dat,1,yWeight[z,:,beam])) ##1 is YY
        return spectrumArr_X,spectrumArr_Y,ints+1 ##TODO: why return ints?
        
