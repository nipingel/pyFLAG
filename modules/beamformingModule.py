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
import glob
import matplotlib.pyplot as pyplot
## class to read in BANK data from FITS; processes each integration to go from raw covariances to beam-formed spectra
class BeamformingModule:
    ## initialize function
    ## grabs the map from covariances to FISHFITS order and assigns path to data
    def __init__ (self, dataPath):
        self.mapVector = np.loadtxt('/users/npingel/FLAG/SpectralFiller/misc/gpuToNativeMap.dat', dtype='int')
        self.dataPath = dataPath
        self.ints = 0
        ## get project ID
        dataPathSplit = dataPath.split('/')
        self.projectID = dataPathSplit[-3] ## always third last
                          
    def progressBar(self,value, endvalue, beam, xID,bar_length=20):

        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))

        sys.stdout.write("\rPercent of integrations filled in beam "+np.str(beam)+",xID "+str(xID)+": [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()    
        
    def getRawCorrelations(self,fitsName):
        hdu = fits.open(fitsName)
        corrData = hdu[1].data.field('DATA')
        return corrData      
    def getWeights(self,numChans, xID): ## TODO: finish weights
        weightFileList = glob.glob(self.dataPath+'weight_files/' + '*FullGrid.FITS')
        for wtFile in weightFileList:
            wtHDU = fits.open(wtFile)
            instID = wtHDU[0].header['XENGINE']
            if instID == xID:
                break
        print(wtFile)
        ## weight array in the form of freqchans, element data, beams
        xWeights = np.zeros([25,40,7], dtype='complex64') 
        yWeights = np.zeros([25,40,7], dtype='complex64')  
        for i in range(0,7):
            ## process XX Pol
            polStr = 'X'
            ## open file
            data = wtHDU[1].data['Beam'+np.str(i)+polStr]
            ## 64 complex pairs per freq channel
            totalElemInChan = 64*2
            ## can drop last 24 correlation pairs as they are unused and zero. First 80 elements are good. 
            for chan in range(0,25):
                weightPairs = data[(totalElemInChan*chan):(totalElemInChan*chan)+80]                 
                xWeights[chan,:,i].real = weightPairs[0::2]
                xWeights[chan,:,i].imag = weightPairs[1::2]
            ## open and process YY Pol
            polStr = 'Y'
            data = wtHDU[1].data['Beam'+np.str(i)+polStr]
            for chan in range(0,25):
                weightPairs = data[(totalElemInChan*chan):(totalElemInChan*chan)+80]                 
                yWeights[chan,:,i].real = weightPairs[0:len(weightPairs):2]
                yWeights[chan,:,i].imag = weightPairs[1:len(weightPairs):2] 
            #for chan in range(0,25): 
            #    for elem in range(1,41):
            #        print(elem-1)
            #        print(self.elemDict[elem])
            #        xWeights[chan, elem-1, i] = xWeights[chan, self.elemDict[elem], i]
            #        yWeights[chan, elem-1, i] = yWeights[chan, self.elemDict[elem], i]
        ## DEBUG
        if self.ints == 1:
            ## diagnositc plotting (weights)
            pyplot.figure()
            pyplot.plot(np.abs(xWeights[9,:,i]), label = 'XX-Pol', linewidth=2)
            pyplot.plot(np.abs(yWeights[9,:,i]), label = 'YY-Pol', linewidth=2)
            pyplot.xlabel('Data Channel')
            pyplot.ylabel('Magnitude')
            pyplot.title('Weight Vector; Coarse Channel: 130')
            pyplot.legend(loc=0)
            pyplot.savefig('/users/npingel/FLAG/2017Reduction/Plots/WeightVector_coarseChan130.pdf')
            pyplot.clf()
            pyplot.close()
            self.ints == 0
        ## DEBUG 
        return xWeights,yWeights
    
    def processPol(self,R,w):
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
           for i in range(0, numFreqs*2112, 2112):
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
                   retCube[row,col, :] = newDataVector[offset::numCorr]
                   if row != col :
                       retCube[col,row,:] = newDataVector[offset::numCorr].conj()
                   offset += 1   
           ##DEBUG
           if self.ints == 1:
               ## diagnostic plotting
               pyplot.figure()
               pyplot.title('Covariance Matrix; Fine Frequency Channel: 650')
               pyplot.xlabel('Data Channel')
               pyplot.ylabel('Data Channel')  
               pyplot.imshow(np.log10(abs(retCube[0:39,0:39,9])), extent = [0,39,0,39])
               pyplot.xlim(0,39)
               pyplot.ylim(0,39)
               pyplot.colorbar()
               pyplot.savefig('/users/npingel/FLAG/2017Reduction/Plots/CoVarMatrix_fineChan650_Int100.pdf')
               pyplot.clf()
               pyplot.close()
               self.ints = 0
           ##DEBUG
           
           ## if in PFB mode (i.e. 160 fine channels), we need to re-stitch the frequency channels as they are
           ## output in the wrong order from the pfb
           if numFreqs == 160:
               ## make 1D vector containing current indices for single 160 channel chunk
               origIdxArr = np.linspace(0, 159, 160, dtype='int32')
               ## reshape into a 32 (rows) x 5 (cols) array
               reshapeIdxArr = np.reshape(origIdxArr, (32,5))
               ## reshape the transpose back into a 1D vector wherein the indices are correctly ordered to
               ## restitch each BANK's 160 freq elements
               stitchIdxArr = np.reshape(reshapeIdxArr.T, (1,160))
               stitchIdxArr = stitchIdxArr.flatten()
               ## finally, loop through 40x40x160 cube,apply fft shift to 32 channel chunks, and reverse indices  
               ## create new cube to sort into and 1D array to hold correct indices
               newCube = np.zeros([numElem, numElem, numFreqs], dtype= 'complex64')
               correctIdxArr = np.zeros([160])
               for idx in range(0, 5):
                   ## get 32 channel chunk to do fftshift on indices
                   chunk = np.fft.fftshift(stitchIdxArr[idx*32:idx*32+32])
                   ## reverse indices in chunk to put in correct order
                   revChunk = chunk[::-1]
                   correctIdxArr[idx*32:idx*32 + 32] = revChunk
               ## finally, loop through cube to re-order freq channels
               for idx in range(0, 160):
                   corrIdx = correctIdxArr[idx]
                   newCube[:,:,idx] = retCube[:,:,corrIdx]
               return newCube
           else:
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
        
    def getSpectralArray(self, fitsName, dataArr, beam, xID, bank):
        ## get number of freq channels
        numFreqs = dataArr.shape[1] / 2112 ## always 2112 complex pairs per frequency channel
        ## get CHANSEL for if we are in PFB mode
        corrHDU = fits.open(fitsName)
        chanSel = np.int(corrHDU[0].header['CHANSEL'])
        ## create bandpass containers. Rows are ints; colums represent integration bandpasses
        spectrumArr_X = np.zeros([len(dataArr[:,0]),numFreqs], dtype='float32')    
        spectrumArr_Y = np.zeros([len(dataArr[:,0]),numFreqs], dtype='float32') 
        ## unpack the weights from the associated binary file
        #if xID == 5:
        #    self.ints=1
        xWeight,yWeight = self.getWeights(numFreqs, xID)
        
        ## DEBUG (build weight bandpass)
        weightXBP = np.zeros([25], dtype='complex64')
        weightYBP = np.zeros([25], dtype='complex64')
        ## DEBUG

        ## loop through integrations and process covariance bandpass to beam-formed spectra
        for ints in range(0,len(dataArr[:,0])):   
            ## TODO: best place for this??
            self.progressBar(ints,len(dataArr[:,0]), beam, xID)            
            ## DEBUG
            #if ints == 99:
            #    self.ints= 1
            ## DEBUG
            
            ## grab a covariance bandpass for a single integration
            dataVector=dataArr[ints,:]
            ## send to get in FISHFITS order before sorting into 'cube' of shape (40,40,freqChans)
            corrCube = self.getCorrelationCube(dataVector, numFreqs)
            
            ## loop through frequency channels to apply weights. Inputs into processPol
            ## are a 40x40 covariance matrix, pol flag, and a length 40 weight vector
            cnt = 0
            absWtIdx = 0
            for z in range(0,numFreqs):
                dat = corrCube[:,:,z]
                if numFreqs == 25:
                    spectrumArr_X[ints,z] = np.abs(self.processPol(dat, xWeight[z,:,beam])) ##0 is XX
                    spectrumArr_Y[ints,z] = np.abs(self.processPol(dat, yWeight[z,:,beam])) ##1 is YY
                    ## DEBUG
                    ## we wish to produce a diagnostic plot which shows a max weight (as determined by the real component)
                    if ints == 0: 
                        weightXBP[z] = np.max(xWeight[z,0,beam])
                        weightYBP[z] = np.max(yWeight[z,20,beam])
                    ## DEBUG
                                           
                elif numFreqs == 160:
                    wtIdx = absWtIdx + (chanSel*5)
                    xWeightIn = xWeight[wtIdx,:,beam]
                    yWeightIn = yWeight[wtIdx,:,beam]
                    spectrumArr_X[ints,z] = np.abs(self.processPol(dat, xWeightIn)) ##0 is XX
                    spectrumArr_Y[ints,z] = np.abs(self.processPol(dat, yWeightIn)) ##1 is YY
                    cnt += 1
                    
                    if cnt == 32:
                        absWtIdx += 1
                        cnt = 0
                    #if ints == 0:
                    #    weightXBP[z] = np.max(np.abs(xWeightIn))
                    #    weightYBP[z] = np.max(np.abs(yWeightIn))
                    ## DEBUG         
        return spectrumArr_X,spectrumArr_Y, weightXBP, weightYBP
        
