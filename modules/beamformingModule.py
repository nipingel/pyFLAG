# -*- coding: utf-8 -*-
"""
11/01/19
This module contains methods/functions for beamforming. The raw 1D correlatoins and weights are read 
in using the passed in paths. The correlations are sorted into 3D matrices (elemXelemXintegrations). 
The beamforming weights are then applied to each plane of this 3D array. The resulting array is passed back
to PAF_Filler.py to sort the frequency channels and collate metadata. The method that drives the beamforming, 
getSpectralArray(), is called from PAF_Filler.py.

Inputs from PAF_Filler.py for task getSpectralArray:
fitsName - The scan name (i.e. time stamp) that is currently being processed
data - correlation matrices in the form of integrationsx(corrs * numFreqs)
beam - beam currently being processed
xID - xID of BANK being processed

__email__ = "Nickolas.Pingel@anu.edu.au"
__status__ = "Production"
"""

##imports
from astropy.io import fits
import numpy as np
import sys
import glob
import matplotlib.pyplot as pyplot

## class to read in BANK data from FITS; processes each integration to go from raw covariances to beam-formed spectra
class BeamformingModule:
    ## function to initialize object
    ## grabs the map from covariances to FISHFITS order and assigns path to data and weight files
    def __init__ (self, dataPath, weightPath):
        self.mapVector = np.loadtxt('/users/npingel/FLAG/SpectralFiller/misc/gpuToNativeMap.dat', dtype='int')
        self.dataPath = dataPath
        self.weightPath = weightPath
        dataPathSplit = dataPath.split('/') ## get project ID
        self.projectID = dataPathSplit[-3] ## always third last
            
    ## function for progress bar that informats user              
    def progressBar(self,value, endvalue, beam, xID,bar_length=20):
        beamName = str(beam)
        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))
        sys.stdout.write("\rPercent of integrations filled in beam "+beamName+",xID "+str(xID)+": [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()    
    
    ## function to return raw data from the bank FITS files
    def getRawCorrelations(self,fitsName):
        hdu = fits.open(fitsName)
        corrData = hdu[1].data.field('DATA')
        return corrData      
    
    ## function to extract beamforming weights. 
    ## Returns: complex weight vector in format that can be applied to corr. matrices (i.e. 2D array with freqChansXdipoles)
    def getWeights(self, numChans, xID, beam): ## TODO: finish weights
        weightFileList = glob.glob(self.weightPath)
        ## check to make sure we have the correct file (in case A does not map to xID 0, B to 1, etc... )
        for wtFile in weightFileList:
            wtHDU = fits.open(wtFile)
            instID = wtHDU[0].header['XENGINE']
            if instID == xID:
                break

        ## weight array in the form of freqchans, element data, beams
        fullDataArr = wtHDU[0].data ## grab data
        xWeights = np.zeros([25, 40], dtype='complex64')  ## define arrays to hold XX and YY weights
        yWeights = np.zeros([25, 40], dtype='complex64')  ## shape is freqChansxDipoles

        ## extract beam data
        beamXStr = 'Beam' + beam + 'X'
        beamYStr = 'Beam' + beam + 'Y'
        wtDataX = wtHDU[1].data[beamXStr]
        wtDataY = wtHDU[1].data[beamYStr]
        
        ## 64 complex pairs per freq channel
        totalElemInChan = 64 * 2
        
        ## can drop last 24 correlation pairs as they are unused and zero. First 80 elements are good. 
        ## fill in weight arrays by looping through freq channels
        for chan in range(0,25):
          weightPairsX = wtDataX[(totalElemInChan*chan):(totalElemInChan*chan) + 80]                 
          xWeights[chan,:].real = weightPairsX[0::2]
          xWeights[chan,:].imag = weightPairsX[1::2]
          weightPairsY = wtDataY[(totalElemInChan*chan):(totalElemInChan*chan) + 80]                 
          yWeights[chan,:].real = weightPairsY[0::2]
          yWeights[chan,:].imag = weightPairsY[1::2] 
        return xWeights, yWeights
    
    """
    function to perform dot product between the Hermition weight vector, correlation matrix, and weight vector.
    Returns: beam formed spectrum for single frequency channel
    """
    def processPol(self, R, w):
        spectrum = np.dot(w.conj().T,np.dot(R,w))
        return spectrum
    
    """    
    function that reorders correlation bandpass to FISHFITS order;
    inputs are a covariance bandpass for a single integration and number of freq chans
    Returns: 3D data cube with dimensions DipoleXDipoleXFreqChans
    """
    def getCorrelationCube(self, dataVector, numFreqs):
      """
      Reorder to FISHFITS
      explicit order is described below, but consists of 820 complex pairs*numFreqs
      """
      numCorr = 820
      newDataVector= np.zeros([820 * numFreqs], dtype='complex64')
       
      ## set starting index for new data buffer
      FITS_strt_idx = 0
      
      ## exapand the map vector to the same size of the 1D correlation vector
      expMapVector = np.array([self.mapVector + 2112*i for i in range(numFreqs)]).flatten()

      newDataVector = dataVector[expMapVector]
      """
      sort correlations into correct FISHFITS order (dp1Xdp1_freq1, dp1Xdp2_freq1,  ... dp1Xdp40_freq1
      dp2xdp2_freq1, ... dp2xdp40_fre1, dp3xdp3_freq1, ... dp40xdp40_freq1, dpxdp1_freq2, ... dp40x4dp0_freq2,
      ... dp40xdp40_freqN, where N is the total number of freq channels. 
      """
      #for i in range(0, numFreqs * 2112, 2112):
        #singleChanData = dataVector[i: i + 2112] ## grab a single freq. channel worth of correlations
           
       # """
       # loop through and assign new order based on the mapVector attribute, which gives the new index of the 
       # raw correlation
       # """ 
        #for z in range(0,len(self.mapVector)):
        #    newDataVector[z+FITS_strt_idx] = singleChanData[self.mapVector[z]]
        #FITS_strt_idx+=820 ## update starting index
       
      ## Once sorted, we can construct a retCube with dim1 = numElements, dim2 = numElements, dim3 = freqChans
      numElem = 40
      retCube = np.zeros([numElem, numElem, numFreqs], dtype=np.complex64)
      offset = 0
       
      """
      we wish to fill out the cube such that 0,0,0 (numElem, numElem, freqChan) is equal to ch1Xch1_freq1.
      Since the correlations are redundant, we can loop through columns while decreasing the number of rows 
      looped through each column iteration by one. (filling in 0,0 ... 39,0; new column: 1,1 ... 39,1). Using the fact that 820 corr pairs 
      make up each freq channel, once on a specific row, col element, we select each 820 pair to fill in the frequency axis. The 'offset', 
      which denotes the correlation pair, must be increased by one after EVERY iteration. The lower index of rows is increased to avoid placing a 
      redundant correlation.
      Returns: correlation cube (of shape DipolesxDipolesxfreqChans)
      """
      for col in range(numElem) :
        for row in range(col, numElem) :
          retCube[row, col, :] = newDataVector[offset::numCorr]
          if row != col :
            retCube[col,row,:] = newDataVector[offset::numCorr].conj()
          offset += 1   
       
      ## if in PFB mode (i.e. 160 fine channels), we need to re-stitch the frequency channels as they are
      ## output in a non-contiguous fashion
      if numFreqs == 160:
        # make 1D vector containing current indices for single 160 channel chunk
        origIdxArr = np.linspace(0, 159, 160, dtype='int32')
           
        ## reshape into a 32 (rows) x 5 (cols) array; each column contains index for one coarse channels worth of data
        reshapeIdxArr = np.reshape(origIdxArr, (32,5))
        
        """   
        reshape the transpose back into a 1D vector wherein the indices are correctly ordered to
        restitch each BANK's 160 freq elements
        """
        stitchIdxArr = np.reshape(reshapeIdxArr.T, (1,160))
        stitchIdxArr = stitchIdxArr.flatten()
        
        """   
        finally, loop through 40x40x160 cube, apply fft shift to 32 channel chunks, and reverse indices  
        create new cube to sort into and 1D array to hold correct indices
        """
        newCube = np.zeros([numElem, numElem, numFreqs], dtype= 'complex64')
        correctIdxArr = np.zeros([160])
        for idx in range(0, 5):
          ## get 32 channel chunk to do fftshift on indices
          chunk = np.fft.fftshift(stitchIdxArr[idx*32:idx*32+32])
               
          ## reverse indices in chunk to put in correct order contiguous order
          revChunk = chunk[::-1]
          correctIdxArr[idx*32:idx*32 + 32] = revChunk
           
        ## finally, loop through cube to re-order freq channels
        #for idx in range(0, 160):
        #  corrIdx = np.int(correctIdxArr[idx])
        #  newCube[:,:,idx] = retCube[:,:,corrIdx]
        newCube = retCube[:, :, correctIdxArr]
        ## return the correlation cube
        return newCube
      else:
        ## return the correlation cube
        return retCube

    """
    Function that was leftover for use from Anish's backend. Interesting to see how the polarizations were
    mapped to the dipole index. I am therefore keep this in as a reference should it ever be useful. 
    """ 
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
        
        xdat = np.empty_like(dat[0:20,0:20])
        ydat = np.empty_like(dat[0:20,0:20])
        xdat = dat[0:19]
        ydat = dat[20:39]
        return xdat, ydat 
    
    """
    This is the main function of this module (i.e., makes all the necessary calls to return a beamformed
    bandpass). This is called from PAF_Filler.py. After obtaining the weights, the bandpasses are beamformed
    by individual integrations. 
    Returns: beamformed spectra for XX and YY polarizations (shape: integrations X freqChans)
    """
    def getSpectralArray(self, fitsName, beam):

      ## bank name is ALWAYS sixth-to-last character in string
      bank = fitsName[-6]
      corrHDU = fits.open(fitsName)

      ## get relevant information from header
      nRows = corrHDU[1].header['NAXIS2']
      dataArr = corrHDU[1].data['DATA']
      xID = np.int(corrHDU[0].header['XID'])
      
      """
      get CHANSEL for if we are in PFB mode. This selects which chunk of coarse channels were
      sent through PFB. For example, if CHANSEL is 0, coarse channels 0-99 were selected; if
      CHANSEL is 1, coarse channels 100-199 were selected.
      """
      chanSel = np.int(corrHDU[0].header['CHANSEL'])

      ## get number of freq channels
      numFreqs = np.int(dataArr.shape[1] / 2112) ## always 2112 complex pairs per frequency channel
        
      ## create bandpass containers. Rows are ints; colums represent bandpasses for that integration
      spectrumArr_X = np.zeros([len(dataArr[:,0]),numFreqs], dtype='float32')    
      spectrumArr_Y = np.zeros([len(dataArr[:,0]),numFreqs], dtype='float32') 
      
      ## unpack the weights from the associated binary file
      xWeight, yWeight = self.getWeights(numFreqs, xID, beam)

      ## loop through integrations and process covariance bandpass to beam-formed spectra
      for ints in range(0, len(dataArr[:,0])):   
        ## update progress bar for BANK file
        self.progressBar(ints, len(dataArr[:,0]), beam, xID)           
        
        ## grab a covariance bandpass for a single integration
        dataVector=dataArr[ints, :]
        
        ## send to get in FISHFITS order before sorting into 'cube' of shape (40, 40, freqChans)
        corrCube = self.getCorrelationCube(dataVector, numFreqs)
        
        ## loop through frequency channels to apply weights. Inputs into processPol
        ## are a 40x40 covariance matrix, and a length 40 weight vector
        cnt = 0
        absWtIdx = 0
        for z in range(0, numFreqs):
          dat = corrCube[:, :, z]
          if numFreqs == 25:
            spectrumArr_X[ints, z] = np.abs(self.processPol(dat, xWeight[z, :])) ##0 is XX
            spectrumArr_Y[ints, z] = np.abs(self.processPol(dat, yWeight[z, :])) ##1 is YY   
            """
            Since every 32 contiguous fine channels corresponds to one coarse channel, the weight corresponding to
            that coarse channel should be applied uniformly over the fine channel bandpass. 
            """             
          elif numFreqs == 160:
            wtIdx = absWtIdx + (chanSel * 5) ## index for selecting correct weight to apply to the next 32 fine channels
            xWeightIn = xWeight[wtIdx, :]
            yWeightIn = yWeight[wtIdx, :]
            spectrumArr_X[ints, z] = np.abs(self.processPol(dat, xWeightIn))
            spectrumArr_Y[ints, z] = np.abs(self.processPol(dat, yWeightIn))
            cnt += 1 ## increase number of fine channels processed

            ## if we have processed 32 channels, select the weight for the next coarse channel 
            if cnt == 32:
              absWtIdx += 1
              cnt = 0 
        
      ## return beamformed bandpasses       
      return spectrumArr_X, spectrumArr_Y
        
