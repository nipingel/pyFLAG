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

def getRawCorrelations(fitsName):
    hdu = fits.open(fitsName)
    corrData = hdu[1].data.field('DATA')
    return corrData
##TODO:
##What is the format of the weight files?        
def getWeights():
    return    

def processPol(corrMatrix):
    R = np.copy(corrMatrix[0:19, 0:19]) ##drop empty data stream
    Reval, Revec = np.linalg.eig(R)
    
    aomega = Revec[:,Reval==max(Reval)] ##TODO: should this be here?
    
    ##TODO:
    ##where do we read weights from?
    w = getWeights() 
    spectrum = np.dot(w.conj().T,np.dot(R,w))
    return spectrum
## Put 1D correlation array into 20x20 matrix (element 20,20 irrelevant)
def unpackCorrelations(dat):
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

def getSpectrum(data): ##data correlations in FISHFITS order for a single frequency channel
    corrMatrix_X,corrMatrix_Y = unpackCorrelations(data)
    xSpec = processPol(corrMatrix_X)
    ySpec = processPol(corrMatrix_Y)
    return xSpec, ySpec