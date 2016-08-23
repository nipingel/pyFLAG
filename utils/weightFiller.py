# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 21:15:01 2016

@author: npingel
"""

import struct
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

with open('/Users/npingel/Desktop/Research/FLAG/data/2016_07_25_04:32:35_xid3_weights.bin', 'rb') as f:
    data=f.read()
##array to hold weights
weightArr = np.zeros([14,64*25*2],dtype='float32')
##bytes(float)*totalElements*freqChan*totBeam*pol*pair
totBytes = 4*64*25*2*7*2
bytesInBeam = 4*64*25*2
for beam in range(0,14):
    weightIdx = 0
    for chunk in range(0,bytesInBeam,512):
        for idx in range(0,256,4):
            absIdx = (beam*bytesInBeam)+chunk+idx 
            weightArr[beam,weightIdx] = struct.unpack('f',data[absIdx:absIdx+4])[0]
            weightIdx+=1
##Metadata
startByte = 179200
offSetArr = np.zeros([2,7],dtype='float32')
for i in range(0,7):
    offsetIdx = startByte+(i*8)
    offSetArr[0,i] = struct.unpack('f',data[offsetIdx:offsetIdx+4])[0]
    offSetArr[1,i] = struct.unpack('f',data[offsetIdx+4:offsetIdx+8])[0]

charStartByte = 179200+(14*4)
calibFilename = struct.unpack('64c',data[charStartByte:charStartByte+64])
beamformMethod = struct.unpack('64c',data[charStartByte+64:charStartByte+128])
xEngineID = struct.unpack('Q',data[-8:])[0]

plt.plot(weightArr[7,:])

##Create FITS file
prihdr = fits.Header()
##TODO: add header parameters
prihdr = fits.Header()
#prihdr.set('CALFILE',calibFilename, 'Calibration Set Filename')
#prihdr.set('BEAMFORM',beamformMethod,'Beamforming algorithm used')
prihdr.set('XENGINE',xEngineID,'X-Engine ID')
prihdu = fits.PrimaryHDU(header=prihdr)



##Make table 
col1 = fits.Column(name='Beam0X', format='1E', array=weightArr[0,:])
col2 = fits.Column(name='Beam1X', format='1E', array=weightArr[1,:])
col3 = fits.Column(name='Beam2X', format='1E', array=weightArr[2,:])
col4 = fits.Column(name='Beam3X', format='1E', array=weightArr[3,:])
col5 = fits.Column(name='Beam4X', format='1E', array=weightArr[4,:])
col6 = fits.Column(name='Beam5X', format='1E', array=weightArr[5,:])
col7 = fits.Column(name='Beam6X', format='1E', array=weightArr[6,:])
col8 = fits.Column(name='Beam0Y', format='1E', array=weightArr[7,:])
col9 = fits.Column(name='Beam1Y', format='1E', array=weightArr[8,:])
col10 = fits.Column(name='Beam2Y', format='1E', array=weightArr[9,:])
col11 = fits.Column(name='Beam3Y', format='1E', array=weightArr[10,:])
col12 = fits.Column(name='Beam4Y', format='1E', array=weightArr[11,:])
col13 = fits.Column(name='Beam5Y', format='1E', array=weightArr[12,:])
col14 = fits.Column(name='Beam6Y', format='1E', array=weightArr[13,:])
col15 = fits.Column(name='BeamOff_AZ', format='1E', array=offSetArr[0,:])
col16 = fits.Column(name='BeamOff_EL', format='1E', array=offSetArr[1,:])
cols = fits.ColDefs([col1, col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16])
tbhdu = fits.new_table(cols)
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto('/Users/npingel/Desktop/Research/FLAG/data/2016_07_25_04:32:35_xid3_weights.fits')
