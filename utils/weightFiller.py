# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 21:15:01 2016

@author: npingel
"""

import struct
import numpy as np
from astropy.io import fits

with open('/Users/npingel/Desktop/weights.bin', 'rb') as f:
    data=f.read()
##Actual Weights
weightArr = np.zeros([14,19*25*2],dtype='float32')
##bytes(float)*totalElements*freqChan*totBeam*pol*pair
totBytes = 4*64*25*2*7*2
bytesInBeam = 4*64*25*2
for beam in range(0,14):
    weightIdx = 0
    for chunk in range(0,bytesInBeam,512):
        for idx in range(0,152,4):
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
print('STOP')