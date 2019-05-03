"""
3/4/19
Script to produce the raw time averged spectrum for (i.e., averages all integartions) for a user provided scan for FLAG 
The spectrum is constructed from a user provided element and has units of 'counts'
Inputs:
projectName - name of your GBT project and session (e.g. AGBT16B_400_14)
timeStamp - time stamp of scan
Element - Dipole element (integer 1 - 19)

Usage:
ipython plotBandpass.py projectName timeStamp Element 
Example:
ipython plotBandpass.py AGBT16B_400_14 2017_05_27_05:13:44 1
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "Nickolas.Pingel@anu.edu.au"
__status__ = "Production"
"""

## imports
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib 
import os
import glob
import sys
import pickle
from scipy.optimize import curve_fit
matplotlib.rc('font', family='serif') 
matplotlib.rc('font', serif='Avant Garde') 
matplotlib.rc('text', usetex='False') 
matplotlib.rcParams.update({'font.size': 14})


# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

## grab user inputs
projectID = sys.argv[1] ## GBO project ID in form of AGBT16B_400
timeStamp = sys.argv[2] ## 'On' timestamp in format 2018_01_01_00:00:00
elem = sys.argv[3] ## needs to be an integer between 1 and 19 

# ERROR HANDLING
"""
First, let's check the correct number of inputs were provided
"""
if not len(sys.argv[1:]) == 3:
	print('Incorrect number of user provided inputs; usage: python plotBandpass.py Project FITSfilename timeStamp Element')

"""
make sure the project string is valid. Any ancillary FITS files will be in /home/gbtdata/projectID
"""    
#if not os.path.exists('/home/gbtdata/' + projectID):
#if not os.path.exists('/home/archive/science-data/17B/' + projectID):
#	print('Incorrect path to project directory. Exiting...')
#	sys.exit(1)

"""
check time stamp format and if it exists
"""
if not (timeStamp[4] == '_') and not (timeStamp[7] == '_') and not (timeStamp[-6] == ':') and not (timeStamp[-3] == ':'):
	print('Incorrect format of the time stamp; should be YYYY_MM_DD_HH:MM:SS')
	sys.exit(1)

#if not os.path.isfile('/home/gbtdata/' + projectID + '/Antenna/' + timeStamp + '.fits'):
#if not os.path.isfile('/home/archive/science-data/17B/' + projectID + '/Antenna/' + timeStamp + '.fits'):
#if not os.path.isfile('/users/npingel/' + timeStamp + '.fits'):
#	print('This is not a recorded scan; exiting...')
#	sys.exit(1)

"""
finally, check formatting of input dipole
"""
intElem = np.int(elem)
if (1 > intElem > 19):
	print('No such dipole exists; range of input must be integer between 1 and 19. Exiting...')
	sys.exit(1)

## Dictionaries to obtain desired element index in 1D covariance data. 
elemMapDict = {'1_X':0, '1_Y': 260, '2_X':3, '2_Y': 263, '3_X': 8, '3_Y': 308, '4_X': 11, '4_Y': 311,
'5_X': 20, '5_Y': 360, '6_X': 23, '6_Y': 363, '7_X': 36, '7_Y': 416, '8_X': 39, '8_Y': 419, '9_X': 56, '9_Y': 476,
'10_X': 59, '10_Y': 479, '11_X': 80, '11_Y': 540, '12_X': 83, '12_Y': 543, '13_X': 108, '13_Y': 608, '14_X': 111, '14_Y': 611,
'15_X': 140, '15_Y': 680, '16_X': 143, '16_Y': 683, '17_X': 176, '17_Y': 756, '18_X': 179, '18_Y': 759, '19_X': 216, '19_Y': 836}

## progress bar
def progressBar(value, endvalue,bar_length=20):

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent of FITS files processed: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()
	  
## grab element Index
elemX = elem + '_X'
elemY = elem + '_Y'
elemIdx_X = elemMapDict[elemX]
elemIdx_Y = elemMapDict[elemY]

## define data directories
dataDir = '/lustre/flag/' + projectID + '/BF/'
#dataDir = '/users/npingel/'

fitsList = glob.glob(dataDir + timeStamp+'*.fits')

## check to see if the FITS files were found
if (len(fitsList) == 0):
	print('No data files were found; check dataDir variable')
	sys.exit(1)

"""
loop through on/off banks to construct and extract bandpass based on XID number. If in calibration mode 
(i.e., 25 total freq channels in a bank), the pattern of frequency channels is then: BANKA contains coarse channels 0-4,100-104,200-204, 300-304, 400-404
BANKB contains coarse channels 5-9, 105-109, 205-209, 305-309, 405-409... If in fine channelized mode (i.e., 3200 total frequency channels per BANK), BANKA 
contains fine channels 0-159, BANKB contains fine channels 160 - 319, ... BANK T contains fine channels 3039 to 3199.
"""
xidArr = np.zeros([len(fitsList)]) ## define array to hold XID value of BANK
for i in range(0, len(fitsList)):
	
	## update progress Bar
	progressBar(i,len(fitsList)-1)

	## open files
	dataHdu = fits.open(fitsList[i])

	## get correlation data
	data = dataHdu[1].data['DATA']
	
	#Only process banks with good data
	if len(data) > 0:

		## get XID from either header; should be static over session
		xid = int(dataHdu[0].header['XID'])
		xidArr[i] = xid ## add to xidArr for book keeping
		
		"""
		get element power and integrations; rows correspond to discrete time samples (i.e. integrations); 
		columns are frequency channels ordered as described above
		"""
		freqChanTimeSeries_XX = data[:, elemIdx_X::2112]
		freqChanTimeSeries_YY = data[:, elemIdx_Y::2112]
		numInts = freqChanTimeSeries_XX.shape[0]
		numChans = freqChanTimeSeries_XX.shape[1]
		
		## create full bandpass array if first loop iteration
		if i == 0:

			## on's
			fullBandpassArr_XX = np.zeros([numInts, numChans*20])
			fullBandpassArr_YY = np.zeros([numInts, numChans*20]) 
			
		## based on the mode, we need to alter how the bank files are indexed. For coarse channel, the numChans will equal 25
		if numChans == 25:
			## loop over integrations to sort into On/Off full bandpass array based on XID. Sorting is done before moving on to next integration
			for j in range(0, numInts):
				
				## starting bandpass index is determined by XID
				bandpassStartChan = xid*5
				bandpassEndChan = bandpassStartChan+5
				for z in range(0,5):
					
					## bankStartChan/bankEndChan are multiples of 5
					bankStartChan = z*5
					bankEndChan = bankStartChan + 5
					fullBandpassArr_XX[j, bandpassStartChan:bandpassEndChan] = np.absolute(freqChanTimeSeries_XX[j,bankStartChan:bankEndChan])
					fullBandpassArr_YY[j, bandpassStartChan:bandpassEndChan] = np.absolute(freqChanTimeSeries_YY[j,bankStartChan:bankEndChan])
					
					## increment bandpassStartChan/bandpassEndChan by 100 for proper position in full bandpass
					bandpassStartChan += 100
					bandpassEndChan = bandpassStartChan+5
			
			## all 25 coarse channels used
			bandWidth = .30318*25*20		
		
		## now process if in fine channelization mode 
		if numChans == 160:

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

        		correctIdxArr = np.zeros([160])
        		for idx in range(0, 5):
            			## get 32 channel chunk to do fftshift on indices
            			chunk = np.fft.fftshift(stitchIdxArr[idx*32:idx*32+32])

				## reverse indices in chunk to put in correct order contiguous order
				revChunk = chunk[::-1]
				correctIdxArr[idx*32:idx*32 + 32] = revChunk
			
			"""
			loop over integrations to average bandpass in time, then sort into On/Off full bandpass array based on XID and CHANSEL
			obtain CHANSEL, or which set of 5 contigious channels were selected for PFB
			"""
			chanSel = int(dataHdu[0].header['CHANSEL'])
			for j in range(0, numInts):

            			newTimeSeries_XX = np.zeros(numChans, dtype= 'complex64')
            			newTimeSeries_YY = np.copy(newTimeSeries_XX)

				## finally, loop through cube to re-order freq channels
                		for idx in range(0, 160):
                			corrIdx = correctIdxArr[idx]
                			newTimeSeries_XX[idx] = freqChanTimeSeries_XX[j, corrIdx]
                			newTimeSeries_YY[idx] = freqChanTimeSeries_YY[j, corrIdx]

				bandpassStartChan = xid*160
				bandpassEndChan = bandpassStartChan + 160
				
				## no need to sort here since all elements in banks are contigious in this mode
				fullBandpassArr_XX[j, bandpassStartChan:bandpassEndChan] = np.absolute(newTimeSeries_XX)
				fullBandpassArr_YY[j, bandpassStartChan:bandpassEndChan] = np.absolute(newTimeSeries_YY)
				
			## only five coarse channels selected in PFB; centralFreq depends on chunk of channels selected
			bandWidth = .30318*5*20

## average bandpasses over integrations
fullBandpassArr_XX = np.mean(fullBandpassArr_XX, axis=0)
fullBandpassArr_YY = np.mean(fullBandpassArr_YY, axis=0)

## frequency axis assuming restfreq of 1450 [MHz]; bandWidth is set in processing above
cenFreq = 1450. ##[MHz]
freqAxis = np.linspace(cenFreq-bandWidth/2, cenFreq+bandWidth/2, num = numChans*20)
## if in PFB mode, shift central Freq
if numChans == 160:
	minXid = np.min(xidArr)
	maxXid = np.max(xidArr)
	nu0_Coarse = cenFreq - 75.795
	nu0_Fine = nu0_Coarse + (chanSel*100)*.30318
	nu1_Fine = nu0_Fine + 100*0.30318
	cenFreq = nu0_Fine + (nu1_Fine - nu0_Fine)/2
	freqAxis = np.linspace(cenFreq-bandWidth/2, cenFreq+bandWidth/2, num = numChans*20)

##inform user that we're plotting
print('\n')
print('Plotting bandpass...')

## Plot data and fit
pyplot.figure()
pyplot.plot(freqAxis, fullBandpassArr_XX, linewidth = 2, color = tableau20[2], label = 'XX-Pol')
pyplot.plot(freqAxis, fullBandpassArr_YY, linewidth = 2, color = tableau20[0], label = 'YY-Pol')
#pyplot.xlim(1418, 1423)
pyplot.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pyplot.ylabel(r'Counts')
pyplot.title('Dipole: ' + np.str(intElem) + '; ' + timeStamp, fontsize='small' )
pyplot.xlabel('Frequency [MHz]')
pyplot.axhline(0, color='black', linewidth=2)
pyplot.legend(loc = 0, prop = {'size':15})
#pyplot.savefig(timeStamp + '_bandpass.pdf')
pyplot.show()
