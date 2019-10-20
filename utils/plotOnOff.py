"""
01/20/18
Script to monitor on/off scans for the GBO-PAF. This procedure extracts the spectrum from a user specified element/polarization pair and Tsys estimate. 
It performs a simple On-Off/Off*Tsys calibration to obtain antenna temperature.
Inputs:
projectName - name of your GBT project and session (e.g. AGBT16B_400_14)
'On' timestamp - time stamp for 'On' scan (e.g. 2018_01_01_00:00:00)
'Off' timestamp - time stamp for 'Off' scan
Element - Dipole element (integer 1 - 19)
Tsys - an estimate of the system temperature in units of K

Usage:
ipython plotOnOff.py projectName OnTimeStamp offTimeStamp Element Tsys
Example:
ipython plotOnOff.py AGBT16B_400_14 2017_05_27_05:13:44 2017_05_27_05:14:33 1 25.0
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "nipingel@mix.wvu.edu"
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
onTimeStamp = sys.argv[2] ## 'On' timestamp in format 2018_01_01_00:00:00
offTimeStamp = sys.argv[3] ## 'Off' timestamp in format 2018_01_01_00:00:00
elem = sys.argv[4] ## needs to be an integer between 1 and 19 
tSys = np.float(sys.argv[5]) ## system temperature in K

# ERROR HANDLING
"""
First, let's check the correct number of inputs were provided
"""
if not len(sys.argv[1:]) == 5:
	print('Incorrect number of user provided inputs; usage: python PAF_Peak.py Project FITSfilename Scan_Direction (either XEL or EL) Element_Pol Cal_Source. Exiting...')
	sys.exit(1)

"""
make sure the project string is valid. Any ancillary FITS files will be in /home/gbtdata/projectID
"""    
if not os.path.exists('/home/gbtdata/' + projectID):
	print('Incorrect path to project directory. Exiting...')
	sys.exit(1)

"""
check ON timestamp format and if it exists
"""
if not (onTimeStamp[4] == '_') and not (onTimeStamp[7] == '_') and not (onTimeStamp[-6] == ':') and not (onTimeStamp[-3] == ':'):
	print('Incorrect format of ON timestamp; should be YYYY_MM_DD_HH:MM:SS')
	sys.exit(1)

if not os.path.isfile('/home/gbtdata/' + projectID + '/Antenna/' + onTimeStamp + '.fits'):
	print('This is not a recorded scan; exiting...')
	sys.exit(1)

"""
check OFF timestamp format and if it exists
"""
if not (offTimeStamp[4] == '_') and not (offTimeStamp[7] == '_') and not (offTimeStamp[-6] == ':') and not (offTimeStamp[-3] == ':'):
	print('Incorrect format of OFF timestamp; should be YYYY_MM_DD_HH:MM:SS')
	sys.exit(1)

if not os.path.isfile('/home/gbtdata/' + projectID + '/Antenna/' + offTimeStamp + '.fits'):
	print('This is not a recorded scan; exiting...')
	sys.exit(1)

"""
check formatting of input dipole
"""
intElem = np.int(elem)
if (1 > intElem > 19):
	print('No such dipole exists; range of input must be integer between 1 and 19. Exiting...')
	sys.exit(1)

"""
finally, check if Tsys is a valid float
"""
if not isinstance(tSys, float):
	print('System temperature input cannot be cast into float value; exiting...')
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

onFitsList = glob.glob(dataDir + onTimeStamp+'*.fits')
offFitsList = glob.glob(dataDir + offTimeStamp+'*.fits')


if (len(onFitsList) == 0) or (len(offFitsList) == 0):
	print('No data files were found; check dataDir variable')
	sys.exit(1)

# determine whether On or Off scans have corresponding banks by grabbing the letter from the FITS name
onBankList = []
offBankList = []

for onName in (onFitsList):
	onBankList.append(onName[-6]) ## sixth last element '.fits' 

for offName in (offFitsList):
       	offBankList.append(offName[-6]) ## sixth last element '.fits'

## bad banks are those which do not have a corresponding partner in the opposing scan type
badBanks = list(set(onBankList) - set(offBankList))

## remove bad banks from which ever list is larger
if len(onBankList) > len(offBankList):
	for badBank in badBanks:
		bankFileName = onFitsList[0]
		
		## replace with bad bank
		bankFileName[-6] = badBank
		onFitsList.remove(bankFileName)
else:
        for badBank in badBanks:
                bankFileName = offFitsList[0]
                
                ## replace with bad bank
                bankFileName[-6] = badBank
                onFitsList.remove(bankFileName)

##sort list to be in same order
onFitsList = sorted(onFitsList)
offFitsList = sorted(offFitsList)

"""
loop through on/off banks to construct and extract bandpass based on XID number. If in calibration mode 
(i.e., 25 total freq channels in a bank), the pattern of frequency channels is then: BANKA contains coarse channels 0-4,100-104,200-204, 300-304, 400-404
BANKB contains coarse channels 5-9, 105-109, 205-209, 305-309, 405-409... If in fine channelized mode (i.e., 3200 total frequency channels per BANK), BANKA 
contains fine channels 0-159, BANKB contains fine channels 160 - 319, ... BANK T contains fine channels 3039 to 3199.
"""

xidArr = np.zeros([len(onFitsList)]) ## define array to hold XID value of BANK
for i in range(0, len(onFitsList)):
	
	## update progress Bar
	progressBar(i,len(onFitsList)-1)

	## open files
	onDataHdu = fits.open(onFitsList[i])
	offDataHdu = fits.open(offFitsList[i])

	## get correlation data
	onData = onDataHdu[1].data['DATA']
	offData = offDataHdu[1].data['DATA']
	
	#Only process banks with good data
	if len(onData) > 0 and len(offData) > 0:

		## get XID from either header; should be static over session
		xid = int(onDataHdu[0].header['XID'])
		xidArr[i] = xid ## add to xidArr for book keeping
		
		"""
		get element power and integrations; rows correspond to discrete time samples (i.e. integrations); 
		columns are frequency channels ordered as described above
		"""
		freqChanTimeSeries_On_XX = onData[:, elemIdx_X::2112]
		freqChanTimeSeries_On_YY = onData[:, elemIdx_Y::2112]
		freqChanTimeSeries_Off_XX = offData[:, elemIdx_X::2112]
		freqChanTimeSeries_Off_YY = offData[:, elemIdx_Y::2112]
		numInts = freqChanTimeSeries_On_XX.shape[0]
		numChans = freqChanTimeSeries_On_XX.shape[1]
		
		## create full bandpass array if first loop iteration
		if i == 0:

			## on's
			fullBandpassArr_On_XX = np.zeros([numInts, numChans*20])
			fullBandpassArr_On_YY = np.zeros([numInts, numChans*20]) 
			
			## off's 
			fullBandpassArr_Off_XX = np.zeros([numInts, numChans*20])
			fullBandpassArr_Off_YY = np.zeros([numInts, numChans*20])
		
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
					fullBandpassArr_On_XX[j, bandpassStartChan:bandpassEndChan] = np.absolute(freqChanTimeSeries_On_XX[j,bankStartChan:bankEndChan])
					fullBandpassArr_On_YY[j, bandpassStartChan:bandpassEndChan] = np.absolute(freqChanTimeSeries_On_YY[j,bankStartChan:bankEndChan])
					fullBandpassArr_Off_XX[j, bandpassStartChan:bandpassEndChan] = np.absolute(freqChanTimeSeries_Off_XX[j,bankStartChan:bankEndChan])
					fullBandpassArr_Off_YY[j, bandpassStartChan:bandpassEndChan] = np.absolute(freqChanTimeSeries_Off_YY[j,bankStartChan:bankEndChan])
					
					## increment bandpassStartChan/bandpassEndChan by 100 for proper position in full bandpass
					bandpassStartChan += 100
					bandpassEndChan = bandpassStartChan+5
			
			## all 25 coarse channels used
			bandWidth = .30318*25*20		
		
		## now process if in fine channelization mode 
		if numChans == 160:
			
			"""
			loop over integrations to average bandpass in time, then sort into On/Off full bandpass array based on XID and CHANSEL
			obtain CHANSEL, or which set of 5 contigious channels were selected for PFB
			"""
			chanSel = int(onDataHdu[0].header['CHANSEL'])
                        for j in range(0, numInts):
                                bandpassStartChan = xid*160
                                bandpassEndChan = bandpassStartChan + 160
				
								## no need to sort here since all elements in banks are contigious in this mode
                                fullBandpassArr_On_XX[j, bandpassStartChan:bandpassEndChan] = np.absolute(freqChanTimeSeries_On_XX[j,:])
                                fullBandpassArr_On_YY[j, bandpassStartChan:bandpassEndChan] = np.absolute(freqChanTimeSeries_On_YY[j,:])
                                fullBandpassArr_Off_XX[j, bandpassStartChan:bandpassEndChan] = np.absolute(freqChanTimeSeries_Off_XX[j,:])
                                fullBandpassArr_Off_YY[j, bandpassStartChan:bandpassEndChan] = np.absolute(freqChanTimeSeries_Off_YY[j,:])
				
			## only five coarse channels selected in PFB; centralFreq depends on chunk of channels selected
			bandWidth = .30318*5*20        	

## average bandpasses over integrations
fullBandpassArr_On_XX = np.mean(fullBandpassArr_On_XX, axis=0)
fullBandpassArr_On_YY = np.mean(fullBandpassArr_On_YY, axis=0)
fullBandpassArr_Off_XX = np.mean(fullBandpassArr_Off_XX, axis=0)
fullBandpassArr_Off_YY = np.mean(fullBandpassArr_Off_YY, axis=0)



## perform calibration
Ta_XX = (fullBandpassArr_On_XX/fullBandpassArr_Off_XX-1) * tSys
Ta_YY = (fullBandpassArr_On_YY/fullBandpassArr_Off_YY-1) * tSys
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
pyplot.plot(freqAxis, Ta_XX, linewidth = 2, color = tableau20[2], label = 'XX-Pol')
pyplot.plot(freqAxis, Ta_YY, linewidth = 2, color = tableau20[0], label = 'YY-Pol')
pyplot.ylabel(r'T$_A$ [K]')
pyplot.title('Dipole: ' + np.str(intElem) + '; On (' + onTimeStamp + '); Off(' + offTimeStamp + ')' ,fontsize='small' )
pyplot.xlabel('Frequency [MHz]')
pyplot.axhline(0, color='black', linewidth=2)
pyplot.legend(loc = 0, prop = {'size':15})
#pyplot.savefig('OnOffExample_FullSize.png')
pyplot.show()
