"""
5/23/17
Script to analyze On/Off scans for the GBO-PAF. This procedure will extracts the spectrum from a user specified element/polarization pair. It performs a simple On-Off/Off*Tsys calibration to obtain antenna temperature 
Usage:
ipython PlotBandpass.py ON_FITSfilename OFF_FITSFilename Element
Example:
ipython PlotBandpass.py 2017_05_27_05:13:44 2017_05_27_05:14:33 1
__author__ = "Nick Pingel"
__version__ = "0.9"
__email__ = "nipingel@mix.wvu.edu"
__status__ = "Beta"
"""
## imports
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib 
import glob
import sys
import pickle
from scipy.optimize import curve_fit
matplotlib.rc('font', family='serif') 
matplotlib.rc('font', serif='Avant Garde') 
matplotlib.rc('text', usetex='False') 
matplotlib.rcParams.update({'font.size': 14})

## base class error
class Error(Exception):
	"""Base class for other exceptions"""
	pass

##define input custom input exception
class WrongDirectionInput(Error):
	"""Raised when input is not 'Az' or 'El'"""
	pass

##grab file timestamp (first and second command line argument) and element number (third command line argument)
OnFileTimeStamp = sys.argv[1]
OffFileTimeStamp = sys.argv[2]
elem = sys.argv[3] ## needs to be integer in ranging from 1-19

##TODO: exception handling

"""
## error handling if scan direction is specified in correctly
while True:
	try:
		if direction != 'Az' and direction != 'El':
			raise WrongDirectionInput
		break
	except WrongDirectionInput:
		print("Input value NOT Az or El!")
	        direction = input("Enter 'Az' or 'El' scan direction: ")	

"""

## Dictionaries to obtain desired element index in 1D covariance data. 
elemMapDict = {'1_X':0, '1_Y': 260, '2_X':3, '2_Y': 263, '3_X': 8, '3_Y': 308, '4_X': 11, '4_Y': 311,
'5_X': 20, '5_Y': 360, '6_X': 23, '6_Y': 363, '7_X': 36, '7_Y': 416, '8_X': 39, '8_Y': 419, '9_X': 56, '9_Y': 476,
'10_X': 59, '10_Y': 479, '11_X': 80, '11_Y': 540, '12_X': 83, '12_Y': 543, '13_X': 108, '13_Y': 608, '14_X': 111, '14_Y': 611,
'15_X': 140, '15_Y': 680, '16_X': 143, '16_Y': 683, '17_X': 176, '17_Y': 756, '18_X': 179, '18_Y': 759, '19_X': 216, '19_Y': 836} 

## function to fit and subtract bandpass structre
def fitBase(xAxis, bandpass, order, x0,x1,x3,x4):
	subSpecArr = bandpass[x0:x1]
	subSpecArr = np.concatenate([subSpecArr, bandpass[x3:x4]])
	subFreqArr = xAxis[x0:x1]
	subFreqArr = np.concatenate([subFreqArr, xAxis[x3:x4]])
	fit = np.polyfit(subFreqArr, subSpecArr, order)
	p = np.poly1d(fit)
	newSpec = bandpass - p(xAxis)
	return newSpec

## smooth function with boxcar
def boxcar(bandpass, width):
	smoothedBandpass = convolve(bandpass, Box1DKernel(width))
	return smoothedBandpass
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
dataDir = '/lustre/projects/flag/AGBT16B_400_03/BF/'

onFitsList = glob.glob(dataDir+OnFileTimeStamp+'*.fits')
offFitsList = glob.glob(dataDir+OffFileTimeStamp+'*.fits')

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


# loop through on/off banks to construct and extract bandpass based on XID number. 
## The pattern of frequency channels is as such: BANKA contains channels 0-4,100-104,200-204, 300-304, 400-404
## BANKB contains channels 5-9, 105-109, 205-209, 305-309, 405-409. 

## if we are in PFB mode (fine channelization), keep track of an xid array to later identify the frequency range

xidArr = np.zeros([len(onFitsList)])
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
		## add to xidArr for book keeping
		xidArr[i] = xid
		## get element power and integrations; rows correspond to discrete time samples (i.e. integrations); 
		## columns are frequency channels ordered as described above
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
		
		## based on the mode, we need to alter how the bank files are indexed. For PAF_CAL mode (coarse channel), the numChans will equal 25
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
		if numChans == 160:
			## loop over integrations to average bandpass in time, then sort into On/Off full bandpass array based on XID and CHANSEL
			## obtain CHANSEL, or which set of 5 contigious channels were selected for PFB
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
##average bandpasses
fullBandpassArr_On_XX = np.mean(fullBandpassArr_On_XX, axis=0)
fullBandpassArr_On_YY = np.mean(fullBandpassArr_On_YY, axis=0)
fullBandpassArr_Off_XX = np.mean(fullBandpassArr_Off_XX, axis=0)
fullBandpassArr_Off_YY = np.mean(fullBandpassArr_Off_YY, axis=0)



## perform calibration
Tsys = 140.4 ## estimate 
Ta_XX = (fullBandpassArr_On_XX/fullBandpassArr_Off_XX-1)*Tsys
Ta_YY = (fullBandpassArr_On_YY/fullBandpassArr_Off_YY-1)*Tsys
## frequency axis assuming restfreq of 1450 [MHz]; bandWidth is set in processing above
cenFreq = 1450. ##[MHz]
freqAxis = np.linspace(cenFreq-bandWidth/2, cenFreq+bandWidth/2,num = numChans*20)

## if in PFB mode, shift central Freq
if numChans == 160:
	minXid = np.min(xidArr)
	maxXid = np.max(xidArr)
	nu0_Coarse = cenFreq - 75.795
	nu0_Fine = nu0_Coarse + (chanSel*100)*.30318
	nu1_Fine = nu0_Fine + 100*0.30318
	cenFreq = nu0_Fine + (nu1_Fine - nu0_Fine)/2
	freqAxis = np.linspace(cenFreq-bandWidth/2, cenFreq+bandWidth/2,num = numChans*20)


## fit polynomial to baseline 
if numChans == 25:
	baseTa_XX = fitBase(freqAxis,Ta_XX, 3, 100, 120, 160, 180)
        baseTa_YY = fitBase(freqAxis, Ta_YY, 3, 50, 140, 160, 450)
elif numChans == 160:
	baseTa_XX = fitBase(freqAxis, Ta_XX, 3, 850,1475,1900,2056) 
	baseTa_YY = fitBase(freqAxis, Ta_YY, 3, 850,1475,1900,2056)

## smooth bandpass
#smoothTa_XX = boxcar(Ta_XX, 2)
#smoothTa_YY = boxcar(Ta_YY, 2)

##inform user that we're plotting
print('\n')
print('Plotting bandpass...')

##save variables if desired
#with open('/users/npingel/FLAG/2017Reduction/PFBTests/M101_centralDipole_fineVars.pickle', 'wb') as f:
#	pickle.dump([freqAxis, baseTa_XX, baseTa_YY], f)


## Plot data and fit
pyplot.figure()
pyplot.plot(freqAxis,baseTa_XX, linewidth = 2, label = 'XX-Pol')
pyplot.plot(freqAxis, baseTa_YY, linewidth = 2, label = 'YY-Pol')
pyplot.ylabel(r'T$_A$ [K]')
pyplot.title('Fine HI Spectrum of M101 (Central Dipole)',fontsize='small' )
pyplot.xlabel('Frequency [MHz]')
pyplot.axhline(0, color='black', linewidth=2)
pyplot.ylim(-.5,5)
pyplot.xlim(1410,1425)
pyplot.legend(loc = 0, prop = {'size':15})
#pyplot.savefig('/users/npingel/FLAG/2017Reduction/PFBTests/Plots/HI_M101_cenralDipole_GBT16B_400_03_Fine_BC.pdf')
pyplot.show()
