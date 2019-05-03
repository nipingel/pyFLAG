"""
3/4/19
Script to analyze pointing scans for the GBO-PAF. This procedure will plot the summed power (in XX and YY polarization) of a user specified element data as a function
of antenna movement in terms of X-El/El pointing offset. Coordinates of Scan must be in an 'Encoder' in order to obtain correct offsets. 
A Gaussian is fit with the fitted mean being the ideal local pointing correction (LPC). In addition, the estimated system temperature is 
written out to the terminal. The FITS file timestamp must be supplied as well as XEL/EL to determine which keyword to grab from the Antenna FITS file. 
The entire inputs are summarized below 
Inputs:
Project - GBO project name 
Timestamp - Time stamp of scan (e.g. 2018_01_18_00:00:00) 
Scan Direction - Either EXL or EL. Necessary to grab the correct information from FITS header
Element - Dipole element ranging from 1 to 19
Cal_Source - the standard calibration source observed. See documentation for full list of available calibrator sources
Usage:
ipython PAF_Peak.py Project FITSfilename Scan_Direction (either XEL or EL) Element_Pol Cal_Source
Example:
ipython PAF_Peak.py AGBT16B_400_14 2017_05_24_22:57:54 Az 1 X 3C295
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "nipingel@mix.wvu.edu"
__status__ = "Production"
"""
## imports
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib 
import glob
import os 
import sys
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


##grab file timestamp (first command line argument) and scan direciton (Az or El; second command line argument) and element number (third command line argument)
projectID = sys.argv[1] ## GBO project ID in form of AGBT16B_400
fileTimeStamp = sys.argv[2] ## GBO timestamp in format 2018_01_01_00:00:00
direction = sys.argv[3] ## must be either XEL or EL
elem = sys.argv[4] ## needs to be an integer between 1 and 19 
src = sys.argv[5] ## standard flux calibrator



## ERROR HANDLING
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
#if not os.path.exists('/home/archive/science-data/17B/' + projectID):
	print('Incorrect path to project directory. Exiting...')
	sys.exit(1)

"""
check timestamp format and if it exists
"""
if not (fileTimeStamp[4] == '_') and not (fileTimeStamp[7] == '_') and not (fileTimeStamp[-6] == ':') and not (fileTimeStamp[-3] == ':'):
	print('Incorrect format of timestamp; should be YYYY_MM_DD_HH:MM:SS')
	sys.exit(1)

if not os.path.isfile('/home/gbtdata/' + projectID + '/Antenna/' + fileTimeStamp + '.fits'):
#f not os.path.isfile('/home/archive/science-data/17B/' + projectID + '/Antenna/' + fileTimeStamp + '.fits'):
	print('This is not a recorded scan; exiting...')
	sys.exit(1)

"""
check direction and exit if not XEL or EL
"""
if not direction in ['XEL', 'EL']:
	print(direction)
	print('Incorrect direction argument; must be either XEL or EL. Exiting...')
	sys.exit(1)
"""
check formatting of input dipole
"""
intElem = np.int(elem)
if (1 > intElem > 19):
	print('No such dipole exists; range of input must be integer between 1 and 19. Exiting...')
	sys.exit(1)

"""
finally, check that we are looking at a valid calibrator source
"""
if not src in ['3C295', '3C48', '3C123', '3C128', '3C147', '3C286', '3C353', 'VirgoA']:
	print('Not a valid calibrator source; exiting...')
	sys.exit(1)

## Dictionaries to obtain desired element index in 1D covariance data. 
elemMapDict = {'1_X':0, '1_Y': 260, '2_X':3, '2_Y': 263, '3_X': 8, '3_Y': 308, '4_X': 11, '4_Y': 311,
'5_X': 20, '5_Y': 360, '6_X': 23, '6_Y': 363, '7_X': 36, '7_Y': 416, '8_X': 39, '8_Y': 419, '9_X': 56, '9_Y': 476,
'10_X': 59, '10_Y': 479, '11_X': 80, '11_Y': 540, '12_X': 83, '12_Y': 543, '13_X': 108, '13_Y': 608, '14_X': 111, '14_Y': 611,
'15_X': 140, '15_Y': 680, '16_X': 143, '16_Y': 683, '17_X': 176, '17_Y': 756, '18_X': 179, '18_Y': 759, '19_X': 216, '19_Y': 836} 

## grab element Index for Y Pol
elemStr_Y = elem + '_Y'
elemStr_X = elem + '_X'
elemIdx_Y = elemMapDict[elemStr_Y]
elemIdx_X = elemMapDict[elemStr_X]

"""
Gaussian function used to fit data; inputs are obs position, amplitude, mean, and sigma
"""
def fitGauss(pos, a, m, s):
	return a*np.exp(-(pos-m)**2/(2*s**2))
"""
Compute flux from calibrator based on eqn 1 and Table 5 of Perley and Butler 2016. Freq must be in GHz. 
We assume the freq is 1.45 GHz.
"""
def computeFlux(a0, a1, a2, a3, a4, a5, freq):
	logS = a0+a1*np.log10(freq)+a2*np.log10(freq)**2+a3*np.log10(freq)**3+a4*np.log10(freq)**4+a5*np.log10(freq)**5
	return 10**logS

## dictionary of PerleyButler coefficients
fluxCoeffDict = {'3C48':[1.3253,-0.7553,-0.1914,0.0498,0,0],'3C123':[1.8017,-0.7884,-.1035,-0.0248,0.009,0],
		 '3C138':[1.0088,-0.4981,-0.155,-0.010,0.022,0], '3C147':[1.4516,-0.6961,-.201,0.064,-0.046,0.029],
		 'VirgoA':[2.4466,-0.8116,-0.048,0,0,0],'3C286':[1.2481,-0.4507,-0.1798,0.0357,0,0],
		 '3C295':[1.4701,-0.7658,-0.2780,-0.0347,0.0399,0,0],'3C353':[1.8627,-0.6938,-0.100,-0.032,0,0],
		 'Cygnus A':[3.3498,-1.0022,-0.225,0.023,0.043,0]}
## dictionary of Source Coordinates
srcCoordDict = {'3C48':[1,2]}

## functiont to return indices of minimum ra/dec from source ra/dec coordinates
def getMinEqInds(sRa, sDec, raArr, decArr):
	distArr = np.zeros([len(raArr)])
	for c in range(0, len(raArr)):
		raVal = raArr[c]
		decVal = decArr[c]
		d = np.sqrt((raVal - sRa)**2 + (decVal - sDec)**2)
		distArr[c] = d
	minDist = np.min(distArr)
	## get index of minimum ra/dec values
	minInd = np.where(distArr == minDist)
	return minInd

## compute source Flux
fluxArr = fluxCoeffDict[src]
srcFlux = computeFlux(fluxArr[0], fluxArr[1], fluxArr[2], fluxArr[3], fluxArr[4], fluxArr[5], 1.45)

## define data directories
antDir = '/home/gbtdata/' + projectID + '/Antenna/'
#antDir = '/home/archive/science-data/17B/' + projectID + '/Antenna/'
dataDir = '/lustre/flag/' + projectID +'/BF/'
#dataDir = '/users/npingel/'
## open all BANK files with timestamp 
fitsList = glob.glob(dataDir + '/' + fileTimeStamp + '*.fits')

## for loop to process all BANK files
for i in range(0, len(fitsList)):
        ## TODO: error handling
	dataHdu = fits.open(fitsList[i])
	## get correlation data
	data = dataHdu[1].data['DATA']
	
	## get element power 
	freqChanTimeSeriesY = data[:, elemIdx_Y::2112]
	freqChanTimeSeriesX = data[:, elemIdx_X::2112]

	if i == 0:
		
		## open Antenna FITS file
		antHdu=fits.open(antDir + fileTimeStamp + '.fits')
		
		## create data arrays
		powerArrY = np.zeros([len(freqChanTimeSeriesY[:,0])])
		powerArrX = np.zeros([len(freqChanTimeSeriesX[:,0])])
		
		## get timestamps
		dataDMJD = dataHdu[1].data['DMJD']
		antDMJD = antHdu[2].data['DMJD']
		antRa = antHdu[2].data['RAJ2000']
		antDec = antHdu[2].data['RAJ2000']

		## get index that corresponds to minimum distance 
		minEqInd = getMinEqInds(212.835495, 52.202770, antRa, antDec)
		
		## grab positions
		if direction == 'XEL':
			
			"""
			Default Peak scans use 'Encoder' coordinates. In this setting, the OBSC_AZ/EL represent the observed
			command of the tracking center. The MNT_AZ/EL are encoder samples (at 10 Hz), and represent the
			relative angle between point on the track and main reflector (position of structure), as opposed to 
			position on sky. In 'Encoder', we can recover offsets from tracking center by subtracting these columns. 
			"""
			antElPos = antHdu[2].data['MNT_EL']
			antAzPos = (antHdu[2].data['MNT_AZ'] - antHdu[2].data['OBSC_AZ'])

                        
			## get positions at center of scan, which is needed to remove pointing model contribution
			sobsc_azArr = antHdu[0].header['SOBSC_AZ']
			smntc_azArr = antHdu[0].header['SMNTC_AZ']
                        
			## get pointing model offsets
			pointingModel_Az = smntc_azArr - sobsc_azArr
                        
			## apply them to final offsets to get source centered
			antAzPos = antAzPos - np.float(pointingModel_Az) 
       
                        
			## loop through positions to apply appropriate cos(el) factor
			antPos = np.zeros([len(antAzPos)])
			for pos in range(0,len(antAzPos)):
				antPos[pos] = antAzPos[pos]*np.cos(np.deg2rad(antElPos[pos]))

		if direction == 'EL':
			antPos = antHdu[2].data['MNT_EL'] - antHdu[2].data['OBSC_EL']
			sobsc_elArr = antHdu[0].header['SOBSC_EL']
			smntc_elArr = antHdu[0].header['SMNTC_EL']
			pointingModel_El = smntc_elArr - sobsc_elArr
			antPos = antPos - np.float(pointingModel_El)
                        
	## if we have data, (i.e. the players haven't stalled), loop through integrations to create 'continuum' data
	if not len(freqChanTimeSeriesY) == 0:
		for j in range(0,len(powerArrX)):
			powerArrY[j] += np.sum(np.real(freqChanTimeSeriesY[j,:]))
			powerArrX[j] += np.sum(np.real(freqChanTimeSeriesX[j,:]))

## interpolate DMJD values to continuum data
dataYInterp = np.interp(antDMJD, dataDMJD, powerArrY)
dataXInterp = np.interp(antDMJD, dataDMJD, powerArrX)

##drop the baselines
meanBaseValueY = (np.mean(dataYInterp[30:40])+np.mean(dataYInterp[-20:-10]))/2
meanBaseValueX = (np.mean(dataXInterp[30:40])+np.mean(dataXInterp[-20:-10]))/2

powerArrYScale = dataYInterp - meanBaseValueY
powerArrXScale = dataXInterp - meanBaseValueX

## determine max total power
maxPowerY = np.max(powerArrYScale)
maxPowerX = np.max(powerArrXScale)

if antPos[0] > 0:
	antPos = antPos[::-1]

## determine Tsys estimate (assuming L-Band gain of 1.86 [K/Jy]
G = 1.86
powerRatioY = maxPowerY/meanBaseValueY
powerRatioX = maxPowerX/meanBaseValueX
Tsys_Y = (srcFlux * G)/(powerRatioY + 1)
Tsys_X = (srcFlux * G)/(powerRatioX + 1)

##Noramlize the arrays
powerArrYNorm = powerArrYScale/np.max(powerArrYScale)
powerArrXNorm = powerArrXScale/np.max(powerArrXScale)

##fit data
popt_Y, pcov_Y = curve_fit(fitGauss, antPos[30:-30]*60,powerArrYNorm[30:-30], p0=[1,0,10])
popt_X, pcov_X = curve_fit(fitGauss, antPos[30:-30]*60,powerArrXNorm[30:-30], p0=[1,0,10])

offset_X = antPos[minEqInd] - popt_X[1] 
offset_Y = antPos[minEqInd] - popt_Y[1]  


## get uncertainties in fit
fitUncert_Y = np.sqrt(np.diag(pcov_Y))
fitUncert_X = np.sqrt(np.diag(pcov_X))

## compute Chi-Squared
num_Y = (powerArrYNorm[30:-30]-fitGauss(antPos[30:-30], *popt_Y))
redChi2_Y = (np.sum(num_Y**2)/popt_Y[2]**2)/len(powerArrYNorm-3)
num_X = (powerArrXNorm[30:-30]-fitGauss(antPos[30:-30], *popt_X))
redChi2_X = (np.sum(num_X**2)/popt_X[2]**2)/len(powerArrXNorm-3)

##print out fit statistics
print('-----Fit Statistics (YY Pol)-----')
print('Mean [arcmin]: '+np.str(popt_Y[1])+'+/-'+np.str(fitUncert_Y[1]))
print('FWHM [arcmin]: '+np.str(np.abs(popt_Y[2]*2.335))+'+/-'+np.str(fitUncert_Y[2]*2.335))
print('Reduced Chi-Sq: '+np.str(redChi2_Y))
print('Peak Power: '+np.str(maxPowerY))
print('Estimated Tsys [K]: '+np.str(Tsys_Y))

print('-----Fit Statistics (XX Pol)-----')
print('Mean [arcmin]: '+np.str(popt_X[1])+'+/-'+np.str(fitUncert_X[1]))
print('FWHM [arcmin]: '+np.str(np.abs(popt_X[2]*2.335))+'+/-'+np.str(fitUncert_X[2]*2.335))
print('Reduced Chi-Sq: '+np.str(redChi2_X))
print('Peak Power: '+np.str(maxPowerX))
print('Estimated Tsys [K]: '+np.str(Tsys_X))

## Plot data and fit
pyplot.figure()
pyplot.plot(antPos[30:-30]*60.,powerArrYNorm[30:-30], linewidth=2, label = 'Data (YY)', color = tableau20[0])
pyplot.plot(antPos[30:-30]*60.,fitGauss(antPos[30:-30]*60,popt_Y[0],popt_Y[1], popt_Y[2]), linewidth=2, label='Fit (YY)', color = tableau20[0], linestyle  = '--')
pyplot.plot(antPos[30:-30]*60.,powerArrXNorm[30:-30], linewidth=2, label = 'Data (XX)', color = tableau20[2])
pyplot.plot(antPos[30:-30]*60.,fitGauss(antPos[30:-30]*60,popt_X[0],popt_X[1], popt_X[2]), linewidth=2, label='Fit (XX)', color = tableau20[2], linestyle  = '--')
pyplot.ylabel('Normalized Power')
pyplot.title('Power vs. Angular Position ('+ src + ' ' + direction + '); ' + fileTimeStamp,fontsize='small' )
pyplot.xlabel('Angular Offset [arcmin]')
pyplot.axvline(antPos[minEqInd], linestyle = '--', color='black', label='Location of Source', linewidth=2)
pyplot.axvline((popt_X[1] + popt_Y[1]) / 2, linestyle='--', color='red', label='Ave. Fitted Mean', linewidth=2)	
pyplot.legend(loc=0, prop={'size':15})
pyplot.savefig('PeakScan_' + src + '_' + direction + '_' + fileTimeStamp + '.pdf')
pyplot.show()
