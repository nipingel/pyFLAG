"""
6/6/18
Script to construct the beam patterns and profiles. This script output two plots: (1) a seven panel plot showing the formed beam patterns on the sky in units of dB
(2) a seven panel plot showing perpandicular profiles with Gaussian fits. The script will report the FWHM of these fits and the esimated area in square arcseconds. 
The script needs access to mat lab files produced by BYU. Specifically, the aggregated grid, weights, and tsys.mat files. These contain the steering vectors for 
each beam and XEL/EL locations. The three inputs are project value (assuming a hard coded location), calibration scan time (either grid or seven), and a path to the 
directory that holds the weight FITS files. 
Inputs are:
-p --project_name - <required> project name and observing session (e.g., AGBT19A_365_03
-c --cal_type - <required> calibration type (i.e., seven or grid)
-w --weights_path - <required> path to weights FITS files


Usage:
ipython plotBeamPatterns.py -p /path/to/steering/vectorXX -c /path/to/steering/vectorYY -w /path/to/weights

Example:
ipython plotBeamPatterns.py -p AGBT16B_400_12 -c grid -w ../../data/AGBT16B_400/AGBT16B_400_12/weight_files/fits_files/scaledWeights/ 
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "nmpingel@wisc.edu"
__status__ = "Production"
"""

## imports
from astropy.io import fits 
import scipy.io
import numpy as np
import argparse
np.seterr(divide='ignore', invalid='ignore')
import glob
import sys
import pickle
from scipy.optimize import curve_fit
from scipy import optimize
from scipy.interpolate import *
import matplotlib.pyplot as pyplot
import matplotlib 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

matplotlib.rc('font', family='sans-serif') 
matplotlib.rc('font', serif='Helvetica Neue') 
matplotlib.rc('text', usetex='false') 
matplotlib.rc('xtick.major.width')
matplotlib.rcParams['contour.negative_linestyle']= 'solid'
matplotlib.rcParams.update({'font.size': 14})
# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib 
# accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

## function for progress bar that informats user              
def progressBar(beamName, chan,  endvalue, bar_length=20):
	percent = float(chan) / endvalue
	arrow = '-' * int(round(percent * bar_length)-1) + '>'
	spaces = ' ' * (bar_length - len(arrow))
	sys.stdout.write("\rPercent of Beam " + np.str(beamName) + " processed: [{0}] {1}%".format(arrow + spaces, np.int(round(percent * 100))))
	sys.stdout.flush()    
"""
1D Gaussian function
"""
def Gauss1D(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

"""
function to select the weight file based on distance from specified location ( 
Xel0, el0). aggWeightFil is the aggregated weight matlab file. 
"""
def getWeights(xel0, el0, aggWeightFile):
	## extract data arrays
	xElArr = aggWeightFile['AZ'] *60
	elArr = aggWeightFile['EL'] *60
	w_agg = aggWeightFile['w_agg']
	dList = [] 
	## loop through beams to compute distance from specified xel and el offsets
	for p in range(0, np.size(xElArr, 0)):
		gridXel = xElArr[p]
		gridEl = elArr[p]
		dXel = (xel0 - gridXel)**2
		dEl =  (el0 - gridEl)**2
		d = np.sqrt(dXel + dEl)
		dList.append(np.float(d))

	"""
	extract the index of the min value and return the weights for the beam. 
	The dimensions of the weight array is dipoles X freq Chans
	"""
	minInd = np.where(dList  == np.min(dList))
	minInd = minInd[0][0]
	weightVals = w_agg[:, minInd, :]
	return weightVals


## parse input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--project_name', help='project name and observing session (e.g., AGBT19A_365_03)', required = True)
parser.add_argument('-c', '--cal_type', help='<required> calibration type (i.e., seven or grid)', required = True)
parser.add_argument('-w', '--weights_path', help='<required> path to weights FITS files', required = True)
args, unknown = parser.parse_known_args()


## make beam dictionary to map from BYU to WVU convention (e.g. BYU 0 -> WVU 1)
wvuBeamDict = {'0':'1', '1':'2', '2':'6', '3':'0', '4':'3', '5':'5','6':'4'}
byuBeamDict = {'1':'0', '2':'1', '6':'2', '0':'3', '3':'4', '5':'5', '4':'6'}
rowLocDict = {0:(0, 1), 1:(0, 3), 2:(2, 0), 3:(2, 2), 4:(2, 4), 5:(4, 1), 6:(4, 3)}
rowLocDict = {0:(2,2), 1:(0,1), 2:(0,3), 3:(2,4), 4:(4,3), 5:(4,1), 6:(2,0)}

## make coordinate dictionary for plotting
yDict = {0:0, 1:1, 2:0, 3:1, 4:0, 5:1,6:0}
xDict = {0:0, 1:0, 2:1, 3:1, 4:2, 5:2, 6:3}

matFilePaths = '/mnt/flag/'



## get project name 
projName = args.project_name

## calibration type 
calType = args.cal_type

## get path to weight files
pathToWeights = args.weights_path

## TODO:
## TAKE THIS OUT --- STUPID ELEMENT MAPPING --- MAKE IT READ FROM TXT FILE
#xElemsIndices = np.array([1, 2, 3, 4, 5, 6, 7, 35, 9, 10, 11, 12, 14, 13, 15, 16, 17, 18, 19]) - 1
#yElemsIndices = np.array([21, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36, 37, 38]) - 1
## load element mapping arrays
proj_str = projName.split('_')[0]+'_'+projName.split('_')[1]
elem_mapping = np.loadtxt('../misc/element_mapping_%s.txt' % (proj_str), dtype = 'int')
xElemsIndices = elem_mapping[0, :]
yElemsIndices = elem_mapping[1, :]

## really, 8, 13, and 14 but need to provide the subsequent inices 
badIndsYY = [9, 14, 15] 

## extract matlab files
## try-catch in case file is not found or file is not matlab binary file
try:
	projNameSplit = projName.split('_')
	aggMatFileXX = scipy.io.loadmat(matFilePaths + '/' + projName + '/BF/mat/' + projName + '_aggregated_grid_X_' + calType + '.mat' )
	aggMatFileYY = scipy.io.loadmat(matFilePaths + '/' + projName + '/BF/mat/' + projName + '_aggregated_grid_Y_' + calType + '.mat' )
	tsysMatFileXX = scipy.io.loadmat(matFilePaths + '/' + projName + '/BF/mat/' + projName + '_Xpol_tsys_' + calType + '.mat' )
	tsysMatFileYY = scipy.io.loadmat(matFilePaths + '/' + projName + '/BF/mat/' + projName + '_Ypol_tsys_' + calType + '.mat' )
	aggWeightFileXX = scipy.io.loadmat(matFilePaths + '/' + projName + '/BF/mat/' + projName + '_aggregated_weights_X_' + calType + '.mat')
	aggWeightFileYY = scipy.io.loadmat(matFilePaths + '/' + projName + '/BF/mat/' + projName + '_aggregated_weights_Y_' + calType + '.mat')
	"""
	get the aggregated steering vector arrays, which
	are shaped as elems X grid points X freq
	"""
	aggArrXX = aggMatFileXX['a_agg']
	aggArrYY = aggMatFileYY['a_agg']

	"""
	get the aggregated steering vector coordinate arrays, basically location
	of grid points.
	"""
	aggElArr = aggMatFileXX['EL'] * 60
	aggXElArr = aggMatFileXX['AZ'] * 60

	"""
	get the Tsys/eta values at each grid point. 
	These arrays are shaped as grid points X freq
	"""
	tsysEtaArrXX = tsysMatFileXX['Tsys_eta']
	tsysEtaArrYY = tsysMatFileYY['Tsys_eta']

except ValueError:
        print('File is not a matlab binary file. Exiting...')
        sys.exit(1)

## extract weights
weightHduList = glob.glob(pathToWeights + '/*.FITS')

"""
create global data buffer to hold the weights, which are shaped as 
elems X beam X freq. 
elems X freq X pol X complex pair
"""
weightArrXX = np.zeros([19, 7, 500], dtype = 'complex64')
weightArrYY = np.zeros([19, 7, 500], dtype = 'complex64')

"""
The weights do not correspond to contigious frequency channels. The 
channels in the BANKA FITS file go as 0-4, 100-104, 200-204, .../, 400-404. The
channels in the BANKB FITS file go as 5-9, 105-109, ..., 405-409. Since the 
steering vectors are ordered in contigious freq. channel order, only weights 
data must be read in and sorted based on the XID number. 
""" 
if True:
	for wtFile in weightHduList:
		print('\n')
		print('Working on weight file: ' + wtFile)
		wtHdu = fits.open(wtFile)
		xid = wtHdu[0].header['XENGINE']

		## save offsets to overplot beam centers
		xElOff = wtHdu[1].data['BeamOff_XEL'][0:7] * 60
		elOff = wtHdu[1].data['BeamOff_EL'][0:7] * 60

		"""
	    loop through beams and XX/YY pols to sort the weights in the correct freq
	    order. 
		"""
		## loop over beams
		for bmNum in range(0, 7):
			bmNumStr = np.str(bmNum)
			wtDataXX = wtHdu[1].data['Beam' + bmNumStr + 'X']
			wtDataYY = wtHdu[1].data['Beam' + bmNumStr + 'Y']

			## 64 complex pairs per freq channel
			totalElemInChan = 64 * 2
	        
			"""
	        can drop last 24 correlation pairs as they are unused and zero. First
	        80 elements are good.
			"""
	        
			"""
	        initialize two counting variables to keep track of total chunks of 
	        contigious channels and the subset of channels within these chunks. 
	        The total count increments after every multiple of
	        five iterations in the subsequent loop, while the subset index resets to
	        zero.
			"""
			totChunkCnt = -1
			chanSubsetIdx = 0
	        ## fill in weight arrays by looping through non-contigious freq channels
			for chan in range(0, 25):

	      		## handle count variable conditions
				if chan % 5 == 0:
					totChunkCnt += 1
					chanSubsetIdx = 0

				initWeightPairsXX = wtDataXX[(totalElemInChan * chan):(totalElemInChan*chan) + 80]
				initWeightPairsYY = wtDataYY[(totalElemInChan * chan):(totalElemInChan*chan) + 80]              

				## make array to hold the full weight vector
				initWeightPairsXX_Complex = np.zeros([40], dtype = 'complex64')
				initWeightPairsYY_Complex = np.zeros([40], dtype = 'complex64')

				## fill the arrays
				initWeightPairsXX_Complex.real = initWeightPairsXX[0::2]
				initWeightPairsXX_Complex.imag = initWeightPairsXX[1::2]
				initWeightPairsYY_Complex.real = initWeightPairsYY[0::2]
				initWeightPairsYY_Complex.imag = initWeightPairsYY[1::2]

				## using imported element maps, select the correct weights. 
				wtVector_XX = initWeightPairsXX_Complex[xElemsIndices]
				wtVector_YY = initWeightPairsYY_Complex[yElemsIndices]

				"""
				based on the XID, sort chunk of five contigious channels in the 
				proper spot in full bandpass array
				"""
				freqChan = totChunkCnt * 100 + (xid * 5) + chanSubsetIdx

				weightArrXX[:, bmNum, freqChan] = initWeightPairsXX_Complex[list(xElemsIndices)]
				weightArrYY[:, bmNum, freqChan] = initWeightPairsYY_Complex[list(yElemsIndices)]
				
				## increment subset channel counter
				chanSubsetIdx += 1

	"""
	initialize array to hold the final pattern with dimensions beams X grid points 
	X freq chans
	"""
	patternArrXX = np.zeros([7, np.size(aggArrXX, 1), np.size(aggArrXX, 2)])
	patternArrYY = np.zeros([7, np.size(aggArrXX, 1), np.size(aggArrXX, 2)])

	"""
	TEST: retrieve the weights as BYU dos
	"""
	for w in range(0, 7):

		## get beam offsets
		xElOffVal = xElOff[w] 
		elOffVal = elOff[w]
		weightArrXX[:, w, :] = getWeights(xElOffVal, elOffVal, aggWeightFileXX)
		weightArrYY[:, w, :] = getWeights(xElOffVal, elOffVal, aggWeightFileYY)

	"""
	Now that the weights are properly sorted, we can proceed with the calculation. 
	This is done by looping through each beam, selecting the corresponding weight 
	vector for each frequency channel, while further selecting a specfic grid point 
	that contains the steering vector of the dipoles. Process the two polarizations 
	together. 
	"""

	## loop over beams
	for bmNum in range(0, 7):
		## select weights and frequency channels for specified beam
		wtArrXX = weightArrXX[:, bmNum, :]
		wtArrYY = weightArrYY[:, bmNum, :]

		## loop over freq channels
		#for chan in range(0, np.size(aggArrXX, 2)):
		for chan in range(100, 101):
			
			## update progress bar 
			if chan % 10 == 0:
				progressBar(bmNum, chan, np.size(aggArrXX, 2)) 

			## loop over grid points
			for pt in range(0, np.size(aggArrXX, 1)):

				## TODO: Mask out low sensitivity areas
				## create weight vector 
				wtVecXX = np.matrix(wtArrXX[:, chan]).T
				wtVecYY = np.matrix(wtArrYY[:, chan]).T

				## if len of these arrays are less than 19, add zeros at bad 
				## index values read in above
				aXX = aggArrXX[:, pt, chan]
				aYY = aggArrYY[:, pt, chan]
				#if len(aXX) < 19:
					#aXX = np.insert(aXX, badIndsXX, 0)
				#if len(aYY) < 19:
					#aYY = np.insert(aYY, badIndsYY, 0)


				## select the steering value of the dipoles at that grid point and 
				## freq channel
				aXX = np.matrix(aXX)
				aYY = np.matrix(aYY)

				## if len of these arrays are less than 19, add zeros at bad 
				## index values read in above


				## select the Tsys/eta value for specific grid point and freq
				tsysEtaXXVal = tsysEtaArrXX[pt, chan]
				tsysEtaYYVal = tsysEtaArrYY[pt, chan]

				## if tsys/eta is > 1e6, mask this out
				if tsysEtaXXVal > 1e6:
					patternArrXX[bmNum, pt, chan] = np.float('nan')
				elif tsysEtaXXVal < 1e6:
					patternArrXX[bmNum, pt, chan] = np.abs(np.matmul(wtVecXX.H, aXX.T))**2
				if tsysEtaYYVal > 1e6:
					patternArrYY[bmNum, pt, chan] = np.float('nan')
				elif tsysEtaYYVal < 1e6:
					patternArrYY[bmNum, pt, chan] = abs(np.matmul(wtVecYY.H, aYY.T))**2					
			
			## normalize patterns to peak
			patternArrXX[bmNum, :, chan] = patternArrXX[bmNum, :, chan] / np.nanmax(patternArrXX[bmNum, :, chan])
			patternArrYY[bmNum, :, chan] = patternArrYY[bmNum, :, chan] / np.nanmax(patternArrYY[bmNum, :, chan])
	## output pickle file of pattern arrays
	saveFile = projName + '_' + calType + '_beamPattern.pkl'
	with open(saveFile, "wb") as f:
		pickle.dump([aggXElArr, aggElArr, patternArrXX, patternArrYY], f)
## load picke file
with open(saveFile, 'rb') as f:
  aggXElArr, aggElArr, patternArrXX, patternArrYY = pickle.load(f)

## formatting for terminal
print('\n')

"""
With the pattern calculated and normalized, we chan check out the beam pattern 
at a specific frequency channel. Loop through to map each beam inividually while
also shifting each beam to the center to create an 'average' beam. 
"""
## set up useful plotting variables
cenFreq = 1450
freqArr = np.linspace(-250 * 0.30318, 250 * 0.30318, 500) + cenFreq
chan = 100
extent = np.array([-0.25, 0.25, -0.25, 0.25]) * 60
xAxis = np.linspace(-0.25 * 60, 0.25 * 60, 200)
yAxis =  np.linspace(0.25 * 60, -0.25 * 60, 200)

## concatenate coordinate arrays
coordVec = np.column_stack((aggXElArr, aggElArr))

## create mesh grid
grid_x, grid_y = np.meshgrid(xAxis, yAxis)

## set plotting parameters
majorYLocFactor = 10
minorYLocFactor = 2

majorYLocator = MultipleLocator(majorYLocFactor)
majorYFormatter = FormatStrFormatter('%d')
minorYLocator = MultipleLocator(minorYLocFactor)

majorXLocator = MultipleLocator(10)
majorXFormatter = FormatStrFormatter('%d')
minorXLocator = MultipleLocator(2)

## make beam pattern plot for both polarizations
for pl in range(0, 2):
	## declare figure object
	fig = pyplot.figure(figsize = (12,12))
	## initialize array to hold interpolated grid
	meanBeamArr = np.zeros([2, 7, 200, 200])
	if pl == 0:
		for cnt in range(0, 7):
			beamIndex = wvuBeamDict[np.str(cnt)]
			## determine row location span
			locSeq = rowLocDict[cnt]
			## select and grid the data
			patternAtChanYY = patternArrYY[cnt, :, chan]
			patternGrid = griddata(coordVec, 10*np.log10(patternAtChanYY), (grid_x, grid_y) , method='linear')

			## find index of maximum value so as to know where to take the cuts
			maxInd = np.where(patternGrid == np.nanmax(patternGrid))
			maxX = maxInd[0][0]
			maxY = maxInd[1][0]

			## extract cuts
			elCut = patternGrid[maxX, :]
			xelCut = patternGrid[:, maxY]
			
			print('Plotting beam ' + beamIndex + ', at frequency %.2f' % freqArr[chan] + ' [MHz]')
			

			##Plot gridded data
			ax = pyplot.subplot2grid((6, 6), locSeq, colspan = 2, rowspan = 2)
			im = ax.imshow(patternGrid, extent=extent, vmin = -40, cmap = 'viridis', aspect = 'equal')
			
			#im = axarr[x, y].imshow(patternGrid, extent=extent, vmin = -40, cmap = 'viridis', aspect = 'equal')

			ax.contour(patternGrid, colors='black', linewidths=2, extent=extent, aspect = 'equal', levels=[-15, -10, -5, -2], origin= 'image')
			ax.scatter(xElOff, elOff, marker = 'x', c = 'red', s = 40, linewidths = 2)
			ax.axvline(xAxis[maxY], color = 'red', linestyle = '--', linewidth = 2)
			ax.axhline(yAxis[maxX], color = 'red', linestyle = '--', linewidth = 2)

			## place pattern in array to compute average later
			meanBeamArr[pl, cnt, :, :] = patternGrid 
			ax.set_title('Beam ' + beamIndex, fontsize = 12)

			ax.yaxis.set_major_locator(majorYLocator)
			ax.yaxis.set_major_formatter(majorYFormatter)
			ax.yaxis.set_minor_locator(minorYLocator)
			ax.xaxis.set_major_locator(majorXLocator)
			ax.xaxis.set_major_formatter(majorXFormatter)
			ax.xaxis.set_minor_locator(minorXLocator)
			ax.tick_params(axis = 'both', which='both', width=2)
		#[left, bottom, width, height]
		fig.subplots_adjust(right=0.9, hspace = 0.3, wspace = 0.3)
		cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
		cb = fig.colorbar(im, cax=cbar_ax, label = '[dB]')
		cb.ax.tick_params(which = 'minor', length = 2)
		cb.ax.tick_params(which = 'major', length = 4)
		cb.ax.minorticks_on()
		fig.text(0.5, 0.04, 'XEL Offset [arcmin]', ha='center')
		fig.text(0.04, 0.5, 'EL Offset [arcmin]', va='center', rotation='vertical')
		pyplot.savefig(projName + '_FormedBeamPatterns_YY.png', bbox_inches='tight')
		pyplot.show(block=True)
		pyplot.clf()
		pyplot.close()
	else:
		for cnt in range(0, 7):
			beamIndex = wvuBeamDict[np.str(cnt)]
			## determine row location span
			locSeq = rowLocDict[cnt]
			## select and grid the data
			patternAtChanXX = patternArrXX[cnt, :, chan]
			patternGrid = griddata(coordVec, 10*np.log10(patternAtChanXX), (grid_x, grid_y) , method='linear')

			## find index of maximum value so as to know where to take the cuts
			maxInd = np.where(patternGrid == np.nanmax(patternGrid))
			maxX = maxInd[0][0]
			maxY = maxInd[1][0]

			## extract cuts
			elCut = patternGrid[maxX, :]
			xelCut = patternGrid[:, maxY]
			
			print('Plotting beam ' + beamIndex + ', at frequency %.2f' % freqArr[chan] + ' [MHz]')
			

			##Plot gridded data
			ax = pyplot.subplot2grid((6, 6), locSeq, colspan = 2, rowspan = 2)
			im = ax.imshow(patternGrid, extent=extent, vmin = -40, cmap = 'viridis', aspect = 'equal')
			
			#im = axarr[x, y].imshow(patternGrid, extent=extent, vmin = -40, cmap = 'viridis', aspect = 'equal')

			ax.contour(patternGrid, colors='black', linewidths=2, extent=extent, aspect = 'equal', levels=[-15, -10, -5, -2], origin= 'image')
			ax.scatter(xElOff, elOff, marker = 'x', c = 'red', s = 40, linewidths = 2)
			ax.axvline(xAxis[maxY], color = 'red', linestyle = '--', linewidth = 2)
			ax.axhline(yAxis[maxX], color = 'red', linestyle = '--', linewidth = 2)

			## place pattern in array to compute average later
			meanBeamArr[pl, cnt, :, :] = patternGrid 
			ax.set_title('Beam ' + beamIndex, fontsize = 12)

			ax.yaxis.set_major_locator(majorYLocator)
			ax.yaxis.set_major_formatter(majorYFormatter)
			ax.yaxis.set_minor_locator(minorYLocator)
			ax.xaxis.set_major_locator(majorXLocator)
			ax.xaxis.set_major_formatter(majorXFormatter)
			ax.xaxis.set_minor_locator(minorXLocator)
			ax.tick_params(axis = 'both', which='both', width=2)
		#[left, bottom, width, height]
		fig.subplots_adjust(right=0.9, hspace = 0.3, wspace = 0.3)
		cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
		cb = fig.colorbar(im, cax=cbar_ax, label = '[dB]')
		cb.ax.tick_params(which = 'minor', length = 2)
		cb.ax.tick_params(which = 'major', length = 4)
		cb.ax.minorticks_on()
		fig.text(0.5, 0.04, 'XEL Offset [arcmin]', ha='center')
		fig.text(0.04, 0.5, 'EL Offset [arcmin]', va='center', rotation='vertical')
		pyplot.savefig(projName + '_FormedBeamPatterns_XX.png', bbox_inches='tight')
		pyplot.show(block=True)
		pyplot.clf()
		pyplot.close()		

"""
Finally, loop through again to create EL/XEL profiles plots of the beams. 
Overplot the series of the four cuts on an individual basis, then overplot all
seven profiles per cut on individual plots
"""

## set plotting parameters
majorYLocFactor = 0.2
minorYLocFactor = 0.2/4

majorYLocator = MultipleLocator(majorYLocFactor)
majorYFormatter = FormatStrFormatter('%.2f')
minorYLocator = MultipleLocator(minorYLocFactor)

majorXLocator = MultipleLocator(5)
majorXFormatter = FormatStrFormatter('%d')
minorXLocator = MultipleLocator(1)

for pl in range(0, 2):
	## declare figure object
	fig = pyplot.figure(figsize = (12,12))
	for cnt in range(0, 7):
		beamIndex = wvuBeamDict[np.str(cnt)]
		## determine row location span
		locSeq = rowLocDict[cnt]

		## extract stored interpolated beam pattern
		beamGridPattern = meanBeamArr[pl, cnt, :, :]

		## find index of maximum value so as to know where to take the cuts
		maxInd = np.where(beamGridPattern == np.nanmax(beamGridPattern))
		maxX = maxInd[0][0]
		maxY = maxInd[1][0]

		## extract EL/XEL cuts
		elCut = beamGridPattern[maxX, :]
		xelCut = beamGridPattern[:, maxY]

		## these are in dB... cast to normalized linear units

		elCut_Linear = 10**(elCut/10)
		xelCut_Linear = 10**(xelCut/10)

		##init guesses
		elInitInd = np.where(elCut_Linear == np.nanmax(elCut_Linear))
		elInitMean = elCut_Linear[elInitInd]
		xelInitInd = np.where(elCut_Linear == np.nanmax(xelCut_Linear))
		xelInitMean = xelCut_Linear[xelInitInd]
		sigma = 9.1
		poptEL, pcovEL = curve_fit(Gauss1D, yAxis, np.nan_to_num(elCut_Linear))#, p0 = [1, elInitMean, sigma])
		poptXEL, pcovXEL = curve_fit(Gauss1D, xAxis, np.nan_to_num(xelCut_Linear))#, p0 = [1, xelInitMean, sigma])

		## compute errors
		perrEL = np.sqrt(np.diag(pcovEL))
		perrXEL = np.sqrt(np.diag(pcovXEL))

		"""
		From the returned fit, estimate the beam area by computing the average FWHM between each XEL/EL profile
		"""

		elFWHM = poptEL[2] * 2.355 ## arcminutes
		xelFWHM = poptXEL[2] * 2.355 ## arcminutes

		## get errors
		elFWHMErr = perrEL[2]
		xelFWHMErr = perrXEL[2]

		aveFWHM = (elFWHM + xelFWHM) / 2.
		aveFWHMErr = (elFWHMErr +xelFWHMErr) / 2.

		## compute area
		beamArea = aveFWHM**2 * 1.1331 * 3600 ## arcseconds
		beamAreaErr = aveFWHMErr**2 * 1.1331 * 3600 ## arcseconds

		## inform user
		print('FWHM of EL fit: %.2f' % elFWHM + '+/-%.2f' % elFWHMErr + '[arcmin]')
		print('FWHM of XEL fit: %.2f' % xelFWHM + '+/-%.2f' % xelFWHMErr + '[arcmin]')
		print('The area of the beam estimated from Gaussian fit of profiles: %.2f' % beamArea + '+/-%2.f' % beamAreaErr + '[sq. arcseconds]')

		ax = pyplot.subplot2grid((6, 6), locSeq, colspan = 2, rowspan = 2)

		## overplot the profiles
		ax.set_title('Beam ' + beamIndex, fontsize = 12)
		ax.plot(yAxis, elCut_Linear, color = tableau20[0], label = 'EL', linewidth = 2)
		ax.plot(xAxis, xelCut_Linear, color = tableau20[2], label = 'XEL', linewidth = 2)
		ax.plot(yAxis, Gauss1D(yAxis, *poptEL), color = tableau20[0], label = 'EL (Fit)', linewidth = 2, linestyle = '--')
		ax.plot(xAxis, Gauss1D(xAxis, *poptXEL), color = tableau20[2], label = 'XEL (Fit)', linewidth = 2, linestyle = '--')
		#ax.set_ylim(0, 1.05)
				
		ax.yaxis.set_major_locator(majorYLocator)
		ax.yaxis.set_major_formatter(majorYFormatter)
		ax.yaxis.set_minor_locator(minorYLocator)
		ax.xaxis.set_major_locator(majorXLocator)
		ax.xaxis.set_major_formatter(majorXFormatter)
		ax.xaxis.set_minor_locator(minorXLocator)
		ax.tick_params(axis = 'both', which='both', width=2)

		if cnt == 3:
			ax.legend(loc=0,prop={'size':10})

	#[left, bottom, width, height]
	fig.subplots_adjust(hspace = 0.4, wspace = 0.4)
	fig.text(0.5, 0.04, 'XEL Offset [arcmin]', ha='center')
	fig.text(0.04, 0.5, 'Normalized Response', va='center', rotation='vertical')
	if pl == 0:
		pyplot.savefig(projName + '_FormedBeamProfiles_YY.png', bbox_inches='tight')
	else:
		pyplot.savefig(projName + '_FormedBeamProfiles_XX.png', bbox_inches='tight')
	pyplot.show(block=True)
	pyplot.clf()
	pyplot.close()


