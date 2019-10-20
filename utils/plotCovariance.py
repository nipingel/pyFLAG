"""
7/20/19
Script to produce time averged (i.e., averages all integartions) and covariance matrix from the 40 data channels for both polarizations. 
The resulting plot is a 40x40 pixel images, where the pixel intensity corresponds to the raw covariance in units of dB.
Inputs:
projectName - name of your GBT project and session (e.g. AGBT16B_400_14)
timeStamp - time stamp of scan

Usage:
ipython plotCorrelations.py projectName timeStamp chanNum
Example:
ipython plotCorrelations.py AGBT16B_400_14 2017_05_27_05:13:44 200
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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import glob
import sys
import pickle
from scipy.optimize import curve_fit
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
projectID = sys.argv[1] ## GBO project ID & session number in form of AGBT16B_400_02
timeStamp = sys.argv[2] ## timestamp in format 2018_01_01_00:00:00

## ERROR HANDLING
"""
First, let's check the correct number of inputs were provided
"""
if not len(sys.argv[1:]) == 2:
    print('Incorrect number of user provided inputs; usage: ipython PlotCorrelations.py PROJECT_SESSION TIMESTAMP; Exiting...')
    sys.exit(1)

"""
make sure the project string is valid. Any ancillary FITS files will be in /home/gbtdata/projectID
"""    
if not os.path.exists('/home/gbtdata/' + projectID):
    print('Project ID and session does not exist. Exiting...')
    sys.exit(1)

"""
check time stamp format and if it exists
"""
if not (timeStamp[4] == '_') and not (timeStamp[7] == '_') and not (timeStamp[-6] == ':') and not (timeStamp[-3] == ':'):
    print('Incorrect format of the time stamp; should be YYYY_MM_DD_HH:MM:SS')
    sys.exit(1)

if not os.path.isfile('/home/gbtdata/' + projectID + '/Antenna/' + timeStamp + '.fits'):
    print('This is not a recorded scan; exiting...')
    sys.exit(1)

## progress bar
def progressBar(value, endvalue,bar_length=20):

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent of integrations processed: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()

"""
Method to determine which covariance BANK file to read in
"""
def getDataVector(fitsList):
    ## determine which mode this scan was taken
    hdu = fits.open(fitsList[0])
    modeName = hdu[0].header['MODENAME']

    """
    determine which bank to open based on channel number. Channel order is mode dependent. 
    If in calcorr mode, bank file with XID 0 contains coarse channels 0-4, 100-104, ..., 400-404
    XID 1: 5-9, 105-109, ...404, ... XID 19: 95-99, 195-199, ... , 495-499

    Select the channel associated with HI line in mode 
    """
    if modeName == 'FLAG_CALCORR_MODE':
        retChanNum = 1
        hdu = fits.open(fitsList[10])
        dataArr = hdu[1].data['DATA']
        return dataArr, retChanNum
    elif modeName == 'FLAG_PFBCORR_MODE':
        retChanNum = 1600
        hdu = fits.open(fitsList[10])
        dataArr = hdu[1].data['DATA']
        return dataArr, retChanNum
      
## define data directories
dataDir = '/lustre/flag/' + projectID + '/BF/'

fitsList = glob.glob(dataDir + timeStamp+'*.fits')
fitsList.sort()

## check to see if the FITS files were found
if (len(fitsList) == 0):
    print('No data files were found; check dataDir variable')
    sys.exit(1)

"""
Read in gpuToNative_Map.dat in ../misc directory. The first element of this 1D vector corresponds to the index of the 
1x1 correlation, the second element is the index for 2x1, the third for 3x1 etc...
"""
mapVector = np.loadtxt('../misc/gpuToNativeMap.dat', dtype='int')

## get data based on channel number
dataArr, retChanNum = getDataVector(fitsList)
numInts = len(dataArr[:, 0])
numElem = 40
corrCube = np.zeros([numElem, numElem, numInts], dtype=np.complex64)
for i in range(0, numInts):
    progressBar(i, numInts)
    singleChanData = dataArr[i, retChanNum*2112:retChanNum*2112 + 2112] ## grab a single freq. channel worth of correlations
    newDataVector= np.zeros([820], dtype='complex64')
    """
    loop through and assign new order based on the mapVector attribute, which gives the new index of the 
    raw correlation
    """ 
    for z in range(0,len(mapVector)):
        newDataVector[z] = singleChanData[mapVector[z]]

    ## Once sorted, we can construct a retCube with dim1 = numElements, dim2 = numElements, dim3 = freqChans
    cnt = 0

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
            corrCube[row, col, i] = newDataVector[cnt]
            if row != col :
               corrCube[col,row, i] = newDataVector[cnt].conj()
            cnt += 1
## now, average along the z axis
corrMatrix = np.mean(corrCube, axis = 2)

## set plotting parameters
np.seterr(divide = 'ignore')
majorYLocFactor = 5
minorYLocFactor = 1
majorXLocator = 5
minorXLocator = 1

majorYLocator = MultipleLocator(majorYLocFactor)
majorYFormatter = FormatStrFormatter('%d')
minorYLocator = MultipleLocator(minorYLocFactor)

majorXLocator = MultipleLocator(majorXLocator)
majorXFormatter = FormatStrFormatter('%d')
minorXLocator = MultipleLocator(minorXLocator)

fig, ax = pyplot.subplots()
im = ax.imshow(10*np.log10(abs(corrMatrix)**2/np.max(abs(corrMatrix)**2)))
ax.set_title('Time-Averaged Correlations: %s' % timeStamp, size = 12)
cb = fig.colorbar(im, label = '[dB]')
cb.ax.tick_params(which = 'minor', length = 2)
cb.ax.tick_params(which = 'major', length = 4)
cb.ax.minorticks_on()
ax.yaxis.set_major_locator(majorYLocator)
ax.yaxis.set_major_formatter(majorYFormatter)
ax.yaxis.set_minor_locator(minorYLocator)
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
ax.xaxis.set_minor_locator(minorXLocator)
ax.set_xlabel('Data Channel')
ax.set_ylabel('Data Channel')
ax.tick_params(axis = 'both', which='both', width=2)
pyplot.show(block = True)
pyplot.clf()
pyplot.close()





