"""
01/16/18
Script to plot and compute average Tsys/eta over FLAG bandpass. An two-column (average Tsys value and 
associated uncertainty) txt file is created. Each row corresponds to one formed beam.
Inputs:
/path/to/mat/file/Xpol - that contains the Tsys/eta for X Pol as a function frequency
/path/to/mat/file/Ypol - that contains the Tsys/eta for Y Pol as a function frequenc)
calType - derived from a full calibation grid or just a seven point calibration (Grid or 7-Pt Cal)
plotFlag - plot and save figures for user. Either True or False is accepted.
Usage:
ipython plotTsys.py /path/to/mat/file/Xpol /path/to/mat/file/Ypol plotFlag
Example:
ipython plotTsys.py /lustre/projects/flag/AGBT16B_400_14/BF/mat/tSys_x.mat /lustre/projects/flag/AGBT16B_400_14/BF/mat/tSys_y.mat Grid True
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "nipingel@mix.wvu.edu"
__status__ = "production"
"""

## imports
import matplotlib.pyplot as pyplot
import numpy as np

import scipy.io
import sys
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

majorYLocator = MultipleLocator(10)
majorYFormatter = FormatStrFormatter('%d')
minorYLocator = MultipleLocator(2)

majorXLocator = MultipleLocator(20)
majorXFormatter = FormatStrFormatter('%d')
minorXLocator = MultipleLocator(5)

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

## make beam dictionary to map from BYU to WVU convention (e.g. BYU 0 -> WVU 1)
wvuBeamDict = {'0':'1', '1':'2', '2':'6', '3':'0', '4':'3', '5':'5','6':'4'}
byuBeamDict = {'1':'0', '2':'1', '6':'2', '0':'3', '3':'4', '5':'5','4':'6'}

## test to ensure correct number of user inputs are given
if len(sys.argv) != 5:
        print('Incorrect number of command line arguments; usage -> python plotTsys.py /path/to/mat/file/XPol /path/to/mat/file/YPol calType plotFlag')
        sys.exit(1)
dataPathX = sys.argv[1]
dataPathY = sys.argv[2]
calType = sys.argv[3]
plotFlag = sys.argv[4]



## test that plotFlag is correctly entered as True or False
if plotFlag != 'True' and plotFlag != 'False':
        print('plotFlag variable is not correctly defined as True or False')
        sys.exit(1)


## open matlabfiles
## try-catch in case file is not found or file is not matlab binary file
try:
        matX = scipy.io.loadmat(dataPathX)
        matY = scipy.io.loadmat(dataPathY)
except ValueError:
        print('File is not a matlab binary file. Exiting...')
        sys.exit(1)

## after error handling, extract project name
pathList = dataPathX.split('/')
projectID = pathList[-3] ## third last (e.g. projectID/projectSession/Tsys/Files)

    ## get beam indexes if calType is grid
if calType == 'grid' and projectID == 'AGBT17B_360_01':
    bmInds = [324, 309, 245, 248, 250, 204, 189]
elif calType == 'grid':
    bmInds = [669, 664, 507, 512, 517, 385, 390]
#elif calType == 'grid' and projectID[:-3] == 'AGBT16B_400':
    #print('here')
    #bmInds = [865, 558, 318, 251, 211, 797, 898]
#elif calType == 'grid' and projectID == 'AGBT16B_400_12':
#    bmInds = [1246, 1256, 892, 882, 872, 559, 549]



freqAxis = matX['freq']
Tsys_eta_X = np.array(matX['Tsys_eta'])
Tsys_eta_Y = np.array(matY['Tsys_eta'])


## numBeamsx500 (beamsxfreqChans); make plot and calculate Tsys_X/Y 
numBeams = 7
## assuming eta=0.60 (Roshi et al. 2017)
eta = 0.60
## set freq chan ranges for first half of bandpass (1400-1416.58 MHz assuming 1450 as center frequency)
startChan1 = 82
endChan1 = 140
## set freq chan range for second half of bandpass (1425.0925 to 1440.28 MHz assuming 1450 as center frequency)
startChan2 = 155
endChan2 = 218

## initialize arrays to store meanTsys and uncertainty values
tSysSigArray_X = np.zeros([numBeams, 2])
tSysSigArray_Y = np.zeros([numBeams, 2])


for i in range(0, 7):
    ## if calType is grid, get correct beam index
    if calType == 'grid':
        bmIndex = bmInds[i]
    else:
        bmIndex = i
	## calculate statistics
	## take mean between 1400.185 and 1416.5875 MHz
    Tsys_eta_X[Tsys_eta_X == np.inf] = np.nan
    Tsys_eta_Y[Tsys_eta_Y == np.inf] = np.nan

    meanTsys_X_1 = np.nanmean(Tsys_eta_X[bmIndex, startChan1:endChan1]) * eta
    ## take mean between 1425.0925 and 1440.28 MHz
    meanTsys_X_2 = np.nanmean(Tsys_eta_X[bmIndex,startChan2:endChan2]) * eta
    meanTsys_X = (meanTsys_X_1 + meanTsys_X_2)/2
    ## calculate standard deviation of two distributions
    stdTsys_X_1 = np.nanstd(Tsys_eta_X[bmIndex, startChan1:endChan1]) * eta
    stdTsys_X_2 = np.nanstd(Tsys_eta_X[bmIndex, startChan2:endChan2]) * eta
    stdTsys_X = np.sqrt(stdTsys_X_1**2 + stdTsys_X_1**2)/2
    meanTsys_Y_1 = np.nanmean(Tsys_eta_Y[bmIndex, startChan1:endChan2]) * eta
    meanTsys_Y_2 = np.nanmean(Tsys_eta_Y[bmIndex, startChan2:endChan2]) * eta
    meanTsys_Y = (meanTsys_Y_1 + meanTsys_Y_2)/2
    stdTsys_Y_1 = np.nanstd(Tsys_eta_Y[bmIndex, startChan1:endChan1]) * eta
    stdTsys_Y_2 = np.nanstd(Tsys_eta_Y[bmIndex, startChan2:endChan2]) * eta
    stdTsys_Y = np.sqrt(stdTsys_Y_1**2 + stdTsys_Y_1**2)/2

    ## make beam dictionary to map from BYU to WVU convention (e.g. BYU 0 -> WVU 1)
    beamDict = {0:1, 1:2, 2:6, 3:0, 4:3, 5:5 ,6:4}

    ## place in arrays
    tSysSigArray_X[beamDict[i], 0] = meanTsys_X
    tSysSigArray_X[beamDict[i], 1]= stdTsys_X
    tSysSigArray_Y[beamDict[i], 0] = meanTsys_Y
    tSysSigArray_Y[beamDict[i], 1] = stdTsys_Y

    ## make beam dictionary to map from BYU to WVU convention (e.g. BYU 0 -> WVU 1)
    beamDict = {0:1, 1:2, 2:6, 3:0, 4:3, 5:5 ,6:4}

    ## if plotFlag set as true, then plot tsys/eta for each pol
    if plotFlag == 'True':
        fig, ax = pyplot.subplots()
        pyplot.xlim(1390, 1500)
        pyplot.ylim(meanTsys_X + 10, meanTsys_X +30)
        pyplot.axvline(1400.185, linestyle = '--', color = 'black')
	pyplot.axvline(1416.5875, linestyle = '--', color = 'black')
	pyplot.axvline(1425.0925, linestyle = '--', color = 'black')
        pyplot.axvline(1440.28, linestyle = '--', color = 'black')	
	pyplot.xlabel('Frequency [MHz]')
        pyplot.ylabel(r'Tsys/$\eta$ [K]')
        pyplot.title('Beam: ' + np.str(beamDict[i]) + ' (' + projectID +' ' + calType + ')')
        pyplot.plot(freqAxis[0], Tsys_eta_X[bmIndex,:], label = 'XX-Pol', color = tableau20[2], linewidth=2)
        pyplot.plot(freqAxis[0], Tsys_eta_Y[bmIndex,:], label = 'YY-Pol', color = tableau20[0], linewidth=2)
        pyplot.legend(loc=2, fontsize=14)

        ## format tick marks
        ax.yaxis.set_major_locator(majorYLocator)
        ax.yaxis.set_major_formatter(majorYFormatter)
        ax.yaxis.set_minor_locator(minorYLocator)
        ax.xaxis.set_major_locator(majorXLocator)
        ax.xaxis.set_major_formatter(majorXFormatter)
        ax.xaxis.set_minor_locator(minorXLocator)
        ax.tick_params(axis = 'both', which='both', width=2)
        pyplot.savefig('Tsys_Eta_' + projectID + '_Bm' + np.str(i) + '_' + calType + '.pdf')
        pyplot.show()
        pyplot.clf()
        pyplot.close()


## print results
for j in range(0, numBeams):
    print('\n')
    print('Beam  (WVU): ' + np.str(j))
    print('Mean Tsys_X: ' + np.str(tSysSigArray_X[j, 0]) + '+/-' + np.str(tSysSigArray_X[j, 1]))
    print('Mean Tsys_Y: ' + np.str(tSysSigArray_Y[j, 0]) + '+/-' + np.str(tSysSigArray_Y[j, 1]))

## save out to two separate text files for each pol
## text files will be two columned: (1) Mean Tsys  (2) Standard Devation. 
## each row represents these values for a given beam

np.savetxt('tSys_X_' + projectID + '_' + calType + '.txt', tSysSigArray_X)
np.savetxt('tSys_Y_' + projectID + '_' + calType + '.txt', tSysSigArray_Y)


