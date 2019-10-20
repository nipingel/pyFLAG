"""
8/3/18
Script to correct Antenna positions of PAF observations. Based on the input file, weight files, and object, the associated Antenna 
and data FITS files will be read in, the Antenna position interpolated to the data (using the respective DMJD values), 
and beam offsets applied. Once the beam offsets are applied, the horiztonal coordinates are transformed to J2000 or Galactic
coordinates. A doppler correction is additionally performed. The user also needs to specify the rest and center frequency. 
User Inputs:
fileName - path to input rawbeamformed SDFITS file
weightPath - path to directory that holds weight FITS files
observedObj - the observed object
restFreq - rest frequency of observed linewidth [Hz]
cenFreq - center TOPOCENTRIC frequency of band [Hz]
Usage:
ipython FixAntPositions.py fileName observedObj weightPath restFreq cenFreq
Eample:
ipython FixAntPositions.py AGBT17B_455_01_G353-4.0_Beam0_1st_seven.fits /home/weights/*.FITS G353-4.0 1420.406e6 1450.0e6
__email__ = "nipingel@mix.wvu.edu"
__status__ = "Production"
"""

from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as pyplot
import pyslalib as pysla
import os 
import numpy as np
import datetime
import collections
import glob
import sys
sys.path.insert(0, '/home/scratch/npingel/FLAG/pros/SpectralFiller/modules/')
import RadVelCorr

"""
Function that generates and returns list of observed objects and the associated timestamp FITS files. 
The GO FITS files are searched and sorted to generate the required lists, which is why the path to these 
files is an input to this function. 
"""
def generateObjAndFitsLists(goFitsPath):        

    goFits = glob.glob(goFitsPath+'/*.fits') ## read in GO FITS files to
    goFits.sort() ## sort to get correct time stamp order
    genObjList = [] ## define empty list to be returned (after casting to list for indexing)
    """
    define empty list. Once we add new object, an empty list will be appended and
    associated FITS files added. For example, if the second object added is '3C147', 
    genFitsList[1] will be a list of all associated timestamp FITS files. 
    """
    genFitsList = [] 
    itr = 0 ## counter variable to check if first iteration
    ## iterate over all GO FITS files to check the observed object
    for goName in goFits:
       goHDU = fits.open(goName)
       if itr == 0: ## upon first iteration, append an empty list to the generated list of FITS files
           genObjList.append(goHDU[0].header['OBJECT'])
           genFitsList.append([])
           genFitsList[0].append(goName[-24:])
           itr+=1
       else: ## now processing for remaining 
           obj = goHDU[0].header['OBJECT'] 
           if not any(x == obj for x in genObjList): ## test if object already in list
               genObjList.append(obj) 
               genFitsList.append([]) ## append empty list
               genFitsList[-1].append(goName[-24:])
           else:
            ind = genObjList.index(obj)
            genFitsList[ind].append(goName[-24:])
    return genObjList, genFitsList

fileName = sys.argv[1] ## get input file name path
weightPath = sys.argv[2] ## get weight directory path
observedObj = sys.argv[3] ## get observed object
restFreq = np.float(sys.argv[4]) 
cenFreq = np.float(sys.argv[5])

hdu = fits.open(fileName, mode = 'update') ## read in FITS file
projectID = fileName[:14] ## get project name
weightFileList = glob.glob(weightPath)
weightHdu = fits.open(weightFileList[0]) ## get weight file


## read in ancilllary FITS files
scanFitsHdu = fits.open('/home/archive/science-data/17B/' + projectID + '/ScanLog.fits')
goFitsPath = '/home/archive/science-data/17B/' + projectID + '/GO'

## generate object and list of FITS files in case user wishes to process all object/timestamps
allObjList, allFitsList = generateObjAndFitsLists(goFitsPath)

objInd = allObjList.index(observedObj)
source = allObjList[objInd] ## get source
fileList = allFitsList[objInd] ## get FITS list for object (just fileneams; no path)

## define lists that will hold data coordinates
indicatedCoordList = []
majList = []
minList = []
raList = []
decList = []
azList = []
elList = []
lstStartList = []
dataDMJDList = []
refractList = []
velList = []
polMotionXList = []
polMotionYList = []
utCorrList = []
tempList = []
pressList = []
humList = []
waveLengthList = []
lapseRateList = []

GBTLAT  = np.deg2rad(38.4331294)
GBTLONG = np.deg2rad(280.160200)
GBTHGT  = 824.36                     # meters above the ellipsoid
c = 299792458.0 ## m/s

#restFreq = 1.4204057517667e9 
#cenFreq = 1449.84841e6
#cenFreq = 1450e6
chanSel = 1
lowEnd = cenFreq - (250 - (chanSel*100))*0.30318*1e6
highEnd = cenFreq - (250 - (chanSel*100 + 100))*0.30318*1e6
freqCrVal = highEnd - (highEnd - lowEnd)/2

def setFreq(fitsHdu):
	fitsHdu[1].data['RESTFREQ'] = restFreq
	fitsHdu[1].data['OBSFREQ'] = freqCrVal
	fitsHdu[1].data['CRVAL1'] = freqCrVal
	return fitsHdu
	

def radVelCorrection(raArr, decArr):

	print('\n')
	print('Applying Doppler Correction to HELIOCENTRIC...')

	radvelcorrObj = RadVelCorr.RadVelCorr() ## initialize object
	
	newHdu = setFreq(hdu) ## ensure center frequency values are correctly set
	
	## collect all relevant data for doppler calculation
	utDate_Time = newHdu[1].data['DATE-OBS']
	cenFreqsArr = newHdu[1].data['CRVAL1']
	restFreqArr = newHdu[1].data['RESTFREQ']
	for velIter in range(0,len(raArr)): ## loop through to calculate each integration/polarization's correciton 
	  raVal = raArr[velIter] ## RA
	  decVal = decArr[velIter] ## Dec 
	  utDateTimeVal = utDate_Time[velIter] ## UT Time
	  cenFreqVal = cenFreqsArr[velIter] ## center frequency value
	  restFreqVal = restFreqArr[velIter] ## restfreq
	  ## put time in necessary format
	  t = Time(utDateTimeVal, format='isot')
	  utdate = t.iso[0:10]
	  uttime = t.iso[11:]
	  radVelCorr_HEL, radVelCorr_LSR = radvelcorrObj.correctVel(utdate,uttime,raVal,decVal) ## calculate correction
	  ## compute optical velocity of ref freq
	  vOpt = c*(1-cenFreqVal/restFreqVal)
	  ## add radial correction
	  newVOpt = vOpt + radVelCorr_HEL
	  ## now convert back to frequency
	  newCenFreq = (1-newVOpt/c)*restFreqVal 
	  ## update reference frequency value
	  cenFreqsArr[velIter] = newCenFreq

	return cenFreqsArr

def az2ra(LST, Az, El, dec, coordSys):
	## convert all angles to radians
	azRad = np.deg2rad(Az) ## Azimuth of NCP is 0 for GBT. This code assumes 180. 
	posInds = np.where(Az > 360)
	negInds = np.where(np.logical_and(Az < 180, Az > 0)) ## need to rotate multiply by negative 1 while in 
	elRad = np.deg2rad(El)
	decRad = np.deg2rad(dec)
	ha = (-1) * np.arccos(1/np.cos(decRad)*(np.sin(elRad)*np.cos(GBTLAT) - np.cos(elRad)*np.cos(azRad)*np.sin(GBTLAT)))
	ha[posInds] = ha[posInds] * (-1)
	ha[negInds] = ha[negInds] * (-1)
	return LST + np.rad2deg(ha)
 
"""  
method to convert Elevation to Declination (geoapparent)
""" 
def el2dec(LST, Az, El):
	## convert all angles to radians
	azRad = np.deg2rad(Az) ## Azimuth of NCP is 0 for GBT. This code assumes 180. 
	elRad = np.deg2rad(El)
	el = np.rad2deg(np.arcsin(np.sin(elRad)*np.sin(GBTLAT) + np.cos(elRad)*np.cos(azRad)*np.cos(GBTLAT)))
	return el

"""
method to transform the beam offsets to Ra/Dec coordinates. The current state of the binary table HDU is passed in. 
That same binary table HDU with the updated coordinates is returned. 
"""
def offsetCorrection(extCoordSys, majList, minList, azList, elList, raList, decList, dataDMJDList, refractList, lstStartSec, velList, polMotionXList, polMotionYList, utCorrList, tempList, pressList, humList, waveLengthList, lapseRateList):
	print('\n')
	print('Correcting coordinates for beam offset...')
	## we have refract value for a single pol
	## extend by 2x to describe both XX and YY pol
	extCoordSys.extend(extCoordSys)
	newMajArr = np.zeros(len(refractList)*2)
	newMinArr = np.zeros(len(refractList)*2)
	newAzArr = np.zeros(len(refractList)*2)
	newElArr = np.zeros(len(refractList)*2)
	newRaArr = np.zeros(len(refractList)*2)
	newDecArr = np.zeros(len(refractList)*2)
	extRefractArr = np.zeros(len(refractList)*2)
	extMajArr = np.zeros(len(refractList)*2)
	extMinArr = np.zeros(len(refractList)*2)
	extAzArr = np.zeros(len(refractList)*2)
	extElArr = np.zeros(len(refractList)*2)
	extRaArr = np.zeros(len(refractList)*2)
	extDecArr = np.zeros(len(refractList)*2)
	extDataDMJDArr = np.zeros(len(refractList)*2)
	extLstArr = np.zeros(len(refractList)*2)
	dataDMJDArr = np.zeros(len(refractList)*2)
	velArr = np.zeros(len(refractList)*2)
	beamXElArr = np.deg2rad(hdu[1].data['FEEDXOFF'])
	beamElArr = np.deg2rad(hdu[1].data['FEEDEOFF'])

	polMotionXArr = np.zeros(len(refractList)*2)
	polMotionYArr = np.zeros(len(refractList)*2)
	utCorrArr = np.zeros(len(refractList)*2)
	tempArr = np.zeros(len(refractList)*2)
	pressArr = np.zeros(len(refractList)*2)
	humArr = np.zeros(len(refractList)*2)
	waveLengthArr = np.zeros(len(refractList)*2)
	lapseRateArr = np.zeros(len(refractList)*2)

	extRefractArr[0::2] = refractList
	extRefractArr[1::2] = refractList

	extMajArr[0::2] = majList
	extMajArr[1::2] = majList

	extMinArr[0::2] = minList
	extMinArr[1::2] = minList

	extAzArr[0::2] = azList
	extAzArr[1::2] = azList

	extElArr[0::2] = elList
	extElArr[1::2] = elList

	extRaArr[0::2] = raList
	extRaArr[1::2] = raList
	extDecArr[0::2] = decList
	extDecArr[1::2] = decList

	extLstArr[0::2] = lstStartSec
	extLstArr[1::2] = lstStartSec

	lstStartDegs = extLstArr/3600*15 ## convert to degs

	dataDMJDArr[0::2] = dataDMJDList
	dataDMJDArr[1::2] = dataDMJDList

	velArr[0::2] = velList
	velArr[1::2] = velList
	polMotionXArr[0::2] = polMotionXList
	polMotionXArr[1::2] = polMotionXList
	polMotionYArr[0::2] = polMotionYList
	polMotionYArr[1::2] = polMotionYList
	utCorrArr[0::2] = utCorrList
	utCorrArr[1::2] = utCorrList
	tempArr[0::2] = tempList
	tempArr[1::2] = tempList
	pressArr[0::2] = pressList
	pressArr[1::2] = pressList
	humArr[0::2] = humList
	humArr[1::2] = humList
	waveLengthArr[0::2] = waveLengthList
	waveLengthArr[1::2] = waveLengthList
	lapseRateArr[0::2] = lapseRateList
	lapseRateArr[1::2] = lapseRateList

	## loop through each coordinate to precess from J2000 (mean place) to geoapparent RA/Dec
	for coordIdx in range(0, len(extRaArr)):
		raVal = extRaArr[coordIdx]
		decVal = extDecArr[coordIdx]
		velVal = velArr[coordIdx]
		dmjdVal = dataDMJDArr[coordIdx]

		deltaUT = utCorrArr[coordIdx]
		polXVal = polMotionXArr[coordIdx]
		polYVal = polMotionYArr[coordIdx]
		tempVal = tempArr[coordIdx]
		pressVal = pressArr[coordIdx]
		humVal = humArr[coordIdx]
		waveVal = waveLengthArr[coordIdx]
		lapseVal = lapseRateArr[coordIdx]


		## convert J2000 -> geoapparent (center of Earth)
		geoRa, geoDec = pysla.slalib.sla_map(np.deg2rad(raVal), np.deg2rad(decVal), 0.0, 0.0, 0.0, velVal, 2000.0, dmjdVal)
		## convert from center-of-earth to observed at GBO
		obsAz, obsZen, obsHA, obsDec, obsRa = pysla.slalib.sla_aop(geoRa, geoDec, dmjdVal, deltaUT, GBTLONG, GBTLAT, GBTHGT, polXVal, polYVal, tempVal, pressVal, humVal, waveVal, lapseVal)
		## convert from obs eq to horizontal
		obsAz, obsEl = pysla.slalib.sla_e2h(obsHA, obsDec, GBTLAT)
		## apply refraction correction

		## apply beam offset corrections
		newElVal = obsEl - beamElArr[coordIdx] + np.deg2rad(extRefractArr[coordIdx])
		newAzVal = obsAz - (beamXElArr[coordIdx] / np.cos(newElVal))

		## place in arrays
		extAzArr[coordIdx] = np.rad2deg(obsAz)
		extElArr[coordIdx] = np.rad2deg(obsEl)
		newAzArr[coordIdx] = np.rad2deg(newAzVal)
		newElArr[coordIdx] = np.rad2deg(newElVal)

		## put the major/minor array in the correct coordinate system based on INDICSYS
		if extCoordSys[coordIdx] == 'AZEL':
			# major and minor are az and el
			newMajArr[coordIdx] = np.rad2deg(newAzVal)
			newMinArr[coordIdx] = np.rad2deg(newElVal)
		else:

			## convert from obs horiz to eq 
			newObsHA, newObsDec = pysla.slalib.sla_h2e(newAzVal, newElVal, GBTLAT)
			newObsRa = np.deg2rad(lstStartDegs[coordIdx]) - newObsHA
			## convert from observed at GBO to geocentric (center of Earth)
			newGeoRa, newGeoDec = pysla.slalib.sla_oap('R', newObsRa, newObsDec, dmjdVal, deltaUT, GBTLONG, GBTLAT, GBTHGT, polXVal, polYVal, tempVal, pressVal, humVal, waveVal, lapseVal)
		
			## convert from geocenteric to mean
			newRa, newDec = pysla.slalib.sla_amp(newGeoRa, newGeoDec, dmjdVal, 2000.0)

			## convert radians to degrees
			newRaArr[coordIdx] = np.rad2deg(newRa)
			newDecArr[coordIdx] = np.rad2deg(newDec)

			newMajArr[coordIdx] = np.rad2deg(newRa)
			newMinArr[coordIdx] = np.rad2deg(newDec)


			if extCoordSys[coordIdx] == 'GALACTIC':
				glon, glat = pysla.slalib.sla_eqgal(newRa, newDec)
				newMajArr[coordIdx] = glon
				newMinArr[coordIdx] = glat
		
	return newMajArr, newMinArr, extMajArr, extMinArr, newAzArr, newElArr, extAzArr, extElArr, newRaArr, newDecArr, extRaArr, extDecArr, lstStartDegs
"""
begin looping over files to get Antenna and data positions, DMJD values, and perfrom correct coord
"""
fileList = fileList[1:]
for flName in fileList:
	print(flName)
	antHdu = fits.open('/home/archive/science-data/17B/' + projectID + '/Antenna/' + flName) 
	dataHdu = fits.open('/lustre/flag/' + projectID + '/BF/' + flName[:-5] + 'A.fits') ## only need one correlator bank FITS file
	goHDU = fits.open(goFitsPath + '/' + flName)

	## extract pointing model information 
	smntAz = antHdu[0].header['SMNTC_AZ']
	sobscAz = antHdu[0].header['SOBSC_AZ']
	sobscEl = antHdu[0].header['SOBSC_EL']
	smntEl = antHdu[0].header['SMNTC_EL']

	indictedCoords = antHdu[0].header['INDICSYS']
	  
	## calculate pointing model
	azPt = smntAz - sobscAz
	elPt = smntEl - sobscEl


	## get DMJD values
	antDMJD = antHdu[2].data['DMJD']
	corrDMJD = dataHdu[1].data['DMJD']

	## get start of scan and convert to mjd
	startTimeArr = antHdu[0].header['DATE-OBS']
	startObj = Time(startTimeArr, format = 'isot', scale = 'utc')
	startMJD = startObj.mjd

	##get Antenna MAJOR and MINOR axis and refraction values for interpolation
	antMaj = antHdu[2].data['MAJOR']
	antMin = antHdu[2].data['MINOR']

	antRa = antHdu[2].data['RAJ2000']
	antDec = antHdu[2].data['DECJ2000']

	## get refraction correction values
	antRefract = antHdu[2].data['REFRACT']
	refractInterp = np.interp(corrDMJD, antDMJD, antRefract)

	## get Azimuth/Elevation and subtract off pointing model
	az = antHdu[2].data['MNT_AZ']# - azPt
	el = antHdu[2].data['MNT_EL']# - elPt
	#az = antHdu[2].data['OBSC_AZ']# - azPt
	#el = antHdu[2].data['OBSC_EL']# - elPt
	#az = antHdu[2].data['MNT_AZ']- antHdu[2].data['OBSC_AZ'] - azPt
	#el = antHdu[2].data['MNT_EL']- antHdu[2].data['OBSC_EL'] - elPt
	azInterp = np.interp(corrDMJD, antDMJD, az)
	elInterp = np.interp(corrDMJD, antDMJD, el)

	lstStart = antHdu[0].header['LSTSTART']
	intLen = np.float(dataHdu[0].header['ACTSTI'])
	lstStartArr = np.zeros([len(azInterp)])
  	lstStartArr = (1.00273790935 * (corrDMJD - startMJD)*86400.0) + lstStart

  	## create list the same length as the coordinate lists and fill with source radial velocity value, polar motions, 
  	## UT1-UTC correction, temperature, pressure, humidity, wavelength, and toposperifc lapse rate; basically, all things
  	## neccessary to go from mean J2000 -> geocentric J2000 -> observered (eq) -> observed (horiz) -> apply offset -> 
  	## modified bserved (horiz) -> modified observed (eq) -> modifed geocentric J2000 -> modified mean J2000
  	velVal = goHDU[0].header['VELOCITY']/1000. ## km/s
  	polMotionX = antHdu[0].header['IERSPMX'] ## radians
  	polMotionY = antHdu[0].header['IERSPMY'] ## radians
  	utCorr = antHdu[0].header['DELTAUTC'] ## s
  	temp = antHdu[0].header['AMBTEMP'] + 273.15 ## K
  	press = antHdu[0].header['AMBPRESS'] ## millibars
  	humidity = antHdu[0].header['AMBHUMID'] ## relative fraction
  	effWaveLength = c / cenFreq/1e6 ## micrometers



  	tempVelList = [velVal] * len(azInterp)
  	tempPolXList = [polMotionX] * len(azInterp)
  	tempPolYList = [polMotionY] * len(azInterp)
  	tempUtCorrList = [utCorr] * len(azInterp)
  	tempTempList = [temp] * len(azInterp)
  	tempPressList = [press] * len(azInterp)
  	tempHumList = [humidity] * len(azInterp)
  	tempEffWaveList = [effWaveLength] * len(azInterp)
  	tempLapseRateList = [0.0065] * len(azInterp)
  	tempIndCoordList = [indictedCoords] * len(azInterp)
	
	## interpolate antenna values to data timestamps
	corrMaj = np.interp(corrDMJD, antDMJD, antMaj)
	corrMin = np.interp(corrDMJD, antDMJD, antMin)

	## interpolate atnenna values to data timestamps (RA/Dec)
	raInterp = np.interp(corrDMJD, antDMJD, antRa)
	decInterp = np.interp(corrDMJD, antDMJD, antDec)

	## place in lists
	indicatedCoordList.extend(tempIndCoordList)
	majList.extend(corrMaj)
	minList.extend(corrMin)
	raList.extend(raInterp)
	decList.extend(decInterp)
	azList.extend(azInterp)
	elList.extend(elInterp)
	lstStartList.extend(lstStartArr)
	dataDMJDList.extend(corrDMJD)
	refractList.extend(refractInterp)
	velList.extend(tempVelList)
	polMotionXList.extend(tempPolXList)
	polMotionYList.extend(tempPolYList)
	utCorrList.extend(tempUtCorrList)
	tempList.extend(tempTempList)
	pressList.extend(tempPressList)
	humList.extend(tempHumList)
	waveLengthList.extend(tempEffWaveList)
	lapseRateList.extend(tempLapseRateList)



newMajArr, newMinArr, dataMajArr, dataMinArr, newAzArr, newElArr, dataAzArr, dataElArr, newRaArr, newDecArr, dataRaArr, dataDecArr, lstStartArr = offsetCorrection(indicatedCoordList, majList, minList, azList, elList, raList, decList, dataDMJDList, refractList, lstStartList, velList, polMotionXList, polMotionYList, utCorrList, tempList, pressList, humList, waveLengthList, lapseRateList)
## inform user about mean offset
print('Mean difference between offset and data RA [arcmin]: %s' % str(np.mean(newRaArr - dataRaArr)*60.))
print('Mean difference between offset and data DEC [arcmin]: %s' % str(np.mean(newDecArr - dataDecArr)*60.))

pyplot.figure()
pyplot.plot(newMajArr, label = 'Computed', linewidth = 2)
pyplot.plot(dataMajArr, label = 'Data', linewidth = 2)
pyplot.xlabel('Element')
pyplot.ylabel('RA (J2000) [deg]')
pyplot.title(fileName[:-5] + ' (' + observedObj + ')' + ' RA Offset Correction', fontsize = 10)
pyplot.legend(loc=0, fontsize=14)
pyplot.savefig(fileName[:-5] + '_' + observedObj + '_Ra_Corr.pdf')
pyplot.show()

pyplot.figure()
pyplot.plot(newMinArr, label = 'Computed', linewidth = 2)
pyplot.plot(dataMinArr, label = 'Data', linewidth = 2)
pyplot.title(fileName[:-5] + ' (' + observedObj + ')' + ' Dec Offset Correction', fontsize = 10)
pyplot.xlabel('Element')
pyplot.ylabel('Dec (J2000) [deg]')
pyplot.legend(loc=0, fontsize=14)
pyplot.savefig(fileName[:-5] + '_' + observedObj + '_Dec_Corr.pdf')
pyplot.show()

"""
pyplot.figure()
#pyplot.scatter(newMajArr, newMinArr, label = 'Computed', linewidth = 2, color='red')
pyplot.xlabel('Az')
pyplot.ylabel('El')
#pyplot.legend(loc=0, fontsize=14)
#pyplot.savefig(fileName[:-5] + '_' + observedObj + '_Dec_Corr.pdf')
pyplot.plot(dataAzArr * np.cos(np.deg2rad(dataElArr)), dataElArr, 'bo', label = 'Computed', linewidth = 2)
#pyplot.legend(loc=0, fontsize=14)
#pyplot.savefig(fileName[:-5] + '_' + observedObj + '_Dec_Corr.pdf')
pyplot.show()
"""

## update values for spatial coordinates
hdu[1].data['CRVAL2'] = newMajArr
hdu[1].data['CRVAL3'] = newMinArr
hdu[1].data['TRGTLONG'] = newMinArr
hdu[1].data['TRGTLAT'] = newMajArr
hdu[1].data['AZIMUTH'] = newAzArr
hdu[1].data['ELEVATIO'] = newElArr

#hdu[1].data['RAJ2000'] = newRaArr
#hdu[1].data['DECJ2000'] = newDecArr

cenFreqsArr = radVelCorrection(newRaArr, newDecArr)

## update values for spectral coordinates
hdu[1].data['CRVAL1'] = cenFreqsArr
hdu[1].data['OBSFREQ'] = cenFreqsArr
hdu[1].data['VELDEF'] = 'OPTI-HEL'

hdu.flush()


