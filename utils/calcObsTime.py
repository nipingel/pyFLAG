"""
7/26/19
This script provides a sensitivity calculator and mapping specific to the FLAG recievier on the GBT. 
It will report the total observing time (excluding overhead) and assumes an integration time of 0.5 [s], SEFD of 10 [Jy], 
and gain of 1.86 K/Jy (Pingel et al. 2018). No input validity checks are performed. 
User inputs:
1 - l_map = length of map in longitude [deg]
2 - b_map = length of map in latitude [deg]
3 - desired rms per spectral channel [K]
4 - final width of spectral channel [kHz]
5 - Long or Lat for mapping along Longitude/Latitude [string]
Usage:
ipython calcObsTime.py A_map [deg^2] sigma [K] chanWidth [kHz]
Example:
ipython calcObsTime.py 
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "Nickolas.Pingel@anu.edu.au"
__status__ = "Production"
"""


## imports
import numpy as np
import sys

## assumed constants
SEFD = 10 ## Jy
t_int = 0.5 ## s
l_dump = 1.67 / 60. ## space between data dumps to ensure adequate spatial Nyquist sampling
l_space = 3/ 60. ## row/column spacing in [deg]
A_beam = 1.1331 * (9.1 / 60)**2 # FWHM of beam will always be ~9.1 arcmin at 1.4 GHz for GBT [deg^2]
N_pols = 2 ## number of polarizations
N_formedBeams = 7 ## number of formed beams
N_offs = 8 ## number of off integrations used to construct reference spectrum 
t_off = t_int * N_offs ## time contributing to the reference spectrum 
gain = 1.86 ## gain [K/Jy]


## get user inputs 
## map area
long_map = np.float(sys.argv[1]) # deg
lat_map = np.float(sys.argv[2]) # deg
## desired spectral channel rms
sigma = np.float(sys.argv[3]) # K 

## final width of spectral channel 
width = np.float(sys.argv[4]) * 1000 ## Hz

## Scan along Longitude/Latitude?
scanDir = sys.argv[5]

## area of map
A_map = long_map * lat_map

## first, determine the number of beams in map by dividing provided map area by area of GBT L-Band beam
N_beams = A_map / A_beam

## inform estimated number of beams per map 
print('Number of beams per map: %.2f' % N_beams)

## determine numer of columns/rows (i.e., scans) and integrations
if scanDir == 'Long': ## means we are scanning along rows of constant latitude
	N_scans = np.int(lat_map / l_space) + 1 
	N_ints = np.int(long_map / l_dump) + 1
else:
	N_scans = np.int(long_map / l_space) + 1 
	N_ints = np.int(lat_map / l_dump) + 1

## compute time per scan and thus time per map
t_scan = N_ints * t_int
t_map = t_scan * N_scans

## compute 'On' signal time
t_on = t_map / N_beams

## compute effective integration time per beam 
t_eff = t_on * t_off / (t_on  + t_off)

## report to user
print('Time per map (excluding overhead): %.2f [sec]' % t_map)
print('Total effective integration time in a single map: %.2f [sec/beam]' % t_eff)

## now, utilze the provided sensitivity per channel to calculate the total mapping time
"""
ideal radiometer equation: 
sigma = SEFD / SQRT(N_pols * N_formedBeams * width * t_eff)

Convert user sigma to Jy using assumed gain 

Solve for t_eff knowing the required sigma to determine how many maps are required
"""
sigma_Jy = sigma  / gain

t_effTarget = (SEFD/sigma_Jy)**2 * 1/(N_pols * N_formedBeams * width)

N_maps = t_effTarget / t_eff
totObsTime = N_maps * t_map / 3600. ## hours

## inform user
print('Necessary effective integration time: %.2f [sec]' % t_effTarget)
print('Total number of maps: %.2f' %  N_maps)
print('TOTAL OBSERVING TIME (EXCLUDING OVERHEAD): %.2f [hours]' % totObsTime)






