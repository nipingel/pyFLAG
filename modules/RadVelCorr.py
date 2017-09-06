"""
09/5/17
Based heavily on Frank Ghingo's online calcultor (all of the logic and computations are his code exactly), this script will take command line input (eventually; now 
just opens a list of Antenna FITS files to get information)
and calculate the velocity corrections needed to go from topocentric to barycenter/LSR radial velocities. This is useful when this calculator needs to be called
many times. 

Usage:
ipython radvelcalc.py UT_Date UT_Time RA Dec
Example:
ipython calcSensitivity.py  
INPUTS:
UT_Date: UT date in the format of YYYY/MM/DD
UT_Time: UT Time in the format of HOUR:MIN:SEC
RA: Right Ascension with the format of HOUR:MIN:SEC
Dec: Declination with the format of HOUR:MIN:SEC 

__author__ = "Nick Pingel"
__version__ = "0.1"
__email__ = "nipingel@mix.wvu.edu"
__status__ = "beta"
"""
# radvel.py - to be calculator for radial velocities
# note: needs /opt/local/bin/python to get slalib
#
import cgi
import sys
import math
import slalib as s
import numpy as np
from astropy.time import Time



class RadVelCorr:
  import cgi
  import sys
  import math
  import slalib as s
  import numpy as np
  from astropy.time import Time
  from astropy.io import fits
  import matplotlib.pyplot as pyplot
  import matplotlib 
  matplotlib.rc('font', family='serif')
  matplotlib.rc('font', serif='Avant Garde')
  matplotlib.rc('text', usetex='False')
  matplotlib.rcParams.update({'font.size': 14})

  def __init__(self):
      import cgi
      import sys
      import math
      import slalib as s
      import numpy as np
      from astropy.time import Time
      from astropy.io import fits
      import matplotlib.pyplot as pyplot
      import matplotlib
      matplotlib.rc('font', family='serif')
      matplotlib.rc('font', serif='Avant Garde')
      matplotlib.rc('text', usetex='False')
      matplotlib.rcParams.update({'font.size': 14})
      return
  
  # longitude and latitude of the GBT
  GBTLAT  = 38.4331294 * math.pi/180.0 # in radians
  GBTLONG = 79.8398397 * math.pi/180.0 # in radians
  GBTHGT  = 824.36                     # meters above the ellipsoid

  # longitude and latitude of the 20-meter
  LAT20M  = 38.436850 * math.pi/180.0
  LONG20M = 79.825518 * math.pi/180.0
  HGT20M  = 835.0

  # 140-foot telescope:
  LAT140  = 38.437823 * math.pi/180.0
  LONG140 = 79.835778 * math.pi/180.0
  HGT140  = 812.6

  global tellat, tellong, telelev
  tellat  = GBTLAT
  tellong = GBTLONG
  telelev = GBTHGT

  #-----------------------------------------------------------------
  # hlst returns the local apparent sidereal time (LAST) in radians
  # input is modified Julian date (MJD)
  #-----------------------------------------------------------------
  def hlst(self,mjd) :
    global tellat, tellong, telelev

    # test 0.2 sec UT1 correction
    # dut1 = (0.2 /3600.0) * math.pi/12.0
    dut1 = 0.0

    last = s.sla_gmst(mjd) - tellong + s.sla_eqeqx(mjd) + dut1
   
    # lmst = s.sla_gmst(mjd) - tellong
    if last < 0.0 : last = last + 2.0*math.pi
    return last 
  #--------------------------------------------------------------


  #---------------------------------------------------------
  # convert sexagesimal string  hh:mm:ss. or +dd:mm:ss.s to hours or degrees.
  # there can be a leading sign (+ or -) on degrees.
  # separators can be colons, commas, or blanks
  # if not separators, just return the number.
  # returns hours or degrees , whatever units are input
  #---------------------------------------------------------
  def sex2deg(self, sstring ) :
    ss1 = sstring.strip()     # remove leading and trailing blanks if any
    ss2 = ss1.split(':')      # can use colon, comma, or blanks
    if len(ss2) <= 1 :
      ss2 = ss1.split(',')
      if len(ss2) <= 1 :
        ss2 = ss1.split(' ')

    # if no delimiters, just return the number
    if len(ss2) <= 1 : return float(ss2[0]), 1

    # print 'ss2=',ss2
    isign = 1.0
    dd1 = int(ss2[0])
    mm1 = int(ss2[1])
    secs = 0.0
    # check for leading sign
    if (ss2[0]).find('-') == 0 : 
      isign = -1.0
      dd1 = -dd1

    # seconds may or may not be there
    if len(ss2) > 2 :
      secs = float(ss2[2])

    angdeg = isign * (dd1 + (mm1 + secs/60.0)/60.0 )

    return angdeg, 0

  #---------------------------------------------------------

  #--------------------------------------------------------------
  # mjdd is modified Julian date and fraction
  # radeg is right ascension in degrees
  # decdeg is declination in degrees
  #--------------------------------------------------------------
  def rvel(self, mjdd, radeg, decdeg) :
    global tellat, tellong, telelev

    last = self.hlst(mjdd)   # apparent LST in radians

    # convert ra,dec to radians
    rarad = radeg * math.pi/180.0
    dcrad = decdeg * math.pi/180.0

    # convert star position to vector
    starvect = s.sla_dcs2c(rarad, dcrad)

    # velocity component in ra,dec due to Earth rotation
    Rgeo = s.sla_rverot( tellat, rarad, dcrad, last)

    # get Barycentric and heliocentric velocity and position of the Earth.
    evp = s.sla_evp(mjdd, 2000.0)
    dvb = evp[0]   # barycentric velocity vector, in AU/sec
    dpb = evp[1]   # barycentric position vector, in AU
    dvh = evp[2]   # heliocentric velocity vector, in AU/sec
    dph = evp[3]   # heliocentric position vector, in AU

    # dot product of vector to object and heliocentric velocity
    # convert AU/sec to km/sec
    vcorhelio = -s.sla_dvdv( starvect, dvh) *149.597870e6
    vcorbary  = -s.sla_dvdv( starvect, dvb) *149.597870e6

    # rvlsrd is velocity component in ra,dec direction due to the Sun's motion with
    # respect to the "dynamical" local standard of rest
    rvlsrd = s.sla_rvlsrd( rarad,dcrad)

    # rvlsrk is velocity component in ra,dec direction due to the Sun's motion with
    # respect to the "kinematic" local standard of rest
    rvlsrk = s.sla_rvlsrk( rarad,dcrad)

    # rvgalc is velocity component in ra,dec direction due to the rotation of the Galaxy.
    rvgalc = s.sla_rvgalc( rarad,dcrad)

    totalhelio = Rgeo + vcorhelio
    totalbary  = Rgeo + vcorbary

    totallsrk = totalhelio + rvlsrk
    totalgal  = totalbary  + rvlsrd + rvgalc

    # ('UTHr', 'LST', 'Geo', 'Helio', 'Total Helio', 'LSRK', 'Galacto')
    #(hour,last*12.0/math.pi,Rgeo,vcorb,totalhelio,totallsrk,totalgal)

    # resulting velocities in km/sec
    return ((mjdd, last*12.0/math.pi, Rgeo, totalhelio,totalbary, totallsrk,totalgal))

  #--------------------------------------------------------------



  #---------------------------------------------------
  # radvellist
  # give it a date as a string 'yyyy/mm/dd'
  # it makes output for 25 hours starting that date
  # rastring and decstring are the RA and DEC of the source
  #  in the form 'hh mm ss.ss', 'dd mm ss.s' (need blank delimiters)
  #---------------------------------------------------
  def radvellist(self, sdate, rastring, decstring) :

    # decode the date, convert to MJD
    ss = sdate.split('-')
    year = int(ss[0])
    mont = int(ss[1])
    day  = int(ss[2])
    mjdv = s.sla_cldj(year,mont,day)
    mjd0 = mjdv[0]
    rstat = mjdv[1]

    raval = s.sla_dafin( rastring, 1)  # convert to radians
    dcval = s.sla_dafin( decstring, 1)

    rarad = raval[1] * 15.
    dcrad = dcval[1]

    # degrees
    radeg  = rarad*180.0/math.pi
    decdeg = dcrad*180.0/math.pi

    print ""
    print ' %3s %7s %7s %7s %13s %13s %13s' %  \
          ('UTHr', 'LST', 'Geo', 'Helio', 'Barycentric', 'LSRK', 'Galacto')

    # loop through 25 hours
    for i in range(0,25) :
      hour = i
      dfrac = hour/24.0
      mjd1 = mjd0 + dfrac

      rr = self.rvel(mjd1, radeg,decdeg)
      #  ((mjdd, last*12.0/math.pi, Rgeo, totalhelio,totalbary, totallsrk,totalgal))

      print ' %3d %7.3f %7.3f %7.3f %13.3f %13.3f %13.3f' % \
         (hour,rr[1],rr[2],rr[3],rr[4],rr[5],rr[6])
     
  #--------------------------------------------------------------
  # date: 'yyyy/mm/dd', time: 'hh:mm:ss' (in UT time)  
  # rastr: 'hh:mm:ss', decstr: 'dd mm dd'  (with colons or blanks)
  #def radvel( sdate, time, rastring, decstring) :
  def radvel(self, sdate, time, radeg, decdeg) :
    # decode the date, convert to MJD
    ss = sdate.split('-')
    year = int(ss[0])
    mont = int(ss[1])
    day  = int(ss[2])
    mjdv = s.sla_cldj(year,mont,day)
    mjd0 = mjdv[0]
    rstat = mjdv[1]

    uttime,tf1 = self.sex2deg(time)           # time in hours

    #radeg,tf1 = sex2deg(rastring) 
    #if tf1 == 0 : radeg = radeg*15.0

    #decdeg,tf1 = sex2deg(decstring)
    # print '<br>RA=%8.4f deg, DEC=%8.4f deg' % (radeg,decdeg)

    mjdt = mjd0 + uttime/24.0

    # print 'Using ', sdate, mjdt, uttime, radeg,decdeg

    rr = self.rvel( mjdt, radeg, decdeg)
    #  ((mjdd, last*12.0/math.pi, Rgeo, totalhelio,totalbary, totallsrk,totalgal))

    return rr
  #--------------------------------------------------------------




  #--------------------------------------------------------------

  def correctVel(self,utdate, uttime, rastr, decstr) :
    #global utdate,uttime,rastr,decstr,vbary,vlsrk
    f0 = 1420405751.7667 ## HI line; [Hz]

    rr = self.radvel( utdate, uttime, rastr, decstr)
    vhelio = rr[4]
    vlsrk  = rr[5]
    vgalac = rr[6]

    #print('vHelio:' + np.str(vhelio))
    #print('vlsrk: ' + np.str(vlsrk))
    #print('vgalac ' + np.str(vgalac))

    ## do frequency conversion
    vlsrk = vlsrk*1000 ## convert to m/s
    vhelio = vhelio*1000 
    
    return vhelio  



