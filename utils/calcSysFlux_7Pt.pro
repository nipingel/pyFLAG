;; define function to compute flux
FUNCTION calcFlux, coeffArr, vG
        a0 = coeffArr[0]
        a1 = coeffArr[1]
        a2 = coeffArr[2]
        a3 = coeffArr[3]
        a4 = coeffArr[4]
        a5 = coeffArr[5]
        logS = a0 + a1*ALOG10(vG) + a2*ALOG10(vG)^2 + a3*ALOG10(vG)^3 + a4*ALOG10(vG)^4 + a5*ALOG10(vG)^5
        S = 10.^logS
        RETURN, S
END
;; function to compute error on flux
FUNCTION calcFluxErr, errCoeffArr, flux, vG
        a0Err = errCoeffArr[0]
        a1Err = errCoeffArr[1]
        a2Err = errCoeffArr[2]
        a3Err = errCoeffArr[3]
        a4Err = errCoeffArr[4]
        a5Err = errCoeffArr[5]
        sqFluxErr = a0Err^2 + a1Err^2*ALOG10(vG)^2 + a2Err^2*ALOG10(vG)^4 + a3Err^2*ALOG10(vG)^5 + a4Err^2*ALOG10(vG)^6 + a5Err^2*ALOG10(vG)^7
        fluxErr = SQRT(sqFluxErr)
        fluxErr = flux*ALOG(10)*fluxErr
        RETURN, fluxErr
        END

PRO calcSysFlux_7Pt, sname, chan0, chan1, chan2, chan3, scan0, pfb
;==============================================================================
; Nick Pingel ----- 08/04/19
; Program to calculate the system equivalent flux density (SEFD) given a known calibrator.
; Ssrc is computed from Equation 1 in Perley & Butler (2017).
; The SEFD is computed by determining the maximum power value in all pointed in the hexagonal beam 
; pattern scans. This is effectively P_on. The uncertainty in this value is the standard deviation given by 
; a Gaussian fit to the distribution of all power values in the bandpass between user provided channels chan0 - chan3
; for the grid scans ONLY. P_off is the fitted mean (and associated standard deviation) from a Gaussian fit
; to the distribution of power values in the bandpass between user provided channels chan0 - chan3 for the 
; two OFF scans ONLY. User needs to provide the eight scans to process.
; The final SEFD is given by the equation: SEFD = Ssrc*P_off/(P_on - P_off). 
; The program also propogates all uncertainties and gives the final uncertainty in the SEFD.   
; Finally, if pfb is set to 1, outlying values caused by 3 dB drop in power are dropped form distribution
; INPUTS: sname => (string) source name (e.g. 3C295)
;   chan0 => first channel denoting the region over which to draw power values
;   chan1 => second channel denoting the region over which to draw power values
;         chan2 => third channel denoting the region over which to draw power values
;         chan3 => fourth channel denoting the region over which to draw power values
;   scan0 => scan of the first off, which is used to create an array of the correct scan numbers
;   pfb => 1 if in pfb mode, 0 if not
; EXAMPLE: calcSysFlux_7pt, '3C295', 125, 135, 155, 180, 5, 1
;===============================================================================
compile_opt idl2
;; set pfb to 0 if not provided
IF NOT KEYWORD_SET(pfb) THEN pfb = 0

;; set constants
kB = 1.3806e-23

;; global variable defined for calibrator flux equation
freqVal = 1.420405752 ;; GHz

scanArr = INDGEN(9)
;; loop to fill scan array with correct scan numbers
FOR i =0, 8 DO BEGIN
  scanArr[i] = scan0 + i 
ENDFOR

print, scanArr

freeze
chan

;; define dictionaries to hold flux
fluxCoeffStruct = create_struct('a', [1.3253,-0.7553,-0.1914,0.0498, 0.0, 0.0], 'b', [1.8017,-0.7884,-.1035,-0.0248,0.009,0.0], 'c', [1.0088,-0.4981,-0.155,-0.010,0.022,0.0], 'd', [1.4516,-0.6961,-.201,0.064,-0.046,0.029], 'e', [2.4466, -0.8116, -0.048, 0.0, 0.0, 0.0], 'f', [1.2481,-0.4507,-0.1798,0.0357, 0.0, 0.0], 'g', [1.4701,-0.7658,-0.2780,-0.0347,0.0399, 0.0, 0.0], 'h', [1.830, -1.05, -0.03, -0.05, -0.03, 0.02], 'i', [1.8627, -0.6938, -0.100, -0.032, 0.0, 0.0], 'j', [3.3498, -1.0022, -0.225, 0.023, 0.043, 0.0])

;; define dictionary to hold flux uncertainties
fluxErrCoeffStruct = create_struct('a', [0.0005, 0.0009, 0.0011, 0.0009, 0, 0], 'b', [0.0007, 0.0012, 0.0023, 0.0013, 0.0013, 0], 'c', [0.0009, 0.0022, 0.003, 0.07, 0.003, 0], 'd', [0.0010, 0.0017, 0.005, 0.004, 0.004, 0.003], 'e', [0.0007, 0.0020, 0.003, 0, 0, 0], 'f', [0.0005, 0.0009, 0.0011, 0.0009, 0, 0], 'g', [0.007, 0.0012, 0.0023, 0.0013, 0.0013, 0], 'h', [0.005, 0.03, 0.06, 0.09, 0.03, 0.04], 'i', [0.0010, 0.0014, 0.005, 0.0005, 0, 0], 'j', [0.0010, 0.0014, 0.006, 0.002, 0.005, 0])

;; define dictionary to hold source coordinates
sourceCoordStruct = create_struct('a', [24.422081, 33.159759], 'b', [69.268230, 29.670505], 'c', [80.291192, 16.639459], 'd', [85.650575, 49.852009],'e', [187.705930, 12.391123], 'f', [202.784533, 30.509155], 'g', [212.835495, 52.202770], 'h', [252.783945, 4.992588], 'i', [260.117326, -0.979617], 'j', [299.868153, 40.733916])

;; 3C48
IF sname EQ '3C48' THEN BEGIN
  fluxCoeffArr = fluxCoeffStruct.a
  fluxEffCoeffArr = fluxErrCoeffStruct.a
    sourceCoords = sourceCoordStruct.a
    sourceMaj = sourceCoords[0]
    sourceMin = sourceCoords[1]
ENDIF
;; 3C123
IF sname EQ '3C123' THEN BEGIN
        fluxCoeffArr = fluxCoeffStruct.b
        fluxEffCoeffArr = fluxErrCoeffStruct.b
        sourceCoords = sourceCoordStruct.b
        sourceMaj = sourceCoords[0]
        sourceMin = sourceCoords[1]
ENDIF
;;3C138
IF sname EQ '3138' THEN BEGIN
        fluxCoeffArr = fluxCoeffStruct.c
        fluxEffCoeffArr = fluxErrCoeffStruct.c
        sourceCoords = sourceCoordStruct.c
        sourceMaj = sourceCoords[0]
        sourceMin = sourceCoords[1]
ENDIF
;; 3C147
IF sname EQ '3C147' THEN BEGIN
        fluxCoeffArr = fluxCoeffStruct.d
        fluxEffCoeffArr = fluxErrCoeffStruct.d
        sourceCoords = sourceCoordStruct.d
        sourceMaj = sourceCoords[0]
        sourceMin = sourceCoords[1]
ENDIF
;; VirgoA
IF sname EQ 'VirgoA' THEN BEGIN
        fluxCoeffArr = fluxCoeffStruct.e
        fluxEffCoeffArr = fluxErrCoeffStruct.e
        sourceCoords = sourceCoordStruct.e
        sourceMaj = sourceCoords[0]
        sourceMin = sourceCoords[1]
ENDIF
;; 3C286
IF sname EQ '3C286' THEN BEGIN
        fluxCoeffArr = fluxCoeffStruct.f
        fluxEffCoeffArr = fluxErrCoeffStruct.f
        sourceCoords = sourceCoordStruct.f
        sourceMaj = sourceCoords[0]
        sourceMin = sourceCoords[1]
ENDIF
;; 3C295
IF sname EQ '3C295' THEN BEGIN
        fluxCoeffArr = fluxCoeffStruct.g
        fluxEffCoeffArr = fluxErrCoeffStruct.g
        sourceCoords = sourceCoordStruct.g
        sourceMaj = sourceCoords[0]
        sourceMin = sourceCoords[1]
ENDIF
;; 3C348
IF sname EQ '3C348' THEN BEGIN
        fluxCoeffArr = fluxCoeffStruct.h
        fluxEffCoeffArr = fluxErrCoeffStruct.h
        sourceCoords = sourceCoordStruct.h
        sourceMaj = sourceCoords[0]
        sourceMin = sourceCoords[1]
ENDIF
;; 3C353
IF sname EQ '3C353' THEN BEGIN
        fluxCoeffArr = fluxCoeffStruct.i
        fluxEffCoeffArr = fluxErrCoeffStruct.i
        sourceCoords = sourceCoordStruct.i
        sourceMaj = sourceCoords[0]
        sourceMin = sourceCoords[1]
ENDIF
;; CygnusA
IF sname EQ 'CygnusA' THEN BEGIN
        fluxCoeffArr = fluxCoeffStruct.j
        fluxEffCoeffArr = fluxErrCoeffStruct.j
        sourceCoords = sourceCoordStruct.j
        sourceMaj = sourceCoords[0]
        sourceMin = sourceCoords[1]
ENDIF


calFlux = calcFlux(fluxCoeffArr, freqVal)
calFluxErr = calcFluxErr(fluxEffCoeffArr, calFlux, freqVal)

PRINT, 'Source flux is: ' + STRING(calFlux) + '+/-' + STRING(calFluxErr)

;; Loop through all scans/integrations. Take signal power value to be the maximum power value at channel 150 -> 1420 MHz. 
;; The scan corresponding to the maximum power value will be used will be used to build up a statistical sample. 
;; The rms will come from a Gaussian fit to this distribution of bandpass power values. 

;; initial power value to be used for comparison
;; loop through signal scans
power0 = 0
FOR i=1, 7 DO BEGIN
    scanNum = scanArr[i]
    print, scanNum
    info=scan_info(scanNum)
    nint = info.n_integrations
    ;; loop through integrations
    FOR j=0, nint -1 DO BEGIN
        gettp, scanNum, int = j, plnum = 1
        ;; get power value at channel 150
        data = getdata()
        powerVal = data[150]
        ;; check if greater than previous power value. If yes, record scan, integration, and set comparison value
        IF powerVal GT power0 THEN BEGIN
            sigScan = scanNum
            maxInt = j
            power0 = powerVal
        ENDIF
    ENDFOR
ENDFOR

print, 'Signal Scan: ', sigScan


;; Now that the signal scan is known, loop through the integrations to build up the statistical distribution
;; Take the signal power value to be the average value over all integrations of the signal scan at channel 150. 

;; clear accum buffers
sclear, 0
sclear, 1

;; select signal scan to get number of integrations
get, scan=sigScan, int=0, plnum=1

;; select the correct scan based on user provided beam number
beamAz = !g.s[0].azimuth
beamEl = !g.s[0].elevation

;; conver to cross-el
beamXEl = beamAz*cos(beamEl*!DTOR)

PRINT, 'Working on scan: ', sigScan
;; get integrations
info=scan_info(sigScan)
nint = info.n_integrations
FOR j=0, nint -1 DO BEGIN
  FOR pl = 0, 1 DO BEGIN
    gettp, sigScan, int=j, plnum=pl, /quiet
    ;; get elevation value for correction
    elVal = !g.s[0].elevation
    ;; compute atmospheric correction factor
    atmosCorr = exp(0.01/sin(elVal*3.14159/180.0))/0.99 ;; numerator is for atmoshperic attenuation; denominator is correction for spillover, ohmic loss, etc...

    bandpass = getdata()
    ;; combine the two bandpass into final on distribution (for each polarization)
    IF pl EQ 0 THEN BEGIN
      onDist_YY = APPEND(onDist_YY, bandpass[chan0:chan1]*atmosCorr)
      onDist_YY = APPEND(onDist_YY, bandpass[chan2:chan3]*atmosCorr)
    ENDIF
    IF pl EQ 1 THEN BEGIN
      onDist_XX = APPEND(onDist_XX, bandpass[chan0:chan1]*atmosCorr)
      onDist_XX = APPEND(onDist_XX, bandpass[chan2:chan3]*atmosCorr)
    ENDIF
  ENDFOR
ENDFOR


;; get average value as 'OnPower' value
onPower_YY = MEAN(onDist_YY)
onPower_XX = MEAN(onDist_XX)

print, 'Ave Total Power (YY): ', onPower_YY
print, 'Ave Total Power (XX): ', onPower_XX

;; fit a gaussian to the 'on' distributions to determine the uncertainty
HISTOGAUSS, onDist_YY, onGaussCoeff_YY
HISTOGAUSS, onDist_XX, onGaussCoeff_XX

;; define the on power values and associated uncertainties
onPowerErr_YY = onGaussCoeff_YY[2]
onPowerErr_XX = onGaussCoeff_XX[2]
;; if histogauss fails to converge, then take straight standard devation of the distribution as
;; the uncertainty
IF ABS(onPowerErr_YY) GT 10 THEN onPowerErr_YY = STDDEV(onDist_YY)
IF ABS(onPowerErr_XX) GT 10 THEN onPowerErr_XX = STDDEV(onDist_XX)

;; decide which off scan is closer to use as the reference scan
;; get coordinates of 1st off
get, scan=scanArr[0], int=0, plnum=1
off1Az = !g.s[0].azimuth
off1El = !g.s[0].elevation

;; conver to cross-el
off1Xel = off1Az*cos(off1El*!DTOR)

;; get coordinates of 2nd off
get, scan=scanArr[8], int=0, plnum=1
off2Az = !g.s[0].azimuth
off2El = !g.s[0].elevation

off2Xel = off2Az*cos(off2El*!DTOR)

;; determine the angular distance between these coordinates
d1 = SQRT((beamXel - off1Xel)^2 + (beamEl - off1El)^2)
d2 = SQRT((beamXel - off2Xel)^2 + (beamEl - off2El)^2)

;; see which is closer, select that as a reference, and construct a reference spectrum
IF d1 LT d2 THEN refScan = scanArr[0]
IF d1 GT d2 THEN refScan = scanArr[8]
info=scan_info(scanArr[0])
nint=info.n_integrations
sclear, 0
sclear, 1

;; build up statistical distribution of power samples from selected reference scan.
;; Take off power value to be the average value over all integrations of reference scan 
;; at channel 150. 

FOR i=0, nint-1 DO BEGIN
  FOR pl=0,1 DO BEGIN
    gettp, refScan, int=i, plnum=pl

    ;; get elevation value for correction
    elVal = !g.s[0].elevation
    ;; compute atmospheric correction factor
    atmosCorr = exp(0.01/sin(elVal*3.14159/180.0))/0.99 ;; numerator is for atmoshperic attenuation; denominator is correction for spillover, ohmic loss, etc...
    bandpass = getdata()
      IF pl EQ 0 THEN BEGIN
        offDist_YY = APPEND(offDist_YY, bandpass[chan0:chan1]*atmosCorr)
        offDist_YY = APPEND(offDist_YY, bandpass[chan2:chan3]*atmosCorr)
      ENDIF
      IF pl EQ 1 THEN BEGIN
        offDist_XX = APPEND(offDist_XX, bandpass[chan0:chan1]*atmosCorr)
        offDist_XX = APPEND(offDist_XX, bandpass[chan2:chan3]*atmosCorr)
      ENDIF
  ENDFOR
ENDFOR

;; take average
offPower_YY = MEAN(offDist_YY)
offPower_XX = MEAN(offDist_XX)

;; fit gaussian to the 'off' distributions to determine the mean and associated uncertainty
HISTOGAUSS, offDist_YY, offGaussCoeff_YY
HISTOGAUSS, offDist_XX, offGaussCoeff_XX

offPowerErr_YY = offGaussCoeff_YY[2]
offPowerErr_XX = offGaussCoeff_XX[2]

;; if histogauss failed to converge, then take straight standard devation of the distribution as
;; the uncertainty
IF ABS(offPowerErr_YY) GT 10 THEN offPowerErr_YY = STDDEV(offDist_YY)
IF ABS(offPowerErr_XX) GT 10 THEN offPowerErr_XX = STDDEV(offDist_XX)

print, 'Off Power (YY): ', offPower_YY
print, 'Off Power (XX): ', offPower_XX

;; compute YY sensitivity at channel 150
Sens_YY = 2*kB/(1e-26*calFlux)*onPower_YY/offPower_YY

;;compute YY sensitivity uncertainty
sig_sens_YY = 2*kB*onPower_YY/(1e-26*calFlux*offPower_YY)*SQRT(onPowerErr_YY^2/onPower_YY^2 + offPowerErr_YY/offPower_YY^2 + calFluxErr^2/calFlux^2)

print, 'YY Sensitivity [m^2 K^-1], ', Sens_YY
print, '+/-', sig_sens_YY

;; compute XX sensitivity at channel 150
Sens_XX = 2*kB/(1e-26*calFlux)*onPower_XX/offPower_XX

;;compute XX sensitivity uncertainty
sig_sens_XX = 2*kB*onPower_XX/(1e-26*calFlux*offPower_XX)*SQRT(onPowerErr_XX^2/onPower_XX^2 + offPowerErr_XX/offPower_XX^2 + calFluxErr^2/calFlux^2)

print, 'XX Sensitivity [m^2 K^-1], ', Sens_XX
print, '+/-', sig_sens_XX

;; calculate YY system flux
Ssys_YY = calFlux*offPower_YY/(onPower_YY - offPower_YY)
;; calculate YY system flux uncertainties
;; first term is from the calibrator flux
firstTerm_YY = ((offPower_YY)/(onPower_YY - offPower_YY))^2*calFluxErr^2
secondTerm_YY = ((calFlux*offPower_YY)/(onPower_YY - offPower_YY)^2)^2*onPowerErr_YY^2
thirdTerm_YY = ((calFlux*onPower_YY)/(onPower_YY - offPower_YY)^2)^2*offPowerErr_YY^2
SsysErr_YY = SQRT(firstTerm_YY + secondTerm_YY + thirdTerm_YY)

;; calculate XX system flux
Ssys_XX = calFlux*offPower_XX/(onPower_XX - offPower_XX)

;; calculate YY system flux uncertainties
;; first term is from the calibrator flux
firstTerm_XX = ((offPower_XX)/(onPower_XX - offPower_XX))^2*calFluxErr^2
secondTerm_XX = ((calFlux*offPower_XX)/(onPower_XX - offPower_XX)^2)^2*onPowerErr_XX^2
thirdTerm_XX = ((calFlux*onPower_XX)/(onPower_XX - offPower_XX)^2)^2*offPowerErr_XX^2
SsysErr_XX = SQRT(firstTerm_XX + secondTerm_XX + thirdTerm_XX)
print, 'YY System Flux Density: ', Ssys_YY
print, '+/-', SsysErr_YY

print, 'XX System Flux Density: ', Ssys_XX
print, '+/-', SsysErr_XX

END


