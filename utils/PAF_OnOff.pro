PRO PAF_OnOff, sname, chanRange = chanRange, Tsys_Y = Tsys_Y, Tsys_X = Tsys_X, order = order, fileout = outfile
;==============================================================================
; Nick Pingel ----- 08/02/17
; Program to calibrate beam-formed spectra from On/off pointings.
; This code is specific to PFB (fine channelized) mode in the GBT PAF (FLAG).
; Calculates Ta = Tsys*{ON - OFF}/OFF
; and then corrects for atmosphere --> Ta* = Ta *
; exp{tau0/sin(el)}/eta.
; User input Tsys is in Kelvin, and derived from the calibration grid
; ASSUMPTION: tau=0.01
;
; INPUTS: sname => (string) source name (e.g. NGC6946)
;         Tsys_Y/_X => (float) system temperature derived from Tsys/eta output from calibratio grid; eta = 0.65
;         chanRange => (float array) channels over which to perform a baseline subtraction after calibration
;         order => (float) order of polynomial to use in baseline fit
;         badScans = > (float array) array of bad scan numbers to remove from analysis
;         fileout => (string) name of output file containing calibrated spectra
; EXAMPLE: PAF_OnOff, 'NGC6946', Tsys_Y=22.56, Tsys_X=20, chanRange=[10, 1000, 2000, 3150], order=3, fileout='AGBT16B_400_09_NGC6946_OnOff_Beam0.fits'
;===============================================================================


allScans = get_scan_numbers(source=sname, procedure='OnOff')
;; remove badscans
IF NOT KEYWORD_SET(badScans) THEN badScans=[0]
goodScans=setdifference(allScans,badScans)


fileout, outfile
chan
freeze
buff = 0
;FOR i=0, n_elements(allScans)-1 DO BEGIN 
for i=4, 5 DO BEGIN
    scanNum = allScans[i] ;; get scans
    get, scan=scanNum, int=0, plnum=0 ;; get first scan/integration for info
    scanNum = allScans[i] ;; get proper scan number
    info = scan_info(scanNum) ;; get info
    nints = info.n_integrations ;; number of integrations
    seqNum = info.procseqn ;; 1 for On, 2 for Off in OnOff scans 
    FOR j=0,nints-1 DO BEGIN  ;; loop through integrations in a scan 
        PRINT, 'Scan: '+ string(scanNum) +'; ' + 'int: ' + string(j)
        FOR p=0,1 DO BEGIN ;; loop through polarizations
            PRINT, 'Scan: '+ string(scanNum) +'; ' + 'int: ' + string(j) + '; ' + 'pol: ' + string(p)
            get, scan=scanNum, int=j, plnum=p ;; get pol
            ;; if ON scan, and YY pol, place in accum buffer 0
            IF (seqNum MOD 2 EQ 1) AND (p EQ 0) THEN BEGIN
                accum, p
            ;; if ON scan and XX pol, place in accum buffer 1
            ENDIF ELSE IF (seqNum MOD 2 EQ 1) AND (p EQ 1) THEN BEGIN
                accum, p
            ;; if OFF scan and YY pol, place in accum buffer 2
            ENDIF ELSE IF (seqNum MOD 2 EQ 0) AND (p EQ 0) THEN BEGIN
                accum, p+2
            ;; if OFF scan and XX pol, place in accum buffer 3
            ENDIF ELSE BEGIN accum, p+2
            ENDELSE
        ENDFOR
    ENDFOR
ENDFOR
    

;; process pols of ON scan
ave, 0 ;; average YY/ON
copy, 0, 1
ave, 1 ;; average XX/ON
copy, 0, 2
sclear, 0 ;; clear accum buffers
sclear, 1
    
;; process pols of OFF scan    
ave, 2 ;; average YY/OFF
copy, 0, 3
ave, 3 ;; average XX/OFF
copy, 0, 4
sclear, 2 ;; clear accum buffers
sclear, 3


;; on-off/off Ypol 
subtract, 1, 3, 0
unfreeze
divide, 0, 3
scale, Tsys_Y
!g.s[0].units='Ta [K]'
nregion = chanRange
;;nfit, order
;;baseline
keep

;;on-off/off xpol
subtract, 2, 4, 0
divide, 0, 4
scale, Tsys_X
!g.s[0].units='Ta [K]'
nregion = chanRange
;;nfit, order
;;baseline
keep

unfreeze

END

 
