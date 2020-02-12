PRO edgeoffkeep_indv
;==============================================================================
; Nick Pingel ----- 02/13/20
; Program to calibrate beam-formed spectra by treating the ends of a map as off 
; source positions. This code is specific to PFB (fine channelized) mode in the 
; GBT PAF.
; Calculates Flux* = SEFD*{ON - <OFF>}/<OFF>
; and then corrects for atmosphere --> Flux* = Flux *
; exp{tau0/sin(el)}/eta.
; User input SEFD is in Jy, and derived from the calibration scabns
; ASSUMPTION: tau=0.01
; 
; INPUTS from python wrapper (in order): 
;         fileName => file to read into GBTIDL instance
;         sname => source name (e.g. NGC6946)
;         fileout => (string) name of output file containing calibrated spectra
;         SEFD_X => SEFD for XX polarization 
;         SEFD_Y => SEFD for YY polarization
;         order => order of polynomial to use in baseline fit
;         chanRange => channels over which to perform a baseline subtraction after calibration
;         raDecLimits => list array contining RA and Dec coordinates (in deg) defining region of companion emission to avoid
;=============================================================================== 

;; set to avoid output from divide by zero exceptions
!EXCEPT = 0

;; unpack user arguments passed from python wrapper
args = COMMAND_LINE_ARGS()
fileName  = args[0]
sname = args[1]
outfile = args[2]
SEFD_X = FLOAT(args[3])
SEFD_Y = FLOAT(args[4])
order = FLOAT(args[5])

;; construct channel range array
FOR i=6, n_elements(args) - 5 DO BEGIN
    chanRange = APPEND(chanRange, UINT(args[i]))
ENDFOR

;; construct raDecLimits array (always last four elements, if given)
FOR i=n_elements(args) - 4, n_elements(args) - 1 DO BEGIN
    raDecLimits = APPEND(chanRange, FLOAT(args[i]))
ENDFOR 

;; check if raDecLimits should be undefined (first and second elements equal to 0)
if raDecLimits[0] EQ 0 and raDecLimits[1] EQ 0 THEN DELVAR, raDecLimits

filein, fileName

;; select all scans of a source that are mapping scans
allScans=get_scan_numbers(source=sname,procedure='*Map')   

;;set output filename
fileout, outfile

;; put plotter into chan mode and freeze
chan
freeze

;; loop through scans to process each integration separately; first need to build the average 'off' 
;;for idx = 0, N_ELEMENTS(allScans)-1 DO BEGIN
FOR idx = 0, 1 DO BEGIN
    goodIdx=allScans(idx)

    get, scan=goodIdx, int=0, plnum=0; call scan to get the number of intergrations
    ;; get scan info
    info=scan_info(goodIdx)
    nint=info.n_integrations    
    seqNum = info.PROCSEQN
    evenodd = seqNum MOD 2
    ; segment to check if any integrations fall into restricted range.
    ; If any integrations fall in the resticted range, we select integrations on the opposite side of the map.
    ; If seqNum is even, Ra/Dec will decrease. If seqNum is odd, raDec will increase.
    ; For example, if a companion exists to the north (east), the seqNum in a DecLat (RaLong) map is odd, and the scan has integrations are within the limits,
    ; we use the first 8 integrations to the south (west). The subsequent seqNum will be even, thus the companion emission exists in the first 8 integrations
    ; meaning we require the last 8 integrations to the south (west). The BREAK commands are to BREAK from the FOR loop with the correct OFFS.

    ;set the default offs
    offs=[0,1,2,3,nint-4,nint-3,nint-2,nint-1]
    IF KEYWORD_SET(raDecLimits) THEN BEGIN
        FOR int=0, nint-1 DO BEGIN
            get, scan=goodIdx, int=int, plnum=0
            ;; get RA/Dec coordinates for integratin
            raDecVal=getradec(!g.s[0]) 
            ;; check to see if integration falls within restricted range
            IF (raDecVal[0] GE raDecLimits[0]) AND (raDecVal[0] LT raDecLimits[1]) AND $ 
            (raDecVal[1] GE raDecLimits[2]) AND (raDecVal[1] LT raDecLimits[3]) THEN BEGIN
            
            ;; if yes, decide which offs to set
            ;; if even and integration is in first half of scan, use last 8 intergrations as off
            IF (evenodd EQ 0) AND (int LT nint/2.) THEN BEGIN
                offs=[nint-8,nint-7,nint-6,nint-5,nint-4,nint-3,nint-2,nint-1]
                BREAK
            ;; if even and during second half of scan, use first 8 integrations
            ENDIF ELSE IF (evenodd EQ 0) AND (int GT nint/2.) THEN BEGIN
                offs=[0,1,2,3,4,5,6,7]
                BREAK
            ;; if odd, and during first half of scan, use last 8 integrations
            ENDIF ELSE IF (evenodd EQ 1) AND (int LT nint/2.) THEN BEGIN
                offs=[nint-8,nint-7,nint-6,nint-5,nint-4,nint-3,nint-2,nint-1]
                BREAK
            ;; finally, if odd and in second half, use first 8 integrations
            ENDIF ELSE BEGIN
                offs=[0,1,2,3,4,5,6,7]
                BREAK
            ENDELSE
            ENDIF
        ENDFOR
    ENDIF
    
    ;; now that we have our offs, we need to construct an average spectrum in each polarization to use for calibration
    FOR e = 0, N_ELEMENTS(offs)-1 DO BEGIN ; loop through off integrations and polarizations
        FOR f = 0,1 DO BEGIN
            gettp, goodIdx, intnum=offs[e], plnum=f, /quiet ; get total power for pol
            accum, f ;; accumlate each polarization to later take average
        ENDFOR
    ENDFOR

    ave, 0,/quiet ; <OFF> for pol 0
    copy, 0, 1    ; copy to data buffer 
    ave, 1,/quiet
    copy, 0, 2      

    sclear, 0 ; clear accumulators 1 and 2
    sclear, 1 

    ;; loop through ons and polarization to calibrate 
    FOR g = 0, nint - 1 DO BEGIN
        FOR f = 0, 1 DO BEGIN
            gettp, goodIdx, int = g, plnum = f, /quiet ;; get on integration
            copy, 0, f + 3
        ENDFOR
    
        ;; apply calibration for both polarizations 
        ;; subtract on-off
        subtract, 3, 1, 5 ; on - <off>, put in data buffer 5 (pol YY/XX)
        subtract, 4, 2, 6 
        
        ;; process YY Pol   
        ;; divide off, put in data buffer 0 for scaling
        divide, 5, 1, 0

        ;; scale to units of Flux*
        scale, SEFD_Y*exp(0.01/sin(!g.s[0].elevation*3.14159/180.0))/0.99 ; Ta*, will default to buffer 0 
        ;; set region for baseline subtraction
        nregion,chanRange
        IF KEYWORD_SET(order) THEN baseline,nfit=order
        x=dcextract(!g.s[0],floor(chanRange[0]),ceil(chanRange[n_elements(chanRange)-1])) ; extracts spectrum over regions of interest
        
        ;; set in data container to keep
        set_data_container,x
        
        ;; free for later use
        data_free,x
        
        ;; set units
        !g.s[0].units='Jy'
        keep 

        ;; process XX Pol
        ; divide off, put in data buffer 0 for scaling
        divide, 6, 2, 0
        
        ; scale to units of Ta*
        scale, SEFD_X*exp(0.01/sin(!g.s[0].elevation*3.14159/180.0))/0.99 ; Ta*, will default to buffer 0
        
        ;; set region for baseline subtraction
        nregion, chanRange
        IF KEYWORD_SET(order) THEN baseline,nfit=order
        x=dcextract(!g.s[0],floor(chanRange[0]),ceil(chanRange[n_elements(chanRange)-1])) ; extracts spectrum over regions of interest
        
        ;; set in data container to keep
        set_data_container,x
        
        ;; free for later use
        data_free,x
        
        ;; set units
        !g.s[0].units='Jy'
        keep
    ENDFOR
ENDFOR

;; unfreeze plotter
unfreeze

exit
END
    





