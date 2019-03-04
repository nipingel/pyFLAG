PRO PAF_chanShift, sname, chanShiftVal = chanShiftVal, fileout=outFile
;==============================================================================
; Nick Pingel ----- 08/14/18
; Program to shift channels channels of the lower LO spectra before combining with a higher Dithered LO spectra. The inputs are the 
; source name (to select the correct scans), value of channels to shift, and the name of the output file. Since the LO is dithered by
; 0.30318/2 = 0.15159 MHz, or 16 fine channels (9.57 kHz wide), and we've smoothed to 24.414 kHz, the number of channels to shift by will
; usually be 16/2.5 = -6.2. It's negative since the lower LO setting is 151.18 kHz LESS than the higher LO setting
; sructure. 
; INPUTS: sname => (string) source name (e.g. NGC6946)
;         chanShiftVal => (double) (e.g. -6.2)
;         fileout => (string) name of output file containing calibrated spectra
; EXAMPLE: PAF_chanShift, 'NGC6946', fileout='AGBT16B_400_12_NGC6946_chanBlank.fits' 
;=============================================================================== 

;; select all scans of a source that are mapping scans
allScans=get_scan_numbers(source=sname,procedure='*Map')
fileout, outfile

freeze
chan
for idx = 0,N_ELEMENTS(allScans)-1 DO BEGIN
    goodIdx=allScans(idx)

    print,'Working on scan: ', goodIdx
    get, scan=goodIdx, int=0, plnum=1; call scan to get the number of intergrations
    ;; get scan info
    info=scan_info(goodIdx)
    nint=info.n_integrations
    for int=0,nint-1 do begin 
        for pl=0,1 do begin ;; only care about XX polarization
           get, scan=goodIdx, int=int, plnum=pl
           dcshift, !g.s[0], chanShiftVal
           keep
        ENDFOR
    ENDFOR
ENDFOR
unfreeze
END
