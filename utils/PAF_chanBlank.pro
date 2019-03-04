PRO PAF_chanBlank, sname, fileout=outFile
;==============================================================================
; Nick Pingel ----- 08/01/18
; Program to blank (i.e., count values are set to NaN) channels affected by PFB scalloping. The only inputs are the 
; source name (to select the correct scans) and the name of the output file. Note that the same channels are 
; blanked each time regardless of rest or central frequency. This script should be run AFTER FixAntennaPositions.py AND PAF_edgeoffkeep
; as the former script updates the header with the correct rest frequency and beam offsets and latter calibrates and fits/removes the baseline
; sructure. 
; INPUTS: sname => (string) source name (e.g. NGC6946)
;         fileout => (string) name of output file containing calibrated spectra
; EXAMPLE: PAF_chanBlank, 'NGC6946', fileout='AGBT16B_400_12_NGC6946_chanBlank.fits' 
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
           data = getdata()
           for i=0, 98 do begin
               strtChan = 28+i*32
               endChan = strtChan + 6
               data[strtChan:endChan] = !Values.F_NAN
           ENDFOR
           setdata, data
           keep
        ENDFOR
    ENDFOR
ENDFOR
unfreeze
END


