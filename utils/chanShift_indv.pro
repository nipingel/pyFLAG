PRO chanShift_indv

;==============================================================================
; Nick Pingel ----- 02/13/20
; Program to shift channels channels of the lower LO spectra before combining with a higher Dithered LO spectra. The inputs from the python 
; wrapper are the file name, source name (to select the correct scans), name of the output file, and value of channels to shift.
; Since the LO is dithered by 0.30318/2 = 0.15159 MHz, or 16 fine channels (9.57 kHz wide), and we've smoothed to 24.414 kHz, the number of channels to shift by will
; usually be 16/2.5 = -6.2. It's negative since the lower LO setting is 151.18 kHz LESS than the higher LO setting
; sructure. 
; INPUTS: 
;         fileName = > name of SDFITS file to be read into GBTIDL instance
;         sname => source name (e.g. NGC6946)
;         fileout => name of output file containing calibrated spectra
;         chanShiftVal => number of channels to shift (e.g. -6.2)
;=============================================================================== 

;; unpack user arguments from python wrapper
args = COMMAND_LINE_ARGS()
fileName = args[0]
sname = args[1]
outfile = args[2]
chanShiftVal = FLOAT(args[3])

;; read in file
filein, fileName

;; select all scans of a source that are mapping scans
allScans=get_scan_numbers(source=sname,procedure='*Map')
fileout, outfile

freeze
chan
for idx = 0,N_ELEMENTS(allScans)-1 DO BEGIN
    goodIdx=allScans(idx)

    get, scan=goodIdx, int=0, plnum=1; call scan to get the number of intergrations
    
    ;; get scan info
    info=scan_info(goodIdx)
    nint=info.n_integrations

    ;; start processing for loop
    for int=0,nint-1 do begin 
        for pl=0,1 do begin ;; only care about XX polarization
           get, scan=goodIdx, int=int, plnum=pl
           dcshift, !g.s[0], chanShiftVal
           keep
        ENDFOR
    ENDFOR
ENDFOR
unfreeze

exit
END
