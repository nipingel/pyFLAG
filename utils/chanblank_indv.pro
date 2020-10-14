PRO chanblank_indv

;==============================================================================
; Nick Pingel ----- 02/13/20
; Program to blank (i.e., count values are set to NaN) channels affected by PFB scalloping. The inputs are the 
; filename, source name (to select the correct scans) and the name of the output file. Note that the same channels are 
; blanked each time regardless of rest or central frequency. 
; INPUTS: 
;         fileName = > name of file to read into GBTIDL instance 
;         sname => source name (e.g. NGC6946)
;         fileout => name of output file containing calibrated spectra
;=============================================================================== 

;; unpack user arguments from python wrapper
args = COMMAND_LINE_ARGS()
fileName = args[0]
sname = args[1]
outfile = args[2]

;; read in file
filein, fileName

;; select all scans of a source that are mapping scans
allScans=get_scan_numbers(source=sname,procedure='*Map')

;; set name of outfile
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
exit

END


