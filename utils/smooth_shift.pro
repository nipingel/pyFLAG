pro smooth_shift,sname,outfile,box=box,kernel=kernel,cshift=cshift
; this code will read in a GBTIDL FITS file then smooth the data either 
; using boxcar or gconvol+resample (to get arbitrary resolution), and then
; will shift the data, with gshift, to line up the final channels with other
; datasets.  The user specifies the source, the amount and type of smoothing,
; and the amount of shift required.  
; The kernel should have the format of [0.7,1.,1.,1.,0.7]/sum
; The user needs to specify the source name to select all relevant scans and
; the name of an output file to keep data in.
; 5/31/12 DJP
; gaussian kernel DOES NOT need to be in the form of [0.7,1,1,0.7]/sum. This normalization now takes place in the IF statement 
; 6/30/16 NMP
allscans=get_scan_numbers(source=sname,procedure='*Map')

fileout,outfile

freeze
for s=0,n_elements(allscans)-1 do begin
    info=scan_info(allscans(s))
    nint=info.n_integrations
    for i=0,nint-1 do begin
	for p=0,1 do begin
	    print, allscans(s),i,p
	    get,scan=allscans(s),int=i,plnum=p
	    if keyword_set(box) then begin
		boxcar,box,/decimate 
	    endif
	    if keyword_set(kernel)then begin
		gconvol,kernel/total(kernel)
		resample,!g.s[0].frequency_interval*total(kernel)
	    endif
	    if keyword_set(cshift) then begin
		gshift,cshift,/spline
	    endif
	    keep
	endfor
   endfor
endfor

unfreeze

end
