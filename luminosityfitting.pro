function lxestimate,cluster,obsid,ra,dec

;; Open Source List
    sourcelist = read_csv('/lustre/scratch/astro/pr83/sourcelists/'+obsid+'.csv')
    sourcenumber = numlines('/lustre/scratch/astro/pr83/sourcelists/'+obsid+'.csv')

;; Open image fits file and header
   image = readfits('/lustre/scratch/astro/pr83/'+obsid+'/sourcelists/'+obsid+'_pn_exp1_bgsub.fits')
   exposure = readfits('/lustre/scratch/astro/pr83/'+obsid+'/images/'+obsid+'-0.50-2.00keV-pn_merged_expmap.fits')
   ignore = readfits('/lustre/scratch/astro/pr83/'+obsid+'/images/'+obsid+'-0.50-2.00keVmerged_img.fits', hdr)

image[where(image lt 0)] = 0

;; Import arcsecond to kpc 
;;   openr, 1, "/lustre/scratch/astro/sb765/luminfit/"+cluster+"/arcsec_to_kpc.dat"
;;   readf, 1, kpcarc

extast, hdr, astr

;; Convert ra,dec to x,y
   ad2xy, ra, dec, astr, x, y

;; Find closest source to cluster
  mindistance = 999999
  for source = 0, sourcenumber-1 do begin
   	distance = sqrt((sourcelist.(0)(source)-x)^2 + (sourcelist.(1)(source)-y)^2)
    	if distance lt mindistance then begin
	     	sourceval = source
        mindistance = distance
      endif
    print, sourceval
  endfor
  print, sourceval

;; Set source to that ra, dec
   x = sourcelist.(0)(sourceval) 
   y = sourcelist.(1)(sourceval)
   print, x, y
;; Remove source from source list
   sourcelist.(0)(sourceval) = 3
   sourcelist.(1)(sourceval) = 3
   sourcelist.(2)(sourceval) = 1
   sourcelist.(3)(sourceval) = 1

;; Remove each source from the image
    for i = 0, sourcenumber-1 do begin
        for j = 0, 511 do begin
            for k = 0, 511 do begin
                if (sourcelist.(0)(i)-j)^2 + (sourcelist.(1)(i)-k)^2 lt (max([sourcelist.(2)(i), sourcelist.(3)(i)]))^2 then begin
                    image(j,k) = 0
                endif
            endfor
        endfor
    endfor
    writefits, "temp.fits", image
;; Create a variable to save the counts(radius)
   countsarr = fltarr(27)

;; Calculate each counts(radius)
  for i = 4, 30 do begin
	counts = 0
    for j = 0, 511 do begin
      for k = 0, 511 do begin
        if (((x-j)^2 + (y-k)^2) lt i^2) then begin
          if exposure(j,k) gt 0 then begin
				  counts = counts + image(j,k)/exposure(j,k)
          endif
        endif
		  endfor
	 endfor
  print, i
	countsarr[i-4] = counts
  endfor
  print, countsarr
  print, "here"

;; Open up LCF table
  lcf = fltarr(27)
  openr, 1, "/lustre/scratch/astro/sb765/lcf/"+cluster+"/lcf_code/lcfpn.dat"
  readf, 1, lcf
  close, 1

;; Calculate Lxs

  LxArr = lcf *countsarr

;; Save counts(radius)

;openw, 1, "/lustre/scratch/astro/sb765/luminfit/lx_results/"+string(ra)+"_"+string(dec)+".csv"
;printf, 1, countsarr
;close, 1

;; Save Lxs

openw, 1, "/lustre/scratch/astro/sb765/luminfit/lx_results/"+cluster+"_"+string(ra)+"_"+string(dec)+".csv"
lxsize=size(LxArr)

for i=0, lxsize(1)-1 do	begin
        printf, 1, LxArr[i]
endfor
;; printf, 1, LxArr
close, 1

end


