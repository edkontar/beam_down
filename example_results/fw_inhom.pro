;subroutine to plot my style images !!! 


pro image_cont_my, a, WINDOW_SCALE = window_scale, ASPECT = aspect, $
	INTERP = interp
;+
; NAME:
;	IMAGE_CONT
;
; PURPOSE:
;	Overlay an image and a contour plot.
;
; CATEGORY:
;	General graphics.
;
; CALLING SEQUENCE:
;	IMAGE_CONT, A
;
; INPUTS:
;	A:	The two-dimensional array to display.
;
; KEYWORD PARAMETERS:
; WINDOW_SCALE:	Set this keyword to scale the window size to the image size.
;		Otherwise, the image size is scaled to the window size.
;		This keyword is ignored when outputting to devices with 
;		scalable pixels (e.g., PostScript).
;
;	ASPECT:	Set this keyword to retain the image's aspect ratio.
;		Square pixels are assumed.  If WINDOW_SCALE is set, the 
;		aspect ratio is automatically retained.
;
;	INTERP:	If this keyword is set, bilinear interpolation is used if 
;		the image is resized.
;
; OUTPUTS:
;	No explicit outputs.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	The currently selected display is affected.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	If the device has scalable pixels, then the image is written over
;	the plot window.
;
; MODIFICATION HISTORY:
;	DMS, May, 1988.
;-

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz[0] lt 2 then message, 'Parameter not 2D'
window, 0, xsize=280, ysize=560
	;set window used by contour
contour,[[0,0],[1,1]],/nodata, xstyle=4, ystyle = 4

px = !x.window * !d.x_vsize	;Get size of window in device units
py = !y.window * !d.y_vsize
swx = px[1]-px[0]		;Size in x in device units
swy = py[1]-py[0]		;Size in Y
six = float(sz[1])		;Image sizes
siy = float(sz[2])
aspi = six / siy		;Image aspect ratio
aspw = swx / swy		;Window aspect ratio
f = aspi / aspw			;Ratio of aspect ratios

if (!d.flags and 1) ne 0 then begin	;Scalable pixels?
  if keyword_set(aspect) then begin	;Retain aspect ratio?
				;Adjust window size
	if f ge 1.0 then swy = swy / f else swx = swx * f
	endif

  

  tvscl,a,px[0],py[0],xsize = swx, ysize = swy, /device

endif else begin	;Not scalable pixels
   if keyword_set(window_scale) then begin ;Scale window to image?
	tvscl,a,px[0],py[0]	;Output image
	swx = six		;Set window size from image
	swy = siy
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
	tv,poly_2d(bytscl(a),$	;Have to resample image
		[[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
		keyword_set(interp),swx,swy), $
		px[0],py[0]
	endelse			;window_scale
  endelse			;scalable pixels

mx = !d.n_colors-1		;Brightest color
colors = [mx,mx,mx,0,0,0]	;color vectors
if !d.name eq 'PS' then colors = mx - colors ;invert line colors for pstscrp
contour,a,/noerase,/xst,/yst,$	;Do the contour
	   pos = [px[0],py[0], px[0]+swx,py[0]+swy],/dev,$
	c_color =  colors, THICK=0, C_thick=0,yrange=[-4,30],$
     xrange=[0,1.1],LEVELS=[10,11],xtitle='velocity v/v0', ytitle='distance x/d'
return
end



pro  fw
; program to plot 1D images of 
;density
;energy density of electrons and plasmons

FileNumber =50

DEVICE, RETAIN=2,  DECOMPOSED =0
LOADCT, 39  ;rainbow +white

myMPEG = OBJ_NEW('IDLgrMPEG',  FILENAME='FW.mpeg') 

for i=0, FileNumber do begin

data_file    = 'fwd' +string(format='(I5.5)', i,/print)+'.dat'
data_file_f  = 'fxv' +string(format='(I5.5)', i,/print)+'.dat'
data_file_w = 'wxv' +string(format='(I5.5)', i,/print)+'.dat'

print,  data_file
openr,1, data_file
openr,2, data_file_f
openr,3, data_file_w

; open data file
;openr,2, 'wxv00100.dat'
;B =fltarr(40,1300)

A = fltarr(4,250)
Fxv =fltarr(50,250)
Wxv=fltarr(50,250)

; describtion of the file content
READF, 1,  A
READF, 2, Fxv
READF, 3, Wxv

CLOSE, 1
CLOSE, 2
CLOSE, 3

;reading data
;READF, 2,  B
A = transpose(A)
; conversion of data
distance               = A(*,0)
electron_density  = A(*,1)
electron_energy  =  A(*,2)
wave_energy        = A(*,3)
;window,0, xsize=400,ysize=400, title='electron density'
;plot,  distance, electron_density, THICK =5, XTITLE = 'distance x/d', $
;YTITLE = 'Electron density, n(x)/nb', YRANGE=[0,1], XRANGE=[-4,160],$
;BACKGROUND=255, COLOR=0
;grabed_image1 = TVRD()

IMAGE_CONT_MY, Fxv, /interp
grabed_image2= TVRD()

IMAGE_CONT_MY, Wxv, /interp
grabed_image3 = TVRD()
window,0, xsize=560, ysize=560
TV, grabed_image2,0
TV, grabed_image3,1
times = string(format='(F4.1, A)', float(i)/10,' sec')
XYOUTS, 275, 545, times, /device
grabed_fw = TVRD()

my_image =OBJ_NEW('IDLgrImage', grabed_fw)
mypalette =OBJ_NEW('IDLgrPalette')
mypalette-> LoadCT, 39
my_image -> SetProperty, Palette=mypalette
write_gif, 'fw_gif.gif', grabed_FW,  /multiple

myMPEG -> Put, my_image, i

end

myMPEG -> Save
;all files saved to disk
PRINT, 'Animation completed  ---- OK '
end
