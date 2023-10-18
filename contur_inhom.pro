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

mx = !d.n_colors-1 ;Brightest color
colors = [mx,mx,mx,0,0,0]	;color vectors
; if !d.name eq 'PS' then colors = mx - colors 
;invert line colors for pstscrp
contour,a,/noerase,/xst,/yst,$	;Do the contour
	   pos = [px[0],py[0], px[0]+swx,py[0]+swy],/dev,$
	c_color =  colors, THICK=2, C_thick=2,yrange=[-4,20],$
     xrange=[0,1.1],nlevels=10, c_labels=[1,1,1,1,1,1,1,1,1],xtitle='velocity v/v0', ytitle='distance x/d'
return
end



pro  contur_inhom
; program to plot 1D images of 
;density
;energy density of electrons and plasmons

Number =1800
;i=70
Ylevels=10
Ymargins=[-4,40]
Ylevels=30
Nx=450
Nv=60
Nk=240

LoadCT,39
DEVICE, RETAIN=2,  DECOMPOSED =0
;FOR i=1, Number  DO BEGIN
;xLOADCT
i=Number

data_file    = 'fwd' +string(format='(I5.5)', i,/print)+'.dat'
data_file_f  = 'fxv' +string(format='(I5.5)', i,/print)+'.dat'
data_file_w = 'wxv' +string(format='(I5.5)', i,/print)+'.dat'
data_file_wf = 'txv' +string(format='(I5.5)', i,/print)+'.dat'

print,  data_file
openr,1, data_file
openr,2, data_file_f
openr,3, data_file_w
openr,4, data_file_wf

;-------------------------------
openr, 5, 'xdf.dat'
XDF = fltarr(5,Nx)
READF,5, XDF
CLOSE, 5

XDF=transpose(XDF)
FREQ=XDF(*,2)
;-------------------------

openr, 10, 'kx_vx.dat'
KV = fltarr(2,Nk)
READF,10, KV
CLOSE, 10

KV=transpose(KV)
Kx=KV(*,0)
Vx=KV(*,1)


; open data file
;openr,2, 'wxv00100.dat'
;B =fltarr(40,1300)

A = fltarr(9,Nx)
Fxv =fltarr(Nv,Nx)
Wxv=fltarr(Nk,Nx)
Txv=fltarr(Nk,Nx)

; describtion of the file content
READF, 1,  A
READF, 2, Fxv
READF, 3, Wxv
READF, 4, Txv

CLOSE, 1
CLOSE, 2
CLOSE, 3
CLOSE, 4

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

;SET_PLOT, 'PS'

;device, /encapsul, finename='fw.ps', xsize=560/100. ,ysize=560/100. 

;IMAGE_CONT_my, Fxv, /interp

window,0, xsize=280, ysize=560

x = 0. + 1.1 *FINDGEN(Nv)/float(Nv)

y = distance

y=freq
; ---------- interpolation to homogeneous grid ---------


FreqHom=Freq(0)-(Freq(0)-Freq(Nx-1))*Transpose(FINDGEN(Nx))/Float(Nx-1)
DistHom =Distance(0)+(Distance(Nx-1)-Distance(0))*FINDGEN(Nx)/Float(Nx-1)


;Wxv=Transpose(Wxv)
FOR j=0, Nk-1 DO BEGIN
Wxv(j,*)=INTERPOL(Wxv(j,*), Distance, DistHom,/spline)
Txv(j,*)=INTERPOL(Txv(j,*), Distance, DistHom,/spline)
END
FOR j=0, Nv-1 DO BEGIN
Fxv(j,*)=INTERPOL(Fxv(j,*), Distance, DistHom,/spline)
END


;contour, fxv,x,FreqHom, BACKGROUND=255,COLOR=1,nlevels=5,xtitle='velocity v/v0',$
; ytitle='Frequency, MHz',  c_labels=[1,1,1,1], xrange=[0,1.1], xticks =6, yticks=Ylevels,$
;C_linestyle =[0,0,0,0], yrange=Ymargins
IMAGE_CONT_MY, Fxv, /interp
grabed_image2= TVRD()

IMAGE_CONT_MY, alog10(1.+Wxv), /interp
grabed_image1= TVRD()

IMAGE_CONT_MY, alog10(1.+Txv), /interp
grabed_image0= TVRD()

;contour,wxv,vx,freqHom, BACKGROUND=255,COLOR=1,nlevels=5,xtitle='wavenumber k*d_De',$
 ;ytitle='Frequency, MHz', c_labels=[1,1,1,1], xrange=[0,1.1], xticks =6, yticks=Ylevels,$
;C_linestyle =[0,0,0,0], yrange=Ymargins
;print, Max(Wxv)


;image_cont, wxv

;TV, image
;DEVICE, /close
;SET_PLOT, 'x'

window,0, xsize=840, ysize=560
TV, grabed_image0,0
TV, grabed_image1,1
TV, grabed_image2,2

;TV, grabed_image4,2
times = string(format='(F4.1, A)', float(i)/100,' s')
XYOUTS, 455, 545, times, /device, color=255

grabed_fw = TVRD()


;write_gif, 'FW_vx.gif', grabed_FW, /multiple

set_plot, 'ps'
DEVICE, /ENCAPSUL, filename='fw.ps', xsize=500/100., ysize=500/100.
TV, grabed_fw
DEVICE, /CLOSE
SET_PLOT, 'x'

;END
print, 'Brghtness Temperature T_F=', Max(Txv)
print, 'Brghtness Temperature T_L=', Max(Wxv)
print, 'OK ----------- !'
end
