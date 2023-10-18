
pro  e
; program to plot 1D images of 
;density
;energy density of electrons and plasmons


Nx=280
FileNumber =458

N_T   = fltarr(FileNumber)
Ee_T = fltarr(FileNumber)
Ew_T  =  fltarr(FileNumber)

DEVICE, RETAIN=2,  DECOMPOSED =0
LOADCT, 39  ;rainbow +white
;myMPEG = OBJ_NEW('IDLgrMPEG',  FILENAME='FW.mpeg') 
;myMPEG = OBJ_NEW('IDLgrMPEG', FILENAME ='ne.mpeg')

for i=0, FileNumber-1 do begin
data_file    = 'fwd' +string(format='(I5.5)', i,/print)+'.dat'

A = fltarr(9,Nx)
; describtion of the file content

print,  data_file
openr,1, data_file
READF, 1,  A
CLOSE, 1

A = transpose(A)
; conversion of data
x  = A(*,0)
N  = A(*,1)
E_b  = A(*,2)
E_w = A(*,3)

window,0, xsize=560, ysize=280
plot,x , N, xrange=[-4,30], xticks=4,Thick=2,yrange=[0,1.0],BACKGROUND=255, COLOR=0, XTITLE ='distance x/d', YTITLE ='n(x)'
image1=TVRD()
window,1, xsize=560, ysize=280
plot,x , E_b, xrange=[-4,30],xticks=4,thick=2,yrange=[0,1.0],BACKGROUND=255, COLOR=8, XTITLE ='velocity', YTITLE ='E(x)'
oplot,x , E_w, thick=2,COLOR=50

image2 =TVRD()
window, 0, xsize=560, ysize=560
TV, image1,0
TV, image2,1
XYOUTS ,500,275,COLOR=0,  string(format='(F4.1, A)', float(i)/10, ' sec', /print ), /device

N_T (i)  = INT_TABULATED( X, N) 
Ee_T (i) = INT_TABULATED( X, E_b)
Ew_T (i) = INT_TABULATED( X, E_w)
end

window,3
PLOT, N_T, Xrange=[0, FileNumber],yrange=[0,2]
oplot, Ee_T/Ee_T(0),line=1
oplot, Ew_T/Ee_T(0),line=2
oplot, (Ew_T+Ee_T)/Ee_T(0), line=3

stop

image = TVRD()
my_image =OBJ_NEW('IDLgrImage', image)
mypalette =OBJ_NEW('IDLgrPalette')
mypalette-> LoadCT, 39
my_image -> SetProperty, Palette=mypalette
write_gif, 'ne.gif', image, /multiple

myMPEG -> Put, my_image, i
myMpeg -> Save

;window, xsize=300, ysize=300
;plot, F_v, vel, BACKGROUND=255, COLOR=0, XTITLE ='time, t', YTITLE ='MAX. VELOCITY, U(t)'


;window,0, xsize=400,ysize=400, title='electron density'
;plot,  distance, electron_density, THICK =5, XTITLE = 'distance x/d', $
;YTITLE = 'Electron density, n(x)/nb', YRANGE=[0,1], XRANGE=[-4,160],$
;BACKGROUND=255, COLOR=0
;grabed_image1 = TVRD()

;IMAGE_CONT_MY, Fxv, /interp
;grabed_image2= TVRD()

;IMAGE_CONT_MY, Wxv, /interp
;grabed_image3 = TVRD()
;window,0, xsize=600, ysize=600
;TV, grabed_image2,0
;TV, grabed_image3,1
;times = string(format='(F4.1, A)', float(i)/10,' sec')
;XYOUTS, 290, 585, times, /device
;grabed_fw = TVRD()

;my_image =OBJ_NEW('IDLgrImage', grabed_fw)
;mypalette =OBJ_NEW('IDLgrPalette')
;mypalette-> LoadCT, 39
;my_image -> SetProperty, Palette=mypalette
;write_gif, 'fw_gif.gif', grabed_FW,  /multiple

;myMPEG -> Put, my_image, i





;myMPEG -> Save
;all files saved to disk
PRINT, 'data saved to tile !'
end
