
pro  velocity
; program to plot 1D images of 
;density
;energy density of electrons and plasmons

FileNumber =3
Nx =217

DEVICE, RETAIN=2,  DECOMPOSED =0
LOADCT, 39  ;rainbow +white

;myMPEG = OBJ_NEW('IDLgrMPEG',  FILENAME='FW.mpeg') 
openw,33, 'speed.dat'



for i=0, FileNumber do begin

data_file    = 'fwd' +string(format='(I5.5)', i,/print)+'.dat'
data_file_f  = 'fxv' +string(format='(I5.5)', i,/print)+'.dat'

print,  data_file
openr,1, data_file
openr,2, data_file_f
;openr,3, data_file_w

; open data file
;openr,2, 'wxv00100.dat'
;B =fltarr(40,1300)

A = fltarr(4,Nx)
Fxv =fltarr(40,Nx)
;Wxv=fltarr(40,1300)

; describtion of the file content
READF, 1,  A
READF, 2, Fxv
;READF, 3, Wxv

CLOSE, 1
CLOSE, 2
;CLOSE, 3

;reading data
;READF, 2,  B
A = transpose(A)
; conversion of data
distance               = A(*,0)
electron_density  = A(*,1)
wave_energy        = A(*,3)
B=transpose(Fxv)
f_x=B(*,5)
print, MAX(f_x, k), k

j=5
WHILE  ((fxv(j,k) GE fxv(5,k)*0.9) AND (j LE 39) ) do j=j+1

step =j

print, step

dens =0.
j=2
WHILE  ( j LE (Nx-1) ) do begin
dens =dens + (distance(j)-distance(j-1))*Electron_density(j)
j=j+1
end

printf, 33,  float(i)/1,  distance(k), (step)*1.1/40., dens
end

close,33
;myMPEG -> Save
;all files saved to disk
PRINT, 'data saved to file !'
end
