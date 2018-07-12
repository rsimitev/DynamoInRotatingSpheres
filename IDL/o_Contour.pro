;####################################################################
;
; program OldContour.pro
; -=-=-=-=-=-=-=-=-=-=-=-
; 
;  Contour plots of velocity and magnetic field in the 
;                    OUTER CORE ONLY
;                    ---------------
;  This program is based on earlier versions of IDL programs that 
;  do basically the same: create contour plots of different observables.
;  This program is designed to work togather with:
;    1) DriveOld.pl  - a driver perl script 
;    2) f_OldContour - a FORTRAN executable that calculates the data to be
;                      ploted
;
;
;  09.06.2001  Radostin Simitev
;
; ####################################################################


; ----- pro make_levels ---------------------------------

pro make_levels, ma,mi,v,vl,vls,vlsC,levelnum

  delta=(ma-mi)/(levelnum-1)
  v=findgen(levelnum)*delta + ceil(mi/delta)*delta
  for i=0,levelnum-1 do begin
  
;    if(v(i) lt 0.0) then vls(i)=2
;    if(v(i) eq 0.0) then vls(i)=1
;    if(v(i) gt 0.0) then vls(i)=0

  ;---Radostin
  ; black&white
  ;
  ;  if(v(i) lt 0.0) then vlsC(i)=5
  ;  if(v(i) eq 0.0) then vlsC(i)=1
  ;  if(v(i) gt 0.0) then vlsC(i)=0

    if(v(i) lt 0.0) then vlsC(i)=0
    if(v(i) eq 0.0) then vlsC(i)=0
    if(v(i) gt 0.0) then vlsC(i)=0
 
  ;  color
  ; 
   if(v(i) lt 0.0) then vls(i)=4
;white - when too much bad resolution
;   if(v(i) eq 0.0) then vls(i)=1
;green
   if(v(i) eq 0.0) then vls(i)=3
   if(v(i) gt 0.0) then vls(i)=2

  endfor
end

; ----- pro circle ---------------------------------------

pro circle, radius, points, plot_color, clipvec

  phi=findgen(points)
  phi=phi/points*2*3.141592654
  xcirc=radius*cos(phi)
  ycirc=radius*sin(phi)
  plots, xcirc(points-1), ycirc(points-1), color=plot_color, clip=clipvec, noclip=0
  plots, xcirc, ycirc, color=plot_color, clip=clipvec, noclip=0, /CONTINUE

end

; ----- pro make_cut ------------------------------------

pro make_cut, cut_type, cutposition, part, field, pipe

  printf, pipe, fix(cut_type)
  printf, pipe, format='(G13.6)', float(cutposition)
  printf, pipe, fix(part)
  printf, pipe, fix(field)

end

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

; ===== main program ====================================

pass=0

;--- some variables
drift=0.0

;--- radii of circles
ri=eta/(1.0-eta)
ro=1.0/(1.0-eta)

case fo of
  1: begin
       set_plot, 'X'
       window, retain=2
       erase, color=-1
     end
  2: begin
     set_plot, 'ps', /copy, /interpolate 
     device, file='plot.ps', /color, bits_per_pixel=8

;       set_plot, 'PS'
;       device,filename='plot.ps',/color

;----------- Radostin
 ;       device,filename='plot.ps',/color,xsize=60,ysize=50
 ;----- sphere - ellipse 
 ;    device,filename='plot.ps',/color,xsize=15,ysize=11
 ;        device,filename='plot.ps',/color,xsize=17,ysize=13
 ;         device,filename='plot.ps',/color,xsize=30,ysize=40
 ;           device,filename='plot.ps',/color,xsize=23,ysize=18

;----------- Radostin
     end
  3: begin
       set_plot,'Z'
       erase, color=-1
     end
  0: stop
endcase

;--- color of plots (0=black)
plot_color=0

circle_points=500

; arrange scales so that circles actually are circular
case fo of
  1: posv=[0.2,0.1,0.8,0.9]
  2: posv=[0.2,0.1,0.767,0.9]
  3: posv=[0.2,0.1,0.767,0.9]
endcase

; plot neither axes nor surr. box:
!x.style=12
!y.style=12

; arrays for the contour defs:
; levelnum=24
v=dblarr(levelnum)
vl=intarr(levelnum)
vls=intarr(levelnum)
vlsC=intarr(levelnum)

if(cut_type le 2) then begin
  if(cut_piece eq 1) then clipvec=[-ro,0,0,ro]
  if(cut_piece eq 2) then clipvec=[0,0,ro,ro]
  if(cut_piece eq 3) then clipvec=[-ro,-ro,0,0]
  if(cut_piece eq 4) then clipvec=[0,-ro,ro,0]
  if(cut_piece eq 5) then clipvec=[-ro,0,ro,ro]
  if(cut_piece eq 6) then clipvec=[-ro,-ro,ro,0]
  if(cut_piece eq 7) then clipvec=[-ro,-ro,0,ro]
  if(cut_piece eq 8) then clipvec=[0,-ro,ro,ro]
  if(cut_piece eq 9) then clipvec=[-ro,-ro,ro,ro]
endif

; =====------------------=====
; ====-  meridional cut  -====
; =====------------------=====

if(cut_type eq 1) then begin

  if(field eq 0) then stop

  field=field + 10*(1 - imean)

; #-#-#-# create plot #-#-#-#

  if(imean eq 1) then begin
    print,'cutting at phi =',phi
  end

  lmax=lmaxf*2+3-1
  ncos=ncosf-1
  z=dblarr(lmax+1,ncos+1)
  y=dblarr(lmax+1,ncos+1)
  x=dblarr(lmax+1,ncos+1)

  openr,3,tmp_x
  for i=0,lmax do for j=0,ncos do begin
    readf,3,tx
    x(i,j)=tx
  end
  close,3

  openr,2,tmp_y
  for i=0,lmax do for j=0,ncos do begin
    readf,2,ty
    y(i,j)=ty
  end
  close,2

  openr,1,tmp_z
  for i=0,lmax do for j=0,ncos do begin
    readf,1,tz
    z(i,j)=tz
  end
  close,1


;define levels for contourplot (including one zero contour!):
;  make_levels, 0.3*abs_max,0.3*abs_min,v,vl,vls,vlsC,levelnum
  make_levels, abs_max,abs_min,v,vl,vls,vlsC,levelnum

  red =   [0,1,1,0,0,1]
  green = [0,1,0,1,0,1]
  blue =  [0,1,0,0,1,0]
  tvlct, 255*red, 255*green, 255*blue

 ;------Radostin
 ; black and white plots
 ;
 ;  contour, z,x,y, /noerase, color=plot_color, clip=clipvec, levels=v, c_linestyle=vls, c_colors=0, c_labels=vl, /follow, pos=posv

 ; color plots 
 ;
; no contour labels
   contour, z,x,y, /noerase, color=plot_color, clip=clipvec, levels=v, c_linestyle=vlsC, c_colors=vls, c_labels=vl, /follow, pos=posv
; contour labels
;   contour, z,x,y, /noerase, color=plot_color, clip=clipvec,   levels=v,  /follow, c_charsize=0.5, c_linestyle=0, c_colors=vls, pos=posv   

; #-#-#-# draw borders #-#-#-#

  circle, ro, circle_points, plot_color, clipvec
  circle, ri, circle_points, plot_color, clipvec

  if((cut_piece eq 1) or (cut_piece eq 3) or (cut_piece eq 5) or (cut_piece eq 6)) then begin
    arrow, -ro,  0, -ri, 0, color=plot_color, /data, hsize=0.
  endif
  if((cut_piece eq 2) or (cut_piece eq 4) or (cut_piece eq 5) or (cut_piece eq 6)) then begin
    arrow, ri,  0, ro, 0, color=plot_color, /data, hsize=0.
  endif
  if((cut_piece eq 1) or (cut_piece eq 2) or (cut_piece eq 7) or (cut_piece eq 8)) then begin
    arrow, 0,  ri, 0, ro, color=plot_color, /data, hsize=0.
  endif
  if((cut_piece eq 3) or (cut_piece eq 4) or (cut_piece eq 7) or (cut_piece eq 8)) then begin
    arrow, 0,  -ri, 0, -ro, color=plot_color, /data, hsize=0.
  endif

  ncos=ncosf
  lmax=lmaxf

endif


; =====------------------=====
; ====-  equatorial cut  -====
; =====------------------=====

if(cut_type eq 2) then begin

  if(field eq 0) then stop

  nangle=m0*nanglef+1-1
  ncos=ncosf-1
  z=dblarr(ncos+1,nangle+1)
  y=dblarr(ncos+1,nangle+1)
  x=dblarr(ncos+1,nangle+1)

  openr,3,tmp_x
  for i=0,ncos do for j=0,nangle do begin
    readf,3,tx
    x(i,j)=tx
  end
  close,3

  openr,2,tmp_y
  for i=0,ncos do for j=0,nangle do begin
    readf,2,ty
    y(i,j)=ty
  end
  close,2

  openr,1,tmp_z
  for i=0,ncos do for j=0,nangle do begin
    readf,1,tz
    z(i,j)=tz
  end
  close,1

;define levels for contourplot (including one zero contour!):
  make_levels, abs_max,abs_min,v,vl,vls,vlsC,levelnum

  red =   [0,1,1,0,0,1]
  green = [0,1,0,1,0,1]
  blue =  [0,1,0,0,1,0]
  tvlct, 255*red, 255*green, 255*blue
 
 ; black&white
 ;
 ; contour, z,x,y, color=plot_color, /noerase, clip=clipvec, levels=v, c_linestyle=vls, c_colors=0, c_labels=vl, /follow, pos=posv
 
 ; color
 ;
   contour, z,x,y, color=plot_color, /noerase, clip=clipvec, levels=v, c_linestyle=vlsC, c_colors=vls, c_labels=vl, /follow, pos=posv

; #-#-#-# draw borders #-#-#-#

  circle, ro, circle_points, plot_color, clipvec
  circle, ri, circle_points, plot_color, clipvec

  if((cut_piece eq 1) or (cut_piece eq 3) or (cut_piece eq 5) or (cut_piece eq 6)) then begin
    arrow, -ro,  0, -ri, 0, color=plot_color, /data, hsize=0.
  endif
  if((cut_piece eq 2) or (cut_piece eq 4) or (cut_piece eq 5) or (cut_piece eq 6)) then begin
    arrow, ri,  0, ro, 0, color=plot_color, /data, hsize=0.
  endif
  if((cut_piece eq 1) or (cut_piece eq 2) or (cut_piece eq 7) or (cut_piece eq 8)) then begin
    arrow, 0,  ri, 0, ro, color=plot_color, /data, hsize=0.
  endif
  if((cut_piece eq 3) or (cut_piece eq 4) or (cut_piece eq 7) or (cut_piece eq 8)) then begin
    arrow, 0,  -ri, 0, -ro, color=plot_color, /data, hsize=0.
  endif

  ncos=ncosf
  nangle=nanglef

endif



; =====-----------------=====
; ====-  spherical cut  -====
; =====-----------------=====

if(cut_type eq 3) then begin

  if(field eq 0) then stop

  if(field ne 9) then begin

  endif else  begin
    r=0.
  endelse

  nangle=m0*nanglef+1-1
  lmax=lmaxf
  z=dblarr(nangle+1,lmax+1)
  lats=dblarr(lmax+1)
  lons=dblarr(nangle+1)

  openr,3,tmp_x
  for i=0,nangle do begin
    readf,3,tx
    lons(i)=tx
  end
  close,3

  openr,2,tmp_y
  for i=0,lmax do begin
    readf,2,ty
    lats(i)=ty
  end
  close,2

  openr,1,tmp_z
  for i=0,nangle do for j=0,lmax do begin
    readf,1,tz
    z(i,j)=tz
  end
  close,1

; sphere
;  map_set, 20, -90, 0, color=plot_color, /noerase, latdel=45, /orthographic, glinestyle=2, glinethick=1, /grid, /noborder, pos=posv

; ellipsoid
  map_set, 0, -90, 0, color=plot_color, /isotropic, /horizon, /noerase, latdel=45, /aitoff, glinestyle=2, glinethick=1, /grid, /noborder, pos=posv

;define levels for contourplot (including one zero contour!):
  make_levels, abs_max,abs_min,v,vl,vls,vlsC,levelnum

  red =   [0,1,1,0,0,1]
  green = [0,1,0,1,0,1]
  blue =  [0,1,0,0,1,0]
  tvlct, 255*red, 255*green, 255*blue

  ;----Radostin
  ;  black and white plots
  ;
  ;  contour, z,lons,lats, /noerase, color=plot_color, levels=v, c_linestyle=vls, c_colors=0,c_labels=vl, /overplot, pos=posv

  ;  color plots
  ;
   contour, z,lons,lats, /noerase, color=plot_color, levels=v, c_linestyle=vlsC, c_colors=vls ,c_labels=vl, /overplot, pos=posv

endif











end





