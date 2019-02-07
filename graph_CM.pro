close, /all

n=3000

t = dblarr(n)
E_l = dblarr(n)

work = dblarr(2)

openr, 2, 'energy_CM.dat'
 for i = 0, n-1 do begin
  readf,2,work
   t(i) = work(0)
   E_l(i) = work(1)
 endfor
close, 2

device, dec=0
loadct, 5  ;5 es el estilo de tabla
;set_plot, 'z'
;p = PLOT(t, M1, COLOR='blue', NAME='M1', xtitle='tiempo', ytitle='Masa', THICK=2)
plot, t, E_l, back=255, thick=2, col=1, xs=1, ys=1, xtitle = 'Tiempo', ytitle='E'
;ps es tipo de puntos, back es el color de background, col es color
; xs y ys ajusta la tabla a los valores de x y y
;p2 = PLOT(t, M2, COLOR='ORANGE', NAME='M2', THICK=2, /OVERPLOT)
;p3 = PLOT(t, M3, COLOR='GREEN', NAME='M3', THICK=2, /OVERPLOT)
;l = LEGEND(TARGET=[p,p2,p3], POSITION=[0.75, 0.75, 0.75], /AUTO_TEXT_COLOR)
;oplot, t, E_l, ps=0, thick=2, col=125 ;naranja
;oplot, t, M3, ps=0, thick=2, col=50 ;azul
;write_png,'M1_20.png',tvrd(/true)
;read_png,'M1_20.png'
 end
