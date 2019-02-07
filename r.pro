close, /all

n=100000

t = dblarr(n)

rsx = dblarr(n)
rsy = dblarr(n)

rpx = dblarr(n)
rpy = dblarr(n)

rmx = dblarr(n)
rmy = dblarr(n)

work = dblarr(7)

openr, 2, 'posiciones.dat'
 for i = 0, n-1 do begin
  readf,2,work
   t(i) = work(0)
   rsx(i) = work(1)
   rsy(i) = work(2)
   rpx(i) = work(3)
   rpy(i) = work(4)
   rmx(i) = work(5)
   rmy(i) = work(6)
 endfor
close, 2

device, dec=0
loadct, 5  ;5 es el estilo de tabla
plot, rpx, rpy, back=255, col=50, xs=1, ys=1, xrange=[-2.0,2.0], yrange=[-2,2]
;ps es tipo de puntos, back es el color de background, col es color
; xs y ys ajusta la tabla a los valores de x y y
oplot, rsx, rsy, ps=0, thick=2, col=150
oplot, rmx, rmy, ps=0, thick=2, col=1
;oplot, t, M3, ps=0, thick=2, col=75

 end
