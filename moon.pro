close, /all

n=100000

t = dblarr(n)

rx = dblarr(n)
ry = dblarr(n)
u = dblarr(n)

work = dblarr(4)

openr, 2, 'pos_moon.dat'
 for i = 0, n-1 do begin
  readf,2,work
   t(i) = work(0)
   rx(i) = work(1)
   ry(i) = work(2)
   u(i) = work(3)
 endfor
close, 2

device, dec=0
loadct, 5  ;5 es el estilo de tabla
;plot, rx, ry, back=255, col=50, xs=1, ys=1, xrange=[-0.1,0.1], yrange=[-0.1,0.1]
plot, t, u, back=255, col=50, xs=1, ys=1

 end
