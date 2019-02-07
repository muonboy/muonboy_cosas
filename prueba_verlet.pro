close, /all

n1=200 ;# de datos en el archivo 

work=dblarr(2)
x=dblarr(n1)
y=dblarr(n1)

openr,2, 'prueba_posterior.dat'
   for i=0, n1-1 do begin
       readf,2,work
      x(i)=work(0)
      y(i)=work(1)
   endfor
 close, 2
  

 device, dec=0 ;quitar tabla de colores
loadct,5  ;escojer tabla de colores 
plot, x, y,back=255, col=0, xs=1,ys=1, xr=[-4.0,16.0]
oplot, x, y, thick=2, col=0
end

