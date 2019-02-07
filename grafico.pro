close,/all
n=8

x=dblarr(n)
y=dblarr(n)
work=dblarr(2)

openr, 5,'resultados_histograma_M_H.dat'

for i=0,n-1 do begin
	readf,5, work	
	x(i)=work(0)
	y(i)=work(1)
endfor

device, dec=0
loadct,5
plot,x,y, back=255, col=0, ps=10, xs=1, ys=1

close,5
end
