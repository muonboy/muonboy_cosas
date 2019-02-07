close,/all
n=200

x=dblarr(n)
y=dblarr(n)
work=dblarr(2)

openr, 5,'mh_1d_posterior.dat'

for i=0,n-1 do begin
	readf,5, work	
	x(i)=work(0)
	y(i)=work(1)
endfor

device, dec=0
loadct,5
plot,x,y, back=255, col=0

close,5
end
