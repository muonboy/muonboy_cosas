close,/all
n=33

u=dblarr(n,n)
work=dblarr(n)

openr, 5,'solve_u_ex1.dat'

for i=0,n-1 do begin
	readf,5, work	
	u(i,*)=work
endfor

device, dec=0
loadct,5
surface,u, back=255, col=0

close,5
end
