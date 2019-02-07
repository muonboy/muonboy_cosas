module histogram

contains

subroutine histograma(n,x,nb, xb,yb)
implicit none

real, dimension(n), intent(in) :: x
real, dimension(n), intent(out) :: xb, yb
real, dimension(n) :: yc
real :: xmin, xmax, delta
integer :: nb, ni, i, n


xmin = x(1)
do i = 2, n
	if(x(i).lt.xmin) xmin = x(i)
enddo

xmax = x(1)
do i = 2, n
	if(x(i).gt.xmax) xmax = x(i)
enddo

delta = (xmax-xmin)/nb

yc = 0.0
do i=1, n
	ni = int((x(i)-xmin)/delta)+1   !int() toma al entero mas proximo
	if(ni.eq.(nb+1)) ni=nb
	yc(ni) = yc(ni)+1
enddo

do i = 1, nb
	xb(i) = xmin + (i-0.5)*delta
	yb(i) = yc(i)/(n)               !yb(i)=yc(i)/(n*delta)
enddo

end subroutine

end module
