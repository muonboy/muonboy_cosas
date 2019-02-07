program prueba_vegas

use random

implicit none

integer :: i,j
integer, parameter :: n=200
real :: dx, dy, pi, intr, intv, sd, chi2a
real, dimension(n) :: x, y
!real, dimension(n,n) :: f
real, dimension(2*n) :: region

pi=acos(-1.0)
sd=0.0
 chi2a=2.0
intr=2.0

x(1)=0
y(1)=0
dx=pi/n
dy=1.5/n
!f(1,1)=0.0
!f(1,2)=1.0

!do i=2, n
!	x(i)=x(i-1)+dx
!	y(i)=y(i-1)+dy
!	f(i,1)=sin(x(i))
!	f(i,2)=1.0	
!enddo

do i=1, 2*n
	if(i.le.n) region(i)=x(i)
	if(i.gt.n) region(i)=y(i-n)
enddo

call vegas(region, n, f, 0, 1000, 5, 0, intv, sd, chi2a)
call vegas(region, n, f, 1, 100000, 1, 0, intv, sd, chi2a)

print*, intv

contains

real function f(x,y)

real :: x, y

f=sin(x)

end function f

end program prueba_vegas
