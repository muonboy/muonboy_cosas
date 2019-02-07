program test_random
use random

implicit none

integer :: j, k, i, a
real :: N, P, Q, xm

a=5
xm=5.0

do i=1, 2
	do j=1, 10
		k=j
		P=gamdev1(a, k)
		Q=gamdev2(a, k)
		N=poidev(xm,k)
		print*, j, P, Q, N
	enddo
	print*, '_________________________________'
enddo

end program
