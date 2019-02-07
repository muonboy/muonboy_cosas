program test_algoritmo_M_H

use random
use histogram

implicit none

real, dimension(1000000) :: x
real, allocatable, dimension(:) :: xb, yb
integer :: N, nb

real*8 :: z
integer :: i

! Prueba de posterior

x=-4.0

open(file='mh_1d_posterior.dat', unit=12, status='unknown')

do i = 1, 200
	write(12,*) z, posterior(z)
	z = z + 0.1
enddo

 close(12)

call M_H()

!histograma

N=1000000
nb=8 !Numero de bins

allocate(xb(n), yb(n))

open(file='M-H.dat', unit=10, status='unknown')
	do i=1, N
		read(10,*) x(i)
	enddo
close(10)

call histograma(N,x,nb, xb,yb)

open(file='resultados_histograma_M_H.dat', unit=11, status='unknown')
	do i=1,nb
		write(11,*) xb(i), yb(i)
	enddo
close(11)

deallocate(xb, yb)

end program
