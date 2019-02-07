program example

use partial_differential_eqs

implicit none

integer, parameter :: N=33
integer :: i, j
real*8 :: x0, y0, L, M, dx, dy
real*8, dimension(N) :: x, y
real*8, dimension(N,N) :: u

L=10.d0
M=10.d0
dx=L/float(N-1)
dy=M/float(N-1)
x0=-5.d0
y0=-5.d0

open(file='grid_ex1.dat', unit=12, status='unknown')
write(11,*) dx, dy

do i=1, N
	x(i)=x0+i*dx
	do j=1, N
		y(j)=y0+j*dy
			if((i.eq.10).and.(j.eq.10)) then
				u(i,j)=1.d0
			else
				u(i,j)=0.d0
			endif

		write(12,*) x(i), y(j), u(i,j)
	enddo
enddo

close(12)

call mglin(u,33,4)

open(file='solve_u_ex1.dat', unit=11, status='unknown')

do i=1, N
	write(11,*) u(i,:)
enddo

 close(11)

end program
