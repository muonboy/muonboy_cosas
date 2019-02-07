program grid

implicit none

integer :: i, j, N, jmax
real*8 :: rjac
real*8 :: r0, theta0, L, dr, dtheta, M
real*8, dimension(8) :: r, theta
real*8, dimension(8,8) ::  a, b, c, d, e, f, u

N=8
L=5.d0
M=2*dacos(-1.d0)
dr=L/float(N-1)
dtheta=M/float(N-1)
r0=0.1d0
theta0=0.d0
jmax=8
rjac=0.9238d0

open(file='test_grid.dat', unit=12, status='unknown')

do i=0, N-1
	r(i+1)=r0+i*dr
	theta(i+1)=theta0+i*dtheta
	write(12,*) r(i+1), theta(i+1)
enddo

close(12)

open(file='theta.dat', unit=10, status='unknown')

do i=1, N
	do j=1, N
		a(i,j)=1.d0
		b(i,j)=1/r(i)
		c(i,j)=1/(r(i))**2
		d(i,j)=0.d0
		e(i,j)=1.d0
		f(i,j)=0.01d0*(10.d0+1.d0*dcos(5.d0*theta(j)))*(r(i)+2.d0)/((r(i)+1.d0)*r(i)**2)
		u(i,j)=(5.d0+3.d0/r(i)+r(i)*dcos(5.d0*theta(j)))/(r(i)+1.d0)
		write(10,*) f(i,j), u(i,j)
	enddo
enddo
close(10)

call sor(a,b,c,d,e,f,u,jmax,rjac)

open(file='solve_u.dat', unit=11, status='unknown')

do i=0, N-1
	do j=0, N-1
		write(11,*) u(i+1,j+1)
	enddo 
enddo

 close(11)



contains

!____________________________________________________________________________________________

subroutine sor(a,b,c,d,e,f,u,jmax,rjac)
implicit none

!Successive overrelaxation solution of equation (19.5.25) with Chebyshev acceleration. a ,
!b , c , d , e , and f are input as the coefficients of the equation, each dimensioned to the
!grid size JMAX × JMAX . u is input as the initial guess to the solution, usually zero, and
!returns with the final value. rjac is input as the spectral radius of the Jacobi iteration,
!or an estimate of it.

integer :: jmax, ipass, j, jsw, l, lsw, n
integer, parameter :: MAXITS=1000
real*8, parameter :: EPS=1.d-5
real*8 :: rjac, anorm, anormf, omega, resid
real*8, dimension(jmax,jmax) :: a, b, c, d, e, f, u

anormf = 0.d0

do j=2, jmax-1
	do l=2, jmax-1
		anormf=anormf+dabs(f(j,l))
	enddo
enddo

omega=1.d0

do n=1, MAXITS
	anorm=0.d0
	jsw=1
	do ipass=1, 2
	lsw=jsw
		do j=2, jmax-1
			do l=lsw+1, jmax-1, 2
				resid=a(j,l)*u(j+1,l)+b(j,l)*u(j-1,l)+c(j,l)*u(j,l+1)&
					+d(j,l)*u(j,l-1)+e(j,l)*u(j,l)-f(j,l)
				anorm=anorm+dabs(resid)
				u(j,l)=u(j,l)-omega*resid/e(j,l)
			enddo
			lsw=3-lsw
		enddo
		jsw=3-jsw
		if(n.eq.1.and.ipass.eq.1) then
			omega=1.d0/(1.d0-0.5d0*rjac**2)
		else
			omega=1.d0/(1.d0-0.25d0*rjac**2*omega)
		endif
	enddo
	if(anorm.lt.EPS*anormf) return
enddo
!stop 'MAXITS exceeded in sor'

end subroutine



end program
