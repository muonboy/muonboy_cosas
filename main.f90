program principal

use statistics
implicit none

integer :: i, mwt
integer, parameter :: ndata=12
real :: sigmax, sigmay, a, b, siga, sigb, chi2, q, adeviation
real, dimension(ndata) :: x, y, sigx, sigy

open(file='datos.dat', unit=11, status='unknown')

do i=1, ndata
	read(11,*) x(i), y(i)
enddo

close(11)

mwt=1
adeviation=1.

call dev(y, ndata, sigmay)
call adev(y, ndata, adeviation)
call dev(x, ndata, sigmax)}

sigx = sigmax
sigy= sigmay

call fit(x,y,ndata,sigy,mwt,a,b,siga,sigb,chi2,q)
!call fit2(x,y,ndata,sigx,sigy,mwt,a,b,siga,sigb,chi2,q)

open(file='outputfit.dat', unit=12, status='unknown')

do i=1, ndata
	write(12,*) x(i), a+b*x(i)
enddo

close(12)

call medfit(x,y,ndata,a,b,adeviation)

open(file='outputmedfit.dat', unit=13, status='unknown')

do i=1, ndata
	write(13,*) x(i), a+b*x(i)
enddo

close(13)

!print*, 'sigmax, sigmay: ', sigmax, sigmay
!print*, 'siga, sigb: ', siga, sigb
!print*, 'chi2: ', chi2
!print*, 'q: ', q

if(q.gt.0.1) print*, 'Modelo confiable'
if((q.gt.0.01).and.(q.lt.0.1)) print*, 'Modelo aceptable'
if(q.lt.0.01) print*,  'Modelo no valido'

end program
