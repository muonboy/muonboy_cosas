program quick_sort

use random

implicit none

integer, dimension(7) :: S
integer :: i

call sort(S,7)

do i=1, 7
	print*, S(i)
enddo

contains

subroutine sort(S, ndim)
implicit none

integer :: n, idum, flag, i, j, ndim
integer, dimension(n) :: S
integer :: mayor, menor, xm, ymayor, ymenor
real :: nx
integer, allocatable, dimension(:) :: S1
integer, allocatable, dimension(:) :: S2

menor=0
mayor=0
flag=0

do i=1, ndim
	do j=i, ndim
		if(S(i).gt.S(j)) flag=flag+1
	enddo
enddo

if(flag.ne.0) then 
	idum=20

	nx=rand2(idum)*ndim

	n=floor(nx)+1

	xm=S(n)

	do i=1, ndim
		if(i.ne.n) then
			if(x(i).lt.xm) menor=menor+1
			if(x(i).gt.xm) mayor=mayor+1	
		endif
	enddo

	allocate(S1(menor))
	allocate(S2(mayor))

	ymenor=1
	ymayor=1

	do i=1, ndim
		if(i.ne.n) then
			if(x(i).lt.xm) then
				S1(ymenor)=x(i)
				ymenor=ymenor+1
			endif
			if(x(i).gt.xm) then
				S2(ymayor)=x(i)
				ymayor=ymayor+1
			endif	
		endif
	enddo

do i=1, ndim
	if(i.le.n) S(i)=S1(i)
	if(i.eq.n) S(i)=S(n)
	if(i.gt.n) S(i)=S2(i-menor-1)
enddo

call sort(S1, menor)
call sort(S2, mayor)

do i=1, ndim
	if(i.le.n) S(i)=S1(i)
	if(i.eq.n) S(i)=S(n)
	if(i.gt.n) S(i)=S2(i-menor-1)
enddo
deallocate(S1, S2)
endif

end subroutine


end program
