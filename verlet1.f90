!molecular dinamic--- algoritmo de verlet
module verlet1_

contains
 subroutine pos_aleatorio(npart, ndim , idum, posicion, a, b)
!a, b ::dimensiones de la caja
implicit none
real*8::a,b, l
integer::npart, ndim, idum, j,i
real*8, dimension(ndim, npart)::posicion
l=b-a

do j=1, npart
 do i=1, ndim
    posicion(i,j)= ran0(idum)*l + a
    idum=idum+1
 enddo
enddo

 end subroutine pos_aleatorio

!____________________________________________________________

subroutine inicializar(npart, ndim, posicion, p, acel,a,b)
implicit none
integer :: npart, ndim, i, j, idum
real*8 :: a, b, Temp, sigma, dbulkx, dbulky
real*8, dimension(ndim, npart) :: posicion, p,acel
real*8, dimension(3) :: vbulk
!!!!!!!
idum=132646
call pos_aleatorio(npart, ndim , idum, posicion,a, b)
!!!!!!!

Temp = 300
sigma = 1.38E-23 * Temp
dbulkx = 0.1d0
dbulky = 0.05d0 
vbulk(1) = 1.d0 + dbulkx
vbulk(2) = 1.d0 + dbulky
vbulk(3) = 1.d0

!do j=1, npart
! do i=1, ndim
!    vel(i,j) = dsqrt(1/(3.141592654d0*sigma**2))*dexp(vbulk(i)**2/(2*sigma**2))		!vbulk es vector que depende de i(dimension)
! enddo
!enddo


do j=1, npart
 do i=1, ndim
    p(i,j)=0.d0
 enddo
enddo

do j=1, npart
 do i=1, ndim
    acel(i,j)=0.d0
 enddo
enddo

end subroutine inicializar

real FUNCTION ran2(idum)
IMPLICIT NONE
	INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
	REAL AM,EPS,RNMX
	PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, IR2=3791)
	PARAMETER (NTAB=32, NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
	!Long period (> 2 × 10 18 ) random number generator of L’Ecuyer with Bays-Durham shuffle
	!and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
	!of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
	!alter idum between successive deviates in a sequence. RNMX should approximate the largest
	!floating value that is less than 1.
	INTEGER idum2,j,k,iv(NTAB),iy
	SAVE iv,iy,idum2
	DATA idum2/123456789/, iv/NTAB*0/, iy/0/
	if (idum.le.0) then !Initialize.
	idum=max(-idum,1)	!Be sure to prevent idum = 0.
	idum2=idum
	do  j=NTAB+8,1,-1 !Load the shuffle table (after 8 warm-ups).
	k=idum/IQ1
	idum=IA1*(idum-k*IQ1)-k*IR1
	if (idum.lt.0) idum=idum+IM1
	if (j.le.NTAB) iv(j)=idum
	enddo 
	iy=iv(1)
	endif
	k=idum/IQ1
	idum=IA1*(idum-k*IQ1)-k*IR1
	if (idum.lt.0) idum=idum+IM1
	k=idum2/IQ2
	idum2=IA2*(idum2-k*IQ2)-k*IR2
	if (idum2.lt.0) idum2=idum2+IM2
	j=1+iy/NDIV
	iy=iv(j)-idum2
	iv(j)=idum
	if(iy.lt.1)iy=iy+IMM1
	ran2=min(AM*iy,RNMX)
	
END FUNCTION ran2

real FUNCTION ran0(idum)
implicit none
	INTEGER:: idum,IA,IM,IQ,IR,MASK
	REAL :: AM
	PARAMETER (IA=16807,IM=2147483647,AM=1./IM, IQ=127773,IR=2836,MASK=123459876)
	!“Minimal” random number generator of Park and Miller. Returns a uniform random deviate
	!between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK )
	!to initialize the sequence; idum must not be altered between calls for successive deviates
	!in a sequence.
	INTEGER k
	idum=ieor(idum,MASK)
	!XORing with MASK allows use of zero and other simple bit patterns for idum.
	k=idum/IQ
	idum=IA*(idum-k*IQ)-IR*k
	!Compute idum=mod(IA*idum,IM) without overflows by Schrage’s method.
	if (idum.lt.0) idum=idum+IM
	ran0=AM*idum !Convert idum to a floating result.
	idum=ieor(idum,MASK)  !Unmask before return.
	
END

end module

