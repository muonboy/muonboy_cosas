module statistics

contains

!________________________________________________________________________________
subroutine mean(x,N, prom)
implicit none
integer, intent(in) :: N
integer :: i
real, dimension(N) :: x
real, intent(out) :: prom

prom=0.

do i=1, N
	prom=prom+x(i)
enddo

prom=prom/N

end subroutine

!________________________________________________________________________________

subroutine var(x,N, varianza)
implicit none
integer, intent(in) :: N
integer :: i
real, dimension(N) :: x
real, intent(out) :: varianza
real :: media

media=1.0
call mean(x,N, media)
varianza=0.0

do i=1, N
	varianza=varianza + (x(i)-media)**2
enddo

varianza=varianza/(N-1.0)
end subroutine

!_______________________________________________________________________________

subroutine dev(x,N, deviation)
implicit none
integer, intent(in) :: N
real, dimension(N) :: x
real, intent(out) :: deviation
real :: varianza

deviation=1.0
varianza=1.0

call var(x,N, varianza)

deviation=sqrt(varianza)

end subroutine

!______________________________________________________________________________

subroutine adev(x,N, adeviation)
implicit none
integer, intent(in) :: N
integer :: i
real, dimension(N) :: x
real, intent(out) :: adeviation
real :: media

media=1.0
call mean(x,N, media)
adeviation=0.0

do i=1, N
	adeviation = adeviation + abs(x(i)-media)
enddo
adeviation = adeviation/N

end subroutine

!_____________________________________________________________________________

subroutine skew(x,N, skewness)
implicit none
integer, intent(in) :: N
integer :: i
real, dimension(N) :: x
real, intent(out) :: skewness
real :: media, sigma, deviation, cond

deviation=1.0

call dev(x,N, deviation)

cond=sqrt(15.0/N)

if(deviation.gt.cond) then

media=1.0
sigma=1.0

call mean(x,N, media)
call var(x,N, sigma)
skewness=0.0

do i=1, N
	skewness = skewness+((x(i)-media)/sigma)**3.0
enddo

skewness = skewness/N

else
	print*, 'Es innecesario calcular skewness'
endif

end subroutine

!_______________________________________________________________________________

subroutine kurt(x,N, kurtosis)
implicit none
integer, intent(in) :: N
integer :: i
real, dimension(N) :: x
real, intent(out) :: kurtosis
real :: media, sigma, deviation, cond

deviation=1.0

call dev(x,N, deviation)

cond=sqrt(96.0/N)

if(deviation.gt.cond) then

media=1.0
call mean(x,N, media)

do i=1, N
	kurtosis = kurtosis + ((x(i)-media)/sigma)**4.0
enddo

kurtosis = kurtosis/N-3.0
else
	print*, 'No es necesario calcular la kurtosis'
endif
end subroutine

!_______________________________________________________________________________________________
real function gammln(xx)
implicit none

real :: xx
!Returns the value ln[Γ( xx )] for xx > 0.
integer :: j
real*8 :: ser, stp, tmp, x, y
real*8, dimension(6) :: cof
!Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
!accuracy is good enough.

cof(1) = 76.18009172947146d0
cof(2) = -86.50532032941677d0
cof(3) = 24.01409824083091d0
cof(4) = -1.231739572450155d0
cof(5) = .1208650973866179d-2
cof(6) = -.5395239384953d-5

stp = 2.5066282746310005d0

x=xx
y=x
tmp=x+5.5d0
tmp=(x+0.5d0)*log(tmp)-tmp
ser=1.000000000190015d0

do j=1,6
	y=y+1.d0
	ser=ser+cof(j)/y
enddo

gammln=tmp+log(stp*ser/x)

end function gammln

!_____________________________________________________________________________________________

function factrl(n)
implicit none

real :: factrl
!Returns the value n ! as a floating-point number.
integer :: j, ntop, n
real, dimension(33) :: a

ntop=0
a(1)=1.

if (n.lt.0) then
	stop 'negative factorial in factrl'
else if (n.le.ntop) then
	factrl=a(n+1)
else if (n.le.32) then
	do j=ntop+1, n
		a(j+1)=j*a(j)
	enddo
	ntop=n
	factrl=a(n+1)
else
	factrl=exp(gammln(n+1.))
endif

end function factrl

!_____________________________________________________________________________________

function factln(n)
implicit none

integer :: n
real :: factln
real, dimension(100) :: a

a=100*(-1.)

if (n.lt.0) stop 'negative factorial in factln'
if (n.le.99) then
	if (a(n+1).lt.0.) a(n+1)=gammln(n+1.)
	factln=a(n+1)
else
	factln=gammln(n+1.)
endif

end function factln

!___________________________________________________________________________________

function bico(n,k)
implicit none

integer :: k, n
real :: bico

bico = nint(exp(factln(n)-factln(k)-factln(n-k)))

end function bico

!___________________________________________________________________________________

function beta(z,w)
implicit none

real :: beta,w,z

beta=exp(gammln(z)+gammln(w)-gammln(z+w))

end function beta

!___________________________________________________________________________________

subroutine gser(gamser,a,x,gln)
implicit none

integer :: n
integer, parameter :: ITMAX=100
real :: a, gamser, gln, x, ap, del, suma
real, parameter ::EPS=3.e-7

gln=gammln(a)
if(x.le.0.) then
	if(x.lt.0.) stop 'x<0 in gser'
	gamser=0.
endif

ap=a
suma=1./a
del=suma

do n=1, ITMAX
	ap=ap+1.
	del=del*x/ap
	suma=suma+del
	if(abs(del).lt.abs(suma)*EPS) goto 1
enddo

stop 'a too large, ITMAX too small in gser'

1 gamser = suma*exp(-x+a*log(x)-gln)

end subroutine gser

!___________________________________________________________________________________

subroutine gcf(gammcf,a,x,gln)
implicit none
!Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction repre-
!sentation as gammcf . Also returns ln Γ(a) as gln .
!Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accu-
!racy; FPMIN is a number near the smallest representable floating-point number.

integer :: i
integer, parameter :: ITMAX=100
real, parameter :: EPS=3.e-7, FPMIN=1.e-30
real :: a, gammcf, gln, x, an, b, c, d, del, h

gln=gammln(a)
b=x+1.-a
c=1./FPMIN
d=1./b
h=d

do i=1, ITMAX
	an=-i*(i-a)
	b=b+2.
	d=an*d+b
	if(abs(d).lt.FPMIN) d=FPMIN
	c=b+an/c
	if(abs(c).lt.FPMIN) c=FPMIN
	d=1./d
	del=d*c
	h=h*del
	if(abs(del-1.).lt.EPS) goto 1
enddo
stop 'a too large, ITMAX too small in gcf'

1 gammcf=exp(-x+a*log(x)-gln)*h

end subroutine gcf

!________________________________________________________________________________________

real function gammp(a,x)
implicit none
!Returns the incomplete gamma function P (a, x).
real :: a, x, gammcf, gamser, gln

if(x.lt.0..or.a.le.0.) stop 'bad arguments in gammp'
if(x.lt.a+1.) then
	call gser(gamser,a,x,gln)
	gammp=gamser
else
	call gcf(gammcf,a,x,gln)
	gammp=1.-gammcf
endif

end function gammp

!_______________________________________________________________________________________

real function gammq(a,x)
implicit none
!Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x).

real :: a, x, gammcf, gamser, gln

if(x.lt.0..or.a.le.0.) stop 'bad arguments in gammq'
if(x.lt.a+1.) then
	call gser(gamser,a,x,gln)
	gammq=1.-gamser
else
	call gcf(gammcf,a,x,gln)
	gammq=gammcf
endif

end function gammq

!______________________________________________________________________________________

subroutine chstwo(bins1, bins2, nbins, knstrn, df, chsq, prob)
implicit none
integer :: knstrn, nbins, j
real :: chsq, df, prob
real, dimension(nbins) :: bins1, bins2

df=nbins-knstrn
chsq=0.
do j=1, nbins
	if(bins1(j).eq.0..and.bins2(j).eq.0.) then
		df=df-1.
	else
		chsq=chsq+(bins1(j)-bins2(j))**2/(bins1(j)+bins2(j))
	endif
enddo

prob=gammq(0.5*df, 0.5*chsq)

end subroutine chstwo

!_____________________________________________________________________________

subroutine fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)

implicit none

!Given a set of data points x(1:ndata) , y(1:ndata) with individual standard deviations
!sig(1:ndata) , fit them to a straight line y = a + bx by minimizing χ 2 . Returned are
!a,b and their respective probable uncertainties siga and sigb , the chi-square chi2 , and
!the goodness-of-fit probability q (that the fit would have χ 2 this large or larger). If mwt=0
!on input, then the standard deviations are assumed to be unavailable: q is returned as 1.0
!and the normalization of chi2 is to unit standard deviation on all points.

integer :: mwt, ndata, i
real :: a, b, chi2, q, siga, sigb, sigdat, ss, st2, sx, sxoss, sy, t, wt
real, dimension(ndata) :: sig, x, y

sx=0.
sy=0.
st2=0.
b=0.

if(mwt.ne.0) then
	ss=0.
	do i=1, ndata
		wt=1./(sig(i)**2)
		ss=ss+wt
		sx=sx+x(i)*wt
		sy=sy+y(i)*wt
	enddo
else
	do i=1, ndata
		sx=sx+x(i)
		sy=sy+y(i)
	enddo
	ss=float(ndata)
endif

sxoss=sx/ss

if(mwt.ne.0) then
	do i=1, ndata
		t=x(i)-sxoss
		st2=st2+t*t
		b=b+t*y(i)
	enddo
endif

b=b/st2
a=(sy-sx*b)/ss
siga=sqrt((1.+sx*sx/(ss*st2))/ss)
sigb=sqrt(1./st2)
chi2=0.
q=1.

if(mwt.eq.0) then
	do i=1, ndata
		chi2=chi2+(y(i)-a-b*x(i))**2
	enddo

	sigdat=sqrt(chi2/(ndata-2))
	siga=siga*sigdat
	sigb=sigb*sigdat

else
	do i=1, ndata
		chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
	enddo
	
	if(ndata.gt.2) q=gammq(0.5*(ndata-2), 0.5*chi2)
endif

end subroutine fit

!________________________________________________________________

subroutine fit2(x,y,ndata,sigx,sigy,mwt,a,b,siga,sigb,chi2,q)

implicit none

!Given a set of data points x(1:ndata) , y(1:ndata) with individual standard deviations
!sigx(1:ndata) and sigy(1:ndata), fit them to a straight line y = a + bx by minimizing χ 2 . Returned are
!a,b and their respective probable uncertainties siga and sigb , the chi-square chi2 , and
!the goodness-of-fit probability q (that the fit would have χ 2 this large or larger). If mwt=0
!on input, then the standard deviations are assumed to be unavailable: q is returned as 1.0
!and the normalization of chi2 is to unit standard deviation on all points.

integer :: mwt, ndata, i
real :: a, b, chi2, q, siga, sigb, sigdat, ss, st2, sx, sxoss, sy, t, wt, sumw
real, dimension(ndata) :: sigx, sigy, x, y

sumw=0.
sx=0.
sy=0.
st2=0.

call fit(x,y,ndata,sigy,mwt,a,b,siga,sigb,chi2,q)

if(mwt.ne.0) then
	ss=0.

	do i=1, ndata
		wt=1./(sigy(i)**2+(b**2)*(sigx(i))**2)
		sumw=sumw+wt
		ss=ss+wt
		sx=sx+x(i)*wt
		sy=sy+y(i)*wt
	enddo
else
	do i=1, ndata
		sx=sx+x(i)
		sy=sy+y(i)
	enddo
	ss=float(ndata)
endif

sxoss=sx/ss

if(mwt.ne.0) then
	do i=1, ndata
		t=x(i)-sxoss
		st2=st2+t*t
		b=b+t*y(i)
	enddo
endif

b=b/st2
a=(sy-sx*b)/(ss*sumw)
siga=sqrt((1.+sx*sx/(ss*st2))/ss)
sigb=sqrt(1./st2)
chi2=0.
q=1.

if(mwt.eq.0) then
	do i=1, ndata
		chi2=chi2+(y(i)-a-b*x(i))**2
	enddo

	sigdat=sqrt(chi2/(ndata-2))
	siga=siga*sigdat
	sigb=sigb*sigdat

else
	do i=1, ndata
		chi2=chi2+((y(i)-a-b*x(i)))**2/(sigy(i)**2+(b**2)*(sigx(i))**2)
	enddo
	
	if(ndata.gt.2) q=gammq(0.5*(ndata-2), 0.5*chi2)
endif

end subroutine fit2

!_____________________________________________________________________________________

real function selection(k, n, arr)
implicit none

integer :: k, n, i, ir, j, l, mid
real :: a, temp
real, dimension(n) :: arr

!Returns the k th smallest value in the array arr(1:n) . The input array will be rearranged
!to have this value in location arr(k) , with all smaller elements moved to arr(1:k-1) (in
!arbitrary order) and all larger elements in arr[k+1..n] (also in arbitrary order).

l=1
ir=n

do while(ir-l.gt.1)
	mid=(l+ir)/2
	temp=arr(mid)
	arr(mid)=arr(l+1)
	arr(l+1)=temp
	if(arr(l).gt.arr(ir)) then
		temp=arr(l)
		arr(l)=arr(ir)
		arr(ir)=temp
	endif
	if(arr(l+1).gt.arr(ir)) then
		temp=arr(l+1)
		arr(l+1)=arr(i)
		arr(ir)=temp
	endif
	if(arr(l).gt.arr(l+1)) then
		temp=arr(l)
		arr(l)=arr(l+1)
		arr(l+1)=temp
	endif
	i=l+1
	j=ir
	a=arr(l+1)

!____________

	do while(arr(i).lt.a)
		i=i+1
	enddo

	do while(arr(j).gt.a)
		j=j-1
	enddo

!________________
	
	do while(j.ge.i)

	do while(arr(i).lt.a)
		i=i+1
	enddo

	do while(arr(j).gt.a)
		j=j-1
	enddo
		temp=arr(i)
		arr(i)=arr(j)
		arr(j)=temp
	enddo

		arr(l+1)=arr(j)
		arr(j)=a
		if(j.ge.k) ir=j-1
		if(j.le.k) l=1

enddo

	if(ir-l.eq.1) then
		if(arr(ir).lt.arr(l)) then
			temp=arr(l)
			arr(l)=arr(ir)
			arr(ir)=temp
		endif
	endif
	selection=arr(k)

end function selection


!_____________________________________________________________________________________

real function rofunc(b)
implicit none

integer, parameter :: NMAX=1000
real, parameter :: EPS=1.e-7
real :: b
integer :: j, ndata 
real :: aa, abdev, d, suma
real, allocatable, dimension(:) :: arr, x, y

!Evaluates the right-hand side of equation (15.7.16) for a given value of b . Communication
!with the program medfit is through a common block.

allocate(arr(NMAX), x(NMAX), y(NMAX))

do j=1, ndata
	arr(j)=y(j)-b*x(j)
enddo

if(mod(ndata,2).eq.0) then
	j=ndata/2
	aa=0.5*(selection(j,ndata,arr)+selection(j+1,ndata,arr))
else
	aa=selection((ndata+1)/2,ndata,arr)
endif
suma=0.
abdev=0.

do j=1, ndata
	d=y(j)-(b*x(j)+aa)
	abdev=abdev+abs(d)
	if(y(j).ne.0.) d=d/abs(y(j))
	if(abs(d).gt.EPS) suma=suma+x(j)*sign(1.0,d)
enddo

rofunc=suma
deallocate(arr, x, y)

end function rofunc

!_____________________________________________________________________________________

subroutine medfit(x,y,ndata,a,b,abdev)
implicit none

integer :: ndata, ndatat, j
integer, parameter :: NMAX=1000
real :: a, b, abdev, aa, abdevt
real, dimension(ndata) :: x, y
real, allocatable, dimension(:) :: xt, yt, arr
real :: b1, b2, bb, chisq, del , f, f1, f2, sigb, sx, sxx, sxy, sy

allocate(xt(NMAX), yt(NMAX), arr(NMAX))

!Fits y = a + bx by the criterion of least absolute deviations. The arrays x(1:ndata)
!and y(1:ndata) are the input experimental points. The fitted parameters a and b are
!output, along with abdev , which is the mean absolute deviation (in y) of the experimental
!points from the fitted line. This routine uses the routine rofunc , with communication via
!a common block.

sx=0.
sy=0.
sxy=0.
sxx=0.

do j=1, ndata
	xt(j)=x(j)
	yt(j)=y(j)
	sx=sx+x(j)
	sy=sy+y(j)
	sxy=sxy+x(j)*y(j)
	sxx=sxx+x(j)**2
enddo

ndatat=ndata
del=ndata*sxx-sxx**2
aa=(sxx*sy-sx*sxy)/del
bb=(ndata*sxy-sx*sy)/del
chisq=0.

do j=1, ndata
	chisq=chisq+(y(j)-(aa+bb*x(j)))**2
enddo

sigb=sqrt(chisq/del)
b1=bb
f1=rofunc(b1)
b2=bb+sign(3.*sigb,f1)
f2=rofunc(b2)

if(b2.eq.b1) then
	a=aa
	b=bb
	abdev=abdevt/ndata
endif

do while(f1*f2.gt.0.)
	bb=b2+1.6*(b2-b1)
	b1=b2
	f1=f2
	b2=bb
	f2=rofunc(b2)
enddo

sigb=0.01*sigb

do while(abs(b2-b1).gt.sigb)
	bb=b1+0.5*(b2-b1)

	if(bb.eq.b1.or.bb.eq.b2) exit

	f=rofunc(bb)
	if(f*f1.ge.0.) then
		f1=f
		b1=bb
	else
		f2=f
		b2=bb
	endif
enddo

a=aa
b=bb
abdev=abdevt/ndata

deallocate(xt, yt, arr)

end subroutine medfit

end module 
