module random
 
contains

!_____________________________________________________________________________________________

real function ran0(idum)
implicit none

integer :: idum, k
integer, parameter :: IA=16807, IM=2147483647, IQ=127773, IR=2836, MASK=123459876
real :: AM

!“Minimal” random number generator of Park and Miller. Returns a uniform random deviate
!between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK )
!to initialize the sequence; idum must not be altered between calls for successive deviates
!in a sequence.

AM=1./IM
idum=ieor(idum, MASK)
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
ran0=AM*idum
idum=ieor(idum, MASK)

end function ran0

!____________________________________________________________________________________________

real function ran1(idum)
implicit none

integer :: idum, j, k, iy
integer, parameter :: IA=16807, IM=2147483647, IQ=127773, IR=2836, NTAB=32
integer , parameter :: NDIV=1+(IM-1)/NTAB
real, parameter :: EPS=1.2E-7, RNMX=1.-EPS
real :: AM
integer, dimension(NTAB) :: iv

save iv, iy
data iv /NTAB*0/, iy /0/
AM=1./IM

if(idum.le.0.or.iy.eq.0) then 
	idum=max(-idum,1)
	do j=NTAB+8,1,-1
		k=idum/IQ
		idum=IA*(idum-k*IQ)-IR*k
		if(idum.lt.0) idum=idum+IM
		if(j.le.NTAB) iv(j)=idum
	enddo
	iy=iv(1)
endif

k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if(idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy, RNMX)

end function ran1

!___________________________________________________________________________________________

real function ran2(idum)
implicit none

integer :: idum, idum2, j, k, iy
integer, parameter :: IM1=2147483563, IM2=2147483399, IA1=40014, IA2=40692
integer, parameter :: IQ1=53668, IQ2=52774, IR1=12211, IR2=3791, NTAB=32, IMM1=IM1-1
real, parameter :: NDIV=1+IMM1/NTAB, EPS=1.2E-7, RNMX=1.-EPS
real :: AM
integer, dimension(NTAB) :: iv

save iv, iy, idum2
data idum2/123456789/, iv/NTAB*0/, iy/0/

AM=1./IM1

if(idum.le.0) then
	idum=max(-idum,1)
	idum2=idum
	do j=NTAB+8, 1,-1
		k=idum/IQ1
		idum=IA1*(idum-k*IQ1)-k*IR1
		if(idum.lt.0) idum=idum+IM1
		if(j.le.NTAB) iv(j)=idum
	enddo
	iy=iv(1)
endif
k=idum/IQ1
idum=IA1*(idum-k*IQ1)-k*IR1
if(idum.lt.0) idum=idum+IM1
k=idum2/IQ2
idum2=IA2*(idum2-k*IQ2)-k*IR2
if(idum2.lt.0) idum2=idum2+IM2
j=1+iy/NDIV
iy=iv(j)-idum2
iv(j)=idum
if(iy.lt.1) iy=iy+IMM1
ran2=min(AM*iy, RNMX)

end function ran2

!___________________________________________________________________________________________

real(kind=8) function m_std_normal_dist()
implicit none

real*8 :: half=0.5d0
real :: s=0.449871, t=-0.386595, a=0.19600, b=0.25472, r1=0.27597, r2=0.27846
real :: u, v, x, y, q

!Crea un numero aleatorio

do 
	call random_number(u) 					!Asigna un aleatorio en u con una semilla del sistema
	call random_number(v)					!u es la semilla para un nuevo aleatorio almacenado en v 

x = u - s
y = abs(v) - t

q = x**2 + y*(a*y - b*x)

if (q.lt.r1) exit
if (q.gt.r2) cycle
if (v**2.lt.-4.0*log(u)*u**2) exit

enddo

m_std_normal_dist = v/u 

end function m_std_normal_dist

!___________________________________________________________________________________________

subroutine random_walk(N, x0, xf, D, Drms)
implicit none

integer :: N, i, j
integer :: idum 
real*8 :: Drms, num, D, step, x0, xf

D=0.d0
xf=x0

do i=1, N
	idum=i
	num=ran2(idum)
	if(num.lt.(0.5d0)) then
		step=-1.d0
	else
		step=1.d0
	endif
	xf=xf+step
enddo

D=dabs(xf-x0)

Drms=sqrt(float(N))

end subroutine

!___________________________________________________________________________________________

real function gamdev1(ia, idum)
implicit none

integer :: ia, idum, j
real :: am, e, s, v1, v2, x, y

!Returns a deviate distributed as a gamma distribution of integer order ia , i.e., a waiting
!time to the ia th event in a Poisson process of unit mean, using ran1(idum) as the source
!of uniform deviates.

if(ia.lt.1) stop 'bad argument in gamdev'
if(ia.lt.6) then
	x=1.0
	do j=1, ia
		x=x*ran1(idum)
	enddo
	x=-log(x)
else
1	v1=ran1(idum)
	v2=2.0*ran1(idum)-1.0
	if(v1**2+v2**2.gt.(1.0)) goto 1
	y=v2/v1
	am=ia-1
	s=sqrt(2.0*am+1.0)
	x=s*y+am
	if(x.le.(0.0)) goto 1
	e=(1.0+y**2)*exp(am*log(x/am)-s*y)
	if(ran1(idum).gt.e) goto 1
endif 
gamdev1=x

end function gamdev1

!_____________________________________________________________________________________________________

real function gamdev2(ia, idum)
implicit none

integer :: ia, idum, j
real :: am, e, s, v1, v2, x, y

!Returns a deviate distributed as a gamma distribution of integer order ia , i.e., a waiting
!time to the ia th event in a Poisson process of unit mean, using ran2(idum) as the source
!of uniform deviates.

if(ia.lt.1) stop 'bad argument in gamdev'
if(ia.lt.6) then
	x=1.0
	do j=1, ia
		x=x*ran2(idum)
	enddo
	x=-log(x)
else
1	v1=ran2(idum)
	v2=2.0*ran2(idum)-1.0
	if(v1**2+v2**2.gt.(1.0)) goto 1
	y=v2/v1
	am=ia-1
	s=sqrt(2.0*am+1.0)
	x=s*y+am
	if(x.le.(0.0)) goto 1
	e=(1.0+y**2)*exp(am*log(x/am)-s*y)
	if(ran2(idum).gt.e) goto 1
endif 
gamdev2=x

end function gamdev2

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

!__________________________________________________________________________________________________

real function poidev(xm, idum)
implicit none

integer :: idum
real :: xm, pi
real :: alxm, em, g, oldm, sq, t, y
save alxm, g, oldm, sq
data oldm /-1./

pi=acos(-1.0) 

!Returns as a floating-point number an integer value that is a random deviate drawn from a
!Poisson distribution of mean xm , using ran1(idum) as a source of uniform random deviates.

if(xm.lt.(12.0)) then
	if(xm.ne.oldm) then
		oldm=xm
		g=exp(-xm)
	endif
	em=-1
	t=1.0
2	em=em+1.0
	t=t*ran1(idum)
	if(t.gt.g) goto 2
else
	if(xm.ne.oldm) then
		oldm=xm
		sq=sqrt(2.0*xm)
		alxm=log(xm)
		g=xm*alxm-gammln(xm+1.0)
	endif
1	y=tan(pi*ran1(idum))
	em=sq*y+xm
	if(em.lt.(0.0)) goto 1
	em=int(em)
	t=0.9*(1.0+y**2)*exp(em*alxm-gammln(em+1.0)-g)
	if(ran1(idum).gt.t) goto 1
end if
poidev=em
end function poidev

!_________________________________________________________________________________________________

subroutine sobseq(n,x)
implicit none
integer :: n, im, ini, ipp, j, k, l, i 
integer, parameter :: MAXBIT=30, MAXDIM=6
real :: fac
real :: x(*)
integer, dimension(MAXDIM) :: ip, ix, mdeg
integer, dimension(MAXBIT*MAXDIM) :: iv
integer, dimension(MAXDIM, MAXBIT) :: iu

!When n is negative, internally initializes a set of MAXBIT direction numbers for each of
!MAXDIM different Sobol’ sequences. When n is positive (but ≤ MAXDIM ), returns as the vector 
!x(1..n) the next values from n of these sequences. ( n must not be changed between
!initializations.)

save ip, mdeg, ix, iv, ini, fac
equivalence(iv,iu)
DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/

if (n.lt.0) then 						!Initialize, don’t return a vector.
	do k=1,MAXDIM
		ix(k)=0
	enddo

	ini=0

	if(iv(1).ne.1) return
		fac=1./2.**MAXBIT
		do k=1,MAXDIM
			do j=1,mdeg(k)				!Stored values only require normalization.
				iu(k,j)=iu(k,j)*2**(MAXBIT-j)
			enddo
			do j=mdeg(k)+1,MAXBIT			!Use the recurrence to get other values.
				ipp=ip(k)
				i=iu(k,j-mdeg(k))
				i=ieor(i,i/2**mdeg(k))
				do l=mdeg(k)-1,1,-1
					if(iand(ipp,1).ne.0) i=ieor(i,iu(k,j-l))
					ipp=ipp/2
				enddo
				iu(k,j)=i
			enddo
		enddo
	else 							!Calculate the next vector in the sequence.
		im=ini
		do j=1,MAXBIT					!Find the rightmost zero bit.
			if(iand(im,1).eq.0)goto 1
			im=im/2
		enddo

		stop 'MAXBIT too small in sobseq'

		1	im=(j-1)*MAXDIM

		do k=1,min(n,MAXDIM)				!XOR the appropriate direction number into each com-
			ix(k)=ieor(ix(k),iv(im+k))		!ponent of the vector and convert to a floating
			x(k)=ix(k)*fac				!number.
		enddo
		ini=ini+1					!Increment the counter.
	endif
return
end subroutine

!__________________________________________________________________________________________________________________


subroutine vegas(region, ndim, fxn, init, ncall, itmx, nprn, tgral, sd, chi2a)
implicit none

integer :: init, itmx, ncall, ndim, nprn
integer, parameter :: NDMX=50, MXDIM=10
real :: tgral, chi2a, sd, fxn
real, dimension(2*ndim) :: region
real, parameter :: ALPH=1.5, TIN=1.e-30
external fxn

!Performs Monte Carlo integration of a user-supplied ndim -dimensional function fxn over
!a rectangular volume specified by region , a 2× ndim vector consisting of ndim “lower
!left” coordinates of the region followed by ndim “upper right” coordinates. The integration
!consists of itmx iterations, each with approximately ncall calls to the function. After each
!iteration the grid is refined; more than 5 or 10 iterations are rarely useful. The input flag
!init signals whether this call is a new start, or a subsequent call for additional iterations
!(see comments below). The input flag nprn (normally 0) controls the amount of diagnostic
!output. Returned answers are tgral (the best estimate of the integral), sd (its standard
!deviation), and chi2a (χ 2 per degree of freedom, an indicator of whether consistent results
!are being obtained). See text for further details.

integer :: i, idum, it, j, k, mds, nd, ndo, ng, npg
integer, dimension(MXDIM) :: ia, kg
real :: calls, dv2g, dxg, f, f2, f2b, fb, rc, ti, tsi, wgt, xjac, xn, xnd, xo
real, dimension(NDMX,MXDIM) :: d, di, xi
real, dimension(MXDIM) :: dt, dx, x
real, dimension(NDMX) :: r, xin
real*8 :: schi, si, swgt
common /ranno/ idum
										!Means for random number initialization.
save										!Best make everything static, allowing restarts.
if(init.le.0)then								!Normal entry. Enter here on a cold start.
	mds=1									!Change to mds=0 to disable stratified sampling, i.e., use im-
	ndo=1									!portance sampling only.
	do j=1,ndim
		xi(1,j)=1.0
	enddo
endif

if (init.le.1)then								!Enter here to inherit the grid from a previous call, but not its
	si=0.d0									!answers.
	swgt=0.d0
	schi=0.d0
endif

if (init.le.2) then								!Enter here to inherit the previous grid and its answers 
	nd=NDMX
	ng=1
	if(mds.ne.0)then							!Set up for stratification.
		ng=(ncall/2.0+0.25)**(1.0/ndim)
		mds=1
		if((2*ng-NDMX).ge.0)then
			mds=-1
			npg=ng/NDMX+1
			nd=ng/npg
			ng=npg*nd
		endif
	endif		
	k=ng**ndim	
	npg=max(ncall/k,2)
	calls=float(npg)*float(k)
	dxg=1.0/ng
	dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-1.0)
	xnd=nd
	dxg=dxg*xnd
	xjac=1./calls

	do j=1,ndim
		dx(j)=region(j+ndim)-region(j)
		xjac=xjac*dx(j)
	enddo

	if(nd.ne.ndo)then							!Do binning if necessary.
		do i=1,max(nd,ndo)
			r(i)=1.
		enddo
		do j=1,ndim
			call rebin(ndo/xnd,nd,r,xin,xi(1,j))
		enddo
		ndo=nd
	endif	

	if(nprn.ge.0) then 
		write(*,200) ndim,calls,it,itmx,nprn,ALPH,mds,nd,(j,region(j),j,region(j+ndim),j=1,ndim)
	endif					
endif

do it=1, itmx
!Main iteration loop. Can enter here (init ≥ 3) to do an additional itmx iterations with all
!other parameters unchanged.
	ti=0.0
	tsi=0.0
	do j=1,ndim
		kg(j)=1
		do i=1,nd
			d(i,j)=0.0
			di(i,j)=0.0
		enddo
	enddo
	10 continue
	fb=0.0
	f2b=0.0
	do k=1,npg
		wgt=xjac
		do j=1,ndim
			xn=(kg(j)-ran2(idum))*dxg+1.
			ia(j)=max(min(int(xn),NDMX),1)
			if(ia(j).gt.1)then
				xo=xi(ia(j),j)-xi(ia(j)-1,j)
				rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
			else
				xo=xi(ia(j),j)
				rc=(xn-ia(j))*xo
			endif
			x(j)=region(j)+rc*dx(j)
			wgt=wgt*xo*xnd
		enddo
		f=wgt*fxn(x,wgt)
		f2=f*f
		fb=fb+f
		f2b=f2b+f2
		do j=1,ndim
			di(ia(j),j)=di(ia(j),j)+f
			if(mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
		enddo
	enddo
	f2b=sqrt(f2b*npg)
	f2b=(f2b-fb)*(f2b+fb)
	if (f2b.le.0.) f2b=TIN
		ti=ti+fb
		tsi=tsi+f2b
		if(mds.lt.0)then						!Use stratified sampling.
		do j=1,ndim
			d(ia(j),j)=d(ia(j),j)+f2b
		enddo
	endif
	do k=ndim,1,-1
		kg(k)=mod(kg(k),ng)+1
		if(kg(k).ne.1) goto 10
	enddo
	tsi=tsi*dv2g								!Compute final results for this iteration.
	wgt=1.0/tsi
	si=si+dble(wgt)*dble(ti)
	schi=schi+dble(wgt)*dble(ti)**2
	swgt=swgt+dble(wgt)
	tgral=si/swgt
	chi2a=max((schi-si*tgral)/(it-0.99d0),0.d0)
	sd=sqrt(1.0/swgt)
	tsi=sqrt(tsi)

	if(nprn.ge.0)then
		write(*,201) it,ti,tsi,tgral,sd,chi2a
		if(nprn.ne.0)then
			do j=1,ndim
				write(*,202) j,(xi(i,j),di(i,j),i=1+nprn/2,nd,nprn)
			enddo
		endif
	endif

	do j=1,ndim 								!Refine the grid. Consult references to understand the subtlety
		xo=d(1,j)  							!of this procedure. The refinement is damped, to avoid
		xn=d(2,j)  							!rapid, destabilizing changes, and also compressed in range
		d(1,j)=(xo+xn)/2.0  						!by the exponent ALPH.
		dt(j)=d(1,j)
		do i=2,nd-1
			rc=xo+xn
			xo=xn
			xn=d(i+1,j)
			d(i,j)=(rc+xn)/3.0
			dt(j)=dt(j)+d(i,j)
		enddo 
		d(nd,j)=(xo+xn)/2.0
		dt(j)=dt(j)+d(nd,j)
	enddo

	do j=1,ndim
		rc=0.0
		do i=1,nd
			if(d(i,j).lt.TIN) d(i,j)=TIN
			r(i)=((1.0-d(i,j)/dt(j))/(log(dt(j))-log(d(i,j))))**ALPH
			rc=rc+r(i)
		enddo
		call rebin(rc/xnd,nd,r,xin,xi(1,j))
	enddo
enddo

200 FORMAT(/' input parameters for vegas: ndim=',i3,' ncall=',f8.0 &
/28x,' it=',i5,' itmx=',i5 &
/28x,' nprn=',i3,' alph=',f5.2/28x,' mds=',i3,' nd=',i4 &
/(30x,'xl(',i2,')= ',g11.4,' xu(',i2,')= ',g11.4))
201 FORMAT(/' iteration no.',I3,': ','integral =',g14.7,'+/- ',g9.2 &
/' all iterations: integral =',g14.7,'+/- ',g9.2, &
' chi**2/it''n =',g9.2)
202 FORMAT(/' data for axis ',I2/' X delta i ', &
' x delta i ',' x delta i ', &
/(1x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))
end subroutine


!_______________________________________________________________________________________________________

real function fxn(x,y)
implicit none

real :: x, y

fxn=sin(x)

end function fxn

!_______________________________________________________________________________________________________


subroutine rebin(rc, nd, r, xin, xi)
implicit none

integer :: nd
real :: rc, r(*), xi(*), xin(*)
integer :: i, k
real :: dr, xn, xo

!Utility routine used by vegas , to rebin a vector of densities xi into new bins defined by
!a vector r .

k=0
xo=0.0
dr=0.0
do i=1,nd-1
	do while(rc.gt.dr)
		k=k+1
		dr=dr+r(k)
	enddo

	if(k.gt.1) xo=xi(k-1)
	xn=xi(k)
	dr=dr-rc
	xin(i)=xn-(xn-xo)*dr/r(k)
enddo

do i=1,nd-1
	xi(i)=xin(i)
enddo

xi(nd)=1.0

end subroutine

!________________________________________________________________________________________________________

subroutine miser(func, region, ndim, npts, dith, ave, var)
implicit none

integer :: ndim, npts
integer, parameter :: MNPT=15, MNBS=4*MNPT, MAXD=10, NSTACK=1000
real :: ave, dith, var, func
real, dimension(2*ndim) :: region
real, parameter :: TIN=1.e-30, BIG=1.e30, PFAC=0.1
external func

!Monte Carlo samples a user-supplied ndim -dimensional function func in a rectangular
!volume specified by region , a 2× ndim vector consisting of ndim “lower-left” coordinates
!of the region followed by ndim “upper-right” coordinates. The function is sampled a total
!of npts times, at locations determined by the method of recursive stratified sampling. The
!mean value of the function in the region is returned as ave ; an estimate of the statistical
!uncertainty of ave (square of standard deviation) is returned as var . The input parameter
!dith should normally be set to zero, but can be set to (e.g.) 0.1 if func ’s active region
!falls on the boundary of a power-of-two subdivision of region .
!Parameters: PFAC is the fraction of remaining function evaluations used at each stage to
!explore the variance of func . At least MNPT function evaluations are performed in any
!terminal subregion; a subregion is further bisected only if at least MNBS function evaluations
!are available. MAXD is the largest value of ndim . NSTACK is the total size of the stack.

integer :: iran, j, jb, jstack, n, naddr, np, npre, nptl, nptr, nptt
real :: avel, fracl, fval, rgl, rgm, rgr, s, sigl, siglb, sigr, sigrb, suma,&
	sumb, summ, summ2, varl
real, dimension(MAXD) :: fmaxl, fmaxr, fminl, fminr, pt, rmid
real, dimension(NSTACK) :: stack
real, dimension(9) :: stf

equivalence (stf(1),avel),(stf(2),varl),(stf(3),jb), (stf(4),nptr),(stf(5),naddr),&
	(stf(6),rgl),(stf(7),rgm),(stf(8),rgr),(stf(9),fracl)
save iran
data iran /0/

jstack=0
nptt=npts
1 continue

if (nptt.lt.MNBS) then									!Too few points to bisect; do straight Monte Carlo.
	np=abs(nptt)
	summ=0.0
	summ2=0.0
	do n=1,np
		call ranpt(pt,region,ndim)
		fval=func(pt)
		summ=summ+fval
		summ2=summ2+fval**2
	enddo
	ave=summ/np
	var=max(TIN,(summ2-summ**2/np)/np**2)
else											!Do the preliminary (uniform) sampling.

	npre=max(int(nptt*PFAC),MNPT)
	do j=1,ndim									!Initialize the left and right bounds for each dimension.
		iran=mod(iran*2661+36979,175000)	
		s=sign(dith,float(iran-87500))
		rmid(j)=(0.5+s)*region(j)+(0.5-s)*region(j+ndim)
		fminl(j)=BIG
		fminr(j)=BIG
		fmaxl(j)=-BIG
		fmaxr(j)=-BIG
	enddo

	do n=1,npre									!Loop over the points in the sample.
		call ranpt(pt,region,ndim)
		fval=func(pt)
		do j=1,ndim								!Find the left and right bounds for each dimension.
			if(pt(j).le.rmid(j))then
				fminl(j)=min(fminl(j),fval)
				fmaxl(j)=max(fmaxl(j),fval)
			else
				fminr(j)=min(fminr(j),fval)
				fmaxr(j)=max(fmaxr(j),fval)
			endif
		enddo
	enddo

	sumb=BIG									!Choose which dimension jb to bisect.
	jb=0
	siglb=1.0
	sigrb=1.0
	do j=1,ndim
		if(fmaxl(j).gt.fminl(j).and.fmaxr(j).gt.fminr(j))then
			sigl=max(TIN,(fmaxl(j)-fminl(j))**(2.0/3.0))
			sigr=max(TIN,(fmaxr(j)-fminr(j))**(2.0/3.0))
			suma=sigl+sigr							!Equation (7.8.24), see text.
			if (suma.le.sumb) then
				sumb=suma
				jb=j
				siglb=sigl
				sigrb=sigr
			endif
		endif
	enddo 

	if (jb.eq.0) jb=1+(ndim*iran)/175000						!MNPT may be too small.
		
	rgl=region(jb)									!Apportion the remaining points between left and right.
	rgm=rmid(jb)
	rgr=region(jb+ndim)
	fracl=abs((rgm-rgl)/(rgr-rgl))
	nptl=MNPT+(nptt-npre-2*MNPT)*fracl*siglb/(fracl*siglb+(1.-fracl)*sigrb) 	!Equation (7.8.23).
	nptr=nptt-npre-nptl
	region(jb+ndim)=rgm								!Set region to left.
	naddr=1										!Push the stack.
	
	do j=1,9
		stack(jstack+j)=stf(j)
	enddo

	jstack=jstack+9
	nptt=nptl

	goto 1										!Dispatch recursive call; will return back here eventually.
	
	10 continue

	avel=ave									!Save left estimates on stack variable.
	varl=var
	region(jb)=rgm 									!Set region to right.
	region(jb+ndim)=rgr
	naddr=2										!Push the stack.
		
	do j=1,9
		stack(jstack+j)=stf(j)
	enddo

	jstack=jstack+9
	nptt=nptr
	goto 1										!Dispatch recursive call; will return back here eventually.
	20 continue
	region(jb)=rgl									!Restore region to original value (so that we don’t
	ave=fracl*avel+(1.0-fracl)*ave							!need to include it on the stack).
	var=fracl**2*varl+(1.0-fracl)**2*var						!Combine left and right regions by equa-
endif											!tion (7.8.11)(1st line)

if (jstack.ne.0) then									!Pop the stack.
	jstack=jstack-9
	do j=1,9
		stf(j)=stack(jstack+j)
	enddo
	goto (10,20),naddr
	stop 'miser: never get here'
endif

end subroutine

!_______________________________________________________________________________________________

subroutine ranpt(pt, region, n)
implicit none

integer :: n,idum
real, dimension(n) :: pt
real, dimension(2*n) :: region

common /ranno/ idum
save /ranno/

!Returns a uniformly random point pt in an n -dimensional rectangular region . Used by
!miser ; calls ran1 for uniform deviates. Your main program should initialize idum , through
!the COMMON block /ranno/ , to a negative seed integer.

integer :: j

do j=1,n
	pt(j)=region(j)+(region(j+n)-region(j))*ran1(idum)
enddo

end subroutine

!________________________________________________________________________________________________

real(kind=8) function posterior(x)
implicit none

real*8 :: x
real*8 :: mu, sigma, pi

real*8 :: num, den

!Funcion gaussiana

pi=dacos(-1.d0)
mu=5.d0
sigma=2.d0

num = dexp(-(x-mu)**2/(2.d0*sigma**2))
den = dsqrt(2.d0*pi*sigma**2)

posterior = num/den

end function posterior

!________________________________Algoritmo de Metropolis-Hastings________________________________
!_______________________basado en la funcion de distribucion ¨posterior¨_________________________

subroutine M_H()
implicit none
integer, parameter :: N = 1e6
integer :: i
real*8 :: delta = 0.1d0
real*8 :: x0, u, p, pn, x, xn

print*, 'x0= '
read*, x0

p = posterior(x0)

open(file='M-H.dat', unit=11, status='unknown')

write(11,*) x0

do i = 1, N
	xn = x0 + m_std_normal_dist()
	pn = posterior(xn)
	
	if (pn.ge.p) then
		x = xn
		p = pn
	else
		u = rand(0)
		if(u.lt.pn/p) then
			p = pn
			x = xn
		endif
	endif

	write(11,*) x
enddo

 close(11)

print*, 'Distancia= ', x-x0

end subroutine M_H

!_______________________________________________________________________________________________

end module
