module partial_differential_eqs

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
real*8, dimension(jmax, jmax) :: a, b, c, d, e, f, u

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
stop 'MAXITS exceeded in sor'

end subroutine

!___________________________________________________________________________________

subroutine rstrct(uc, uf, nc)
implicit none

integer :: nc, ic, it, jc, jf
real*8, dimension(nc,nc) :: uc
real*8, dimension(2*nc-1,2*nc-1) :: uf

!Half-weighting restriction. nc is the coarse-grid dimension. The fine-grid solution is input
!in uf(1:2*nc-1,1:2*nc-1) , the coarse-grid solution is returned in uc(1:nc,1:nc) .

do jc=2,nc-1
	jf=2*jc-1
	do ic=2, nc-1
		it=2*ic-1
		uc(ic,jc)=0.5d0*uf(it,jf)+0.125d0*(uf(it+1,jf)+&
			  uf(it-1,jf)+uf(it,jf+1)+uf(it,jf-1))
	enddo
enddo

do ic=1,nc
	uc(1,jc)=uf(1,2*jc-1)
	uc(nc,jc)=uf(2*nc-1,2*jc-1)
enddo

end subroutine 

!_______________________________________________________________________________

subroutine interp(uf,uc,nf)
implicit none

integer :: nf, ic, jc, jf, nc, it
real*8, dimension(nf,nf) :: uf
real*8, dimension(nf/2+1, nf/2+1) :: uc

!Coarse-to-fine prolongation by bilinear interpolation. nf is the fine-grid dimension. The
!coarse-grid solution is input as uc(1:nc,1:nc) , where nc = nf /2 + 1. The fine-grid
!solution is returned in uf(1:nf,1:nf) .

nc=nf/2+1

do jc=1,nc
	jf=2*jc-1
	do ic=1,nc
		uf(2*ic-1,jf)=uc(ic,jc)
	enddo
enddo

do jf=1,nf,2
	do it=2, nf-1, 2
		uf(it,jf)=0.5d0*(uf(it+1,jf)+uf(it-1,jf))
	enddo
enddo

do jf=2,nf-1,2
	do it=1, nf
		uf(it,jf)=0.5d0*(uf(it,jf+1)+uf(it,jf-1))
	enddo
enddo
end subroutine

!____________________________________________________________________

subroutine addint(uf, uc,res,nf)
implicit none

integer :: nf, i, j
real*8, dimension(nf,nf) :: res, uf
real*8, dimension(nf/2+1, nf/2+1) :: uc

!Does coarse-to-fine interpolation and adds result to uf . nf is the fine-grid dimension. The
!coarse-grid solution is input as uc(1:nc,1:nc) , where nc = nf /2 + 1. The fine-grid
!solution is returned in uf(1:nf,1:nf) . res(1:nf,1:nf) is used for temporary storage.

call interp(res,uc,nf)

do j=1,nf
	do i=1,nf
		uf(i,j)=uf(i,j)+res(i,j)
	enddo
enddo
end subroutine

!______________________________________________________________________

subroutine slvsml(u,rhs)
implicit none

real*8, dimension(3,3) :: rhs, u
real*8 :: h

!Solution of the model problem on the coarsest grid, where h = 12 . The right-hand side is
!input in rhs(1:3,1:3) and the solution is returned in u(1:3,1:3) .

call fill0(u,3)
h=0.5d0
u(2,2)=-h*h*rhs(2,2)/4.d0

end subroutine

!______________________________________________________________________

subroutine relax(u,rhs,n)
implicit none

integer :: n, i ,ipass, isw, j, jsw
real*8 :: h, h2
real*8, dimension(n,n) :: rhs, u

!Red-black Gauss-Seidel relaxation for model problem. The current value of the solution
!u(1:n,1:n) is updated, using the right-hand side function rhs(1:n,1:n) .

h=1.d0/(n-1)
h2=h*h
jsw=1

do ipass=1,2
	isw=jsw
	do j=2,n-1
		do i=isw+1,n-1,2
			u(i,j)=0.25d0*(u(i+1,j)+u(i-1,j)+u(i,j+1)+&
				u(i,j-1)-h2*rhs(i,j))
		enddo
		isw=3-isw
	enddo
	jsw=3-jsw
enddo

end subroutine

!____________________________________________________________________

subroutine resid(res,u,rhs,n)
implicit none

integer :: n, i, j
real*8 :: h, h2i
real*8, dimension(n,n) :: res, rhs, u

!Returns minus the residual for the model problem. Input quantities are u(1:n,1:n) and
!rhs(1:n,1:n) , while res(1:n,1:n) is returned.

h=1.d0/(n-1)
h2i=1.d0/(h*h)

do j=2,n-1
	do i=2,n-1
		res(i,j)=-h2i*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-&
			 4.d0*u(i,j))+rhs(i,j)
	enddo
enddo

do i=1, n
	res(i,1)=0.d0
	res(i,n)=0.d0
	res(1,i)=0.d0
	res(n,i)=0.d0
enddo

end subroutine

!__________________________________________________________________

subroutine copy(aout, ain, n)
implicit none

integer :: n, i, j
real*8, dimension(n,n) :: ain, aout

!Copies ain(1:n,1:n) to aout(1:n,1:n)

do i=1,n
	do j=1,n
		aout(j,i)=ain(j,i)
	enddo
enddo

end subroutine

!__________________________________________________________________

subroutine fill0(u,n)
implicit none

integer :: n, i, j
real*8, dimension(n,n) :: u

!Fills u(1:n,1:n) with zeros.

do j=1,n
	do i=1,n
		u(i,j)=0.d0
	enddo
enddo

end subroutine

!_________________________________________________________________

integer function maloc(length)
implicit none

integer :: length, mem
real*8 :: z
!For mglin
!integer, parameter :: NG=5, MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-100/3

!For mgfas
integer, parameter :: NG=5, MEMLEN=17*2**(2*NG)/3+18*2**NG+10*NG-86/3

!Dynamical storage allocation. Returns integer pointer to the starting position for len array
!elements in the array z . The preceding array element is filled with the value of len , and
!the variable mem is updated to point to the last element of z that has been used.

common /memory/ z(MEMLEN), mem

if(mem+length+1.gt.MEMLEN) stop 'insufficient memory in maloc'

z(mem+1)=length
maloc=mem+2
mem=mem+length+1

end function maloc

!____________________________________________________________________________________

subroutine mglin(u,n,ncycle)
implicit none

integer, parameter :: NG=5, MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-100/3
integer, parameter :: NPRE=1, NPOST=1
integer :: n, ncycle, j, jcycle, jj, jpost, jpre, nf, ngrid, nn, mem
integer, dimension(NG) :: ires, irho, irhs, iu
real*8, dimension(n,n) :: u
real*8 :: z
common /memory/ z(MEMLEN), mem


!Full Multigrid Algorithm for solution of linear elliptic equation, here the model problem
!(19.0.6). On input u(1:n,1:n) contains the right-hand side ρ, while on output it returns
!the solution. The dimension n is related to the number of grid levels used in the solution,
!NG below, by n = 2** NG + 1. ncycle is the number of V-cycles to be used at each level.
!Parameters: NG is the number of grid levels used; MEMLEN is the maximum amount of
!memory that can be allocated by calls to maloc ; NPRE and NPOST are the number of
!relaxation sweeps before and after the coarse-grid correction is computed.

mem=0
nn=n/2+1
ngrid=NG-1
irho(ngrid)=maloc(nn*2)

call rstrct(z(irho(ngrid)),u,nn)

do while(nn.gt.3)
	nn=nn/2+1
	ngrid=ngrid-1
	irho(ngrid)=maloc(nn**2)
	call rstrct(z(irho(ngrid)), z(irho(ngrid+1)), nn)
enddo

nn=3
iu(1)=maloc(nn**2)
irhs(1)=maloc(nn**2)
call slvsml(z(iu(1)), z(irho(1)))
ngrid=NG

do j=2, ngrid
	nn=2*nn-1
	iu(j)=maloc(nn**2)
	irhs(j)=maloc(nn**2)
	ires(j)=maloc(nn**2)
	call interp(z(iu(j)), z(iu(j-1)),nn)
	if(j.ne.ngrid) then
		call copy(z(irhs(j)), z(irho(j)), nn)
	else
		call copy(z(irhs(j)), u, nn)
	endif
	do jcycle=1, ncycle
		nf=nn
		do jj=j,2,-1
			do jpre=1, NPRE
				call relax(z(iu(jj)), z(irhs(jj)),nf)
			enddo
			call resid(z(ires(jj)), z(iu(jj)), z(irhs(jj)),nf)
			nf=nf/2+1
			call rstrct(z(irhs(jj-1)), z(ires(jj)), nf)
			
			call fill0(z(iu(jj-1)), nf)
		enddo
		call slvsml(z(iu(1)),z(irhs(1)))
		nf=3
		do jj=2,j
			nf=2*nf-1
			call addint(z(iu(jj)),z(iu(jj-1)),z(ires(jj)),nf)

			do jpost=1,NPOST
				call relax(z(iu(jj)), z(irhs(jj)), nf)
			enddo 
		enddo 
	enddo 			
enddo 	

call copy(u,z(iu(ngrid)),n)

end subroutine

!_______________________________________________________________________

subroutine relax2(u,rhs,n)
implicit none

integer :: n, ipass, isw, i, j, jsw
real*8 :: foh2, h, h2i, res
real*8, dimension(n,n) :: rhs, u

!Red-black Gauss-Seidel relaxation for equation (19.6.44). The current value of the solution
!u(1:n,1:n) is updated, using the right-hand side function rhs(1:n,1:n) .

h=1.d0/(n-1)
h2i=1.d0/(h*h)
foh2=-4.d0*h2i
jsw=1

do ipass=1, 2
	isw=jsw
	do j=2, n-1
		do i=isw+1,n-1,2
			res=h2i*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-&
				4.d0*u(i,j))+u(i,j)**2-rhs(i,j)
			u(i,j)=u(i,j)-res/(foh2+2.d0*u(i,j))
		enddo
		isw=3-isw
	enddo
	jsw=3-jsw
enddo

end subroutine

!___________________________________________________________________

subroutine slvsm2(u,rhs)
implicit none

real*8 :: disc, fact, h
real*8, dimension(3,3) :: rhs, u

!Solution of equation (19.6.44) on the coarsest grid, where h = 12 . The right-hand side is
!input in rhs(1:3,1:3) and the solution is returned in u(1:3,1:3) .

call fill0(u,3)
h=0.5d0
fact=2.d0/h**2
disc=dsqrt(fact**2+rhs(2,2))
u(2,2)=-rhs(2,2)/(fact+disc)

end subroutine

!____________________________________________________________________

subroutine lop(out1,u,n)
implicit none

integer :: n, i, j
real*8 :: h, h2i
real*8, dimension(n,n) :: out1, u 

!Given u(1:n,1:n) , returns L h (u h ) for equation (19.6.44) in out(1:n,1:n) .

h=1.d0/(n-1)
h2i=1.d0/(h*h)

do j=2, n-1
	do i=2,n-1
		out1=h2i*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-&
			4.d0*u(i,j))+u(i,j)**2
	enddo
enddo

do i= 1, n
	out1(i,1)=0.d0
	out1(i,n)=0.d0
	out1(1,i)=0.d0
	out1(n,1)=0.d0
enddo
end subroutine

!_________________________________________________________________

subroutine matadd(a,b,c,n)
implicit none

integer :: i, j, n
real*8, dimension(n,n) :: a, b, c

!Adds a(1:n,1:n) to b(1:n,1:n) and returns result in c(1:n,1:n) .

do j=1,n
	do i=1,n
		c(i,j)=a(i,j)+b(i,j)
	enddo
enddo

end subroutine

!_____________________________________________________________

subroutine matsub(a,b,c,n)
implicit none

integer :: n, i, j
real*8, dimension(n,n) :: a, b, c

!Subtracts b(1:n,1:n) from a(1:n,1:n) and returns result in c(1:n,1:n) .

do j=1, n
	do i=1, n
		c(i,j)=a(i,j)-b(i,j)
	enddo
enddo

end subroutine 

!_________________________________________________________________________

real*8 function anorm2(a,n)
implicit none

integer :: n, i, j
real*8 :: suma
real*8, dimension(n,n) :: a

!Returns the Euclidean norm of the matrix a(1:n,1:n) .

suma=0.d0

do j=1,n
	do i=1,n
		suma=suma+a(i,j)**2
	enddo
enddo

anorm2=dsqrt(suma)/n

end function anorm2

!__________________________________________________________

subroutine mgfas(u,n,maxcyc)
implicit none

integer, parameter :: NG=5, MEMLEN=17*2**(2*NG)/3+18*2**NG+10*NG-86/3
integer, parameter :: NPRE=1, NPOST=1
real*8, parameter :: ALPHA=0.33d0
real*8, dimension(n,n) :: u
integer :: maxcyc, n, j, jcycle, jj, jm1, jpost, jpre, mem, nf, ngrid, nn
integer, dimension(NG) :: irho, irhs, itau, itemp, iu
real*8 :: res, trerr, z
common /memory/ z(MEMLEN), mem

!Full Multigrid Algorithm for FAS solution of nonlinear elliptic equation, here equation
!(19.6.44). On input u(1:n,1:n) contains the right-hand side ρ, while on output it re-
!turns the solution. The dimension n is related to the number of grid levels used in the
!solution, NG below, by n = 2** NG + 1. maxcyc is the maximum number of V-cycles to
!be used at each level.
!Parameters: NG is the number of grid levels used; MEMLEN is the maximum amount of
!memory that can be allocated by calls to maloc ; NPRE and NPOST are the number of
!relaxation sweeps before and after the coarse-grid correction is computed; ALPHA relates
!the estimated truncation error to the norm of the residual.

mem=0
nn=n/2+1
ngrid=NG-1
irho(ngrid)=maloc(nn**2)

call rstrct(z(irho(ngrid)), u, nn)

do while(nn.gt.3)
	nn=nn/2+1
	ngrid=ngrid-1
	irho(ngrid)=maloc(nn**2)
	call rstrct(z(irho(ngrid)), z(irho(ngrid+1)), nn)
enddo

nn=3
iu(1)=maloc(nn**2)
irhs(1)=maloc(nn**2)
itau(1)=maloc(nn**2)
itemp(1)=maloc(nn**2)
call slvsm2(z(iu(1)),z(irho(1)))
ngrid=NG

do j=2, ngrid
	nn=2*nn-1
	iu(j)=maloc(nn**2)
	irhs(j)=maloc(nn**2)
	itau(j)=maloc(nn**2)
	itemp(j)=maloc(nn**2)
	call interp(z(iu(j)), z(iu(j-1)), nn)
	if(j.ne.ngrid) then
		call copy(z(irhs(j)), z(irho(j)), nn)
	else
		call copy(z(irhs(j)), u, nn)
	endif

!_______________
	jcycle=1
	do while((res.gt.trerr).and.(jcycle.le.maxcyc))
		nf=nn
		do jj=j,2,-1
			do jpre=1, NPRE
				call relax2(z(iu(jj)), z(irhs(jj)),nf)
			enddo
			call lop(z(itemp(jj)), z(iu(jj)), nf)
			nf=nf/2+1
			jm1=jj-1
			call rstrct(z(itemp(jm1)), z(itemp(jj)), nf)
			call rstrct(z(iu(jm1)), z(iu(jj)), nf)
			call lop(z(itau(jm1)), z(iu(jm1)), nf)
			call matsub(z(itau(jm1)), z(itemp(jm1)), z(itau(jm1)), nf)

			if(jj.eq.j) trerr=ALPHA*anorm2(z(itau(jm1)),nf)
			call rstrct(z(irhs(jm1)), z(irhs(jj)), nf)
			call matadd(z(irhs(jm1)), z(itau(jm1)), z(irhs(jm1)), nf)
		enddo
		call slvsm2(z(iu(1)), z(irhs(1)))
		nf=3
		do jj=2,j
			jm1=jj-1
			call rstrct(z(itemp(jm1)), z(iu(jj)), nf)
			call matsub(z(iu(jm1)), z(itemp(jm1)), z(itemp(jm1)),nf)
			nf=2*nf-1
			call interp(z(itau(jj)), z(itemp(jm1)),nf)
			call matadd(z(iu(jj)), z(itau(jj)), z(iu(jj)), nf)
			do jpost=1, NPOST
				call relax2(z(iu(jj)), z(irhs(jj)), nf)
			enddo
		enddo
		call lop(z(itemp(j)), z(iu(j)), nf)
		call matsub(z(itemp(j)), z(irhs(j)), z(itemp(j)), nf)
		res=anorm2(z(itemp(j)), nf)
		jcycle=jcycle+1
	enddo
!______________

enddo
call copy(u,z(iu(ngrid)),n)

end subroutine

end module
