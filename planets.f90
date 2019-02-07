program planets_dynamics

implicit none

integer, parameter :: ndim=2, npart=3, steps = 100000
real*8 :: pot, kin, xm, vm, w
real*8, parameter :: dt = 0.0001d0, Ms=1.d0, Mp=0.1d0, Mm=0.01d0
real*8, dimension(steps) :: E, Ki, U, t
integer :: i, j, k
real*8, dimension(npart) :: masa
real*8, dimension (npart,ndim) :: posicion, posicion_CM, f, p
real*8, dimension(2) :: r_moon

masa(1) = Ms
masa(2) = Mp
masa(3) = Mm

!SOl
posicion(1,1) = 0.d0
posicion(1,2) = 0.d0

p(1,1) = 0.d0
p(1,2) = 0.d0

!PLANETA
posicion(2,1) = 1.d0
posicion(2,2) = 0.d0

p(2,1) = 0.d0
p(2,2) = 1.d0*Mp

!LUNA
xm = 1.075d0
vm = 3.d0

posicion(3,1) = xm
posicion(3,2) = 0.d0

p(3,1) = 0.d0
p(3,2) = vm*Mm

call fuerzas(npart,ndim,posicion,p,masa, f, pot,kin)

posicion_CM = 0.d0

call centro(npart,ndim,masa,posicion, posicion_CM)

open(file='posiciones.dat', unit=11, status='unknown')
open(file='posiciones_CM.dat', unit=10, status='unknown')
open(file='energy.dat', unit=7, status='unknown')
open(file='pos_moon.dat', unit=17, status='unknown')


do i=1,steps
	t(i)=(i-1)*dt

	call pos_moon(posicion_CM, r_moon, w)

	write(11,*) t(i), posicion(1,1), posicion(1,2), posicion(2,1), posicion(2,2), posicion(3,1), posicion(3,2)

	write(10,*) t(i), posicion_CM(1,1), posicion_CM(1,2), posicion_CM(2,1), posicion_CM(2,2), posicion_CM(3,1), posicion_CM(3,2)

	write(17,*) t(i), r_moon(1), r_moon(2), w

!_________________________Libre__________________________________

	call fuerzas(npart,ndim,posicion,p,masa, f,pot,kin)

	Ki(i) = kin
	U(i) = pot
	E(i) = kin + pot

	call leapfrog(npart, ndim, posicion, p, f, masa, dt)

!__________________________CM___________________________________

	call centro(npart,ndim,masa,posicion, posicion_CM)

	write(7,*) t(i), E(i)-E(1) 
enddo

 close(17)
 close(7)
 close(10)
 close(11)


contains

!__________________________________________________________________________________________

subroutine distancia(ndim,r1,r2, despl,dist) !desp=desplazamiento
implicit none
 integer::ndim
 real*8, dimension(ndim) :: r1, r2, despl
 real*8::dist
 integer::i
 do i=1, ndim
 	despl(i) = r2(i) - r1(i)
 enddo

dist = 0.d0

do i=1,ndim
	dist = dist + despl(i)*despl(i)
enddo
dist=dsqrt(dist)

end subroutine distancia

!________________________________________________________________________________________

subroutine fuerzas(npart,ndim,posicion,p,masa,  f,pot,kin) 		!npart=numero de particulas que rodean, vel=velocidades, f=fuerza, 

implicit none
integer :: i,j,k, ndim, npart
real*8, dimension(npart, ndim) :: f, posicion, p
real*8, dimension(ndim) :: rji
real*8, dimension(npart) :: masa
real*8 :: pot,kin, d, r, G

G = 4.d0
pot=0.d0
kin=0.d0

!!!!!!!!!!!!!comienzo el calculo de energia potencial y de fuerzas
f = 0.d0
do i=1, npart
   	do j=1, npart
     		if(i.ne.j) then
      			call distancia(ndim, posicion(i,:), posicion(j,:), rji, d)

			pot = pot - G*masa(j)/d
			
			r = dsqrt(rji(1)**2.d0 + rji(2)**2.d0)			
			
    			f(i,1) = f(i,1) + G*(masa(i)*masa(j)/(d*d*d))*rji(1)
		        f(i,2) = f(i,2) + G*(masa(i)*masa(j)/(d*d*d))*rji(2)
     		endif	
   	enddo

    	do k=1, ndim
        	 kin = kin + 0.5d0*p(k,i)*p(k,i)/masa(i)
    	enddo

enddo

end subroutine fuerzas

!______________________________________________________________________________________________________

subroutine leapfrog(npart, ndim, posicion, p, f, masa, dt)
implicit none

!input-output:: posicion, velocidad, aceleracion
!para paralelizar filas y columnas, se lo realiza desde fuera; es decir, primero columnas y luego filas

integer::npart, ndim,i,j
real*8, dimension(npart, ndim) :: f, posicion, p
real*8 :: rmasa,dt
real*8, dimension(npart) :: masa

do j=1, npart
	rmasa = 1.d0/masa(j)
    do i=1, ndim
	p(j,i) = p(j,i) + 0.5d0*f(j,i)*dt
        posicion(j,i) = posicion(j,i) + rmasa*p(j,i)*dt
    enddo
enddo

end subroutine leapfrog

!_____________________________________________________________________________

subroutine centro(npart,ndim,masa,posicion_in, posicion_out)
implicit none
integer :: ndim, npart
real*8, dimension(npart, ndim) :: posicion_in, posicion_out
real*8, dimension(ndim) :: rcm
real*8, dimension(npart) :: masa
real*8 :: rx,ry,mass
integer :: j,k

mass = 0.d0

do j=1,npart
	mass = mass + masa(j)
enddo

rx = 0.d0
ry = 0.d0

do j=1,npart
	rx = rx + posicion_in(j,1)*masa(j)
	ry = ry + posicion_in(j,2)*masa(j)
enddo
	rcm(1) = rx/mass
	rcm(2) = ry/mass

do k=1,npart
	do j=1,ndim
		posicion_out(k,j) = posicion_in(k,j) - rcm(j) 
  	enddo
enddo

end subroutine

!_____________________________________________________

subroutine pos_moon(r_in, r_moon, u)
implicit none
real*8, dimension(3, 2) :: r_in
real*8, dimension(2) :: r_moon
real*8 :: u
integer :: i,j

u = 0.d0
r_moon = 0.d0

r_moon(1) = r_in(3,1) - r_in(2,1)
r_moon(2) = r_in(3,2) - r_in(2,2)
u = dsqrt(r_moon(1)**2 + r_moon(2)**2)

end subroutine

end program 
