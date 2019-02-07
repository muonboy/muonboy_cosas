program planetas

implicit none

integer, parameter :: ndim=2, npart=3, steps = 30000
real*8 :: pot, kin, xm, vm, pot_CM, kin_CM, w
real*8, parameter :: dt = 0.0001d0, Ms=1.d0, Mp=0.1d0, Mm=0.01d0
real*8, dimension(steps) :: E, Ki, U, t, E_CM, Ki_CM, U_CM
integer :: i, j, k
real*8, dimension(npart) :: masa
real*8, dimension (npart,ndim) :: posicion, posicion_CM, f,f_CM, p, p_CM
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
xm = 1.2d0
vm = 3.0d0

posicion(3,1) = xm
posicion(3,2) = 0.d0

p(3,1) = 0.d0
p(3,2) = vm*Mm

call fuerzas(npart,ndim,posicion,p,masa, f, pot,kin)

!print*, posicion(1,1), posicion(1,2)
!print*, posicion(2,1), posicion(2,2)
!print*, posicion(3,1), posicion(3,2)

posicion_CM = 0.d0

call centro(npart,ndim,masa,posicion, posicion_CM)

!print*, posicion_CM(1,1), posicion_CM(1,2)
!print*, posicion_CM(2,1), posicion_CM(2,2)
!print*, posicion_CM(3,1), posicion_CM(3,2)

call centro_vel(npart,ndim,masa,p, p_CM)

call fuerzas(npart,ndim,posicion_CM,p_CM,masa, f_CM,pot_CM,kin_CM)

open(file='posiciones.dat', unit=11, status='unknown')
open(file='posiciones_CM.dat', unit=10, status='unknown')
!open(file='p.dat', unit=9, status='unknown')
!open(file='p_CM.dat', unit=8, status='unknown')
open(file='energy.dat', unit=7, status='unknown')
open(file='energy_CM.dat', unit=16, status='unknown')
open(file='pos_moon.dat', unit=17, status='unknown')


do i=1,steps
	t(i)=(i-1)*dt

	call pos_moon(posicion_CM, r_moon, w)

	write(11,*) t(i), posicion(1,1), posicion(1,2), posicion(2,1), posicion(2,2), posicion(3,1), posicion(3,2)
!	write(9,*) t(i), p(1,1), p(1,2), p(2,1), p(2,2), p(3,1), p(3,2) 

	write(10,*) t(i), posicion_CM(1,1), posicion_CM(1,2), posicion_CM(2,1), posicion_CM(2,2), posicion_CM(3,1), posicion_CM(3,2)
!	write(8,*) t(i), p_CM(1,1), p_CM(1,2), p_CM(2,1), p_CM(2,2), p_CM(3,1), p_CM(3,2)

	write(17,*) t(i), r_moon(1), r_moon(2), w

!_________________________Libre__________________________________

	call fuerzas(npart,ndim,posicion,p,masa, f,pot,kin)

	Ki(i) = kin
	U(i) = pot
	E(i) = kin + pot

	call leapfrog(npart, ndim, posicion, p, f, masa, dt)

!__________________________CM___________________________________

	call centro(npart,ndim,masa,posicion, posicion_CM)
	call centro_vel(npart,ndim,masa,p, p_CM)

	call fuerzas(npart,ndim,posicion_CM,p_CM,masa, f_CM,pot_CM,kin_CM)

	Ki_CM(i) = kin_CM
	U_CM(i) = pot_CM
	E_CM(i) = kin_CM + pot_CM

	call leapfrog(npart, ndim, posicion_CM, p_CM, f_CM, masa, dt)

	write(7,*) t(i), E(i)-E(1) 
	write(16,*) t(i), E_CM(i)-E_CM(1) 
enddo

 close(17)
 close(16)
 close(7)
! close(8)
! close(9)
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
	!print*, j, ' ______', rx, ry

enddo
	rcm(1) = rx/mass
	rcm(2) = ry/mass

!print*, ' ______________', rcm(1), rcm(2)

do k=1,npart
	do j=1,ndim
		posicion_out(k,j) = posicion_in(k,j) - rcm(j) 
  	enddo
enddo

end subroutine

!______________________________________________________


subroutine centro_vel(npart,ndim,masa,p_in, p_out)
implicit none
integer :: ndim, npart
real*8, dimension(npart, ndim) :: p_in, p_out
real*8, dimension(ndim) :: vcm
real*8, dimension(npart) :: masa
real*8 :: px,py,mass
integer :: i,j

mass=0.d0
px=0.d0
py=0.d0

do j=1,npart
	mass = mass + masa(j)
	px = px + p_in(j,1)
	py = py + p_in(j,2)
enddo

vcm(1) = px/mass
vcm(2) = py/mass

do i=1, npart
	do j=1, ndim
		p_out(i,j) = p_in(i,j) - vcm(j)*masa(i)
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
