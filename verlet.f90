!molecular dinamic--- algoritmo de verlet
module verlet_

contains
subroutine distancia(ndim,r1,r2,    despl,dist ) !desp=desplazamiento
implicit none
 integer::ndim
 real*8, dimension(ndim)::r1, r2, despl
 !real*8::r1(ndim),r2(ndim)
 !real*8::despl(ndim), dist
 real*8::dist
 integer::i
 do i=1, ndim
    despl(i)=r1(i)- r2(i)
 enddo

dist=0.d0
do i=1,ndim
  dist= dist + despl(i)*despl(i)
enddo
dist=dsqrt(dist)
end subroutine distancia

subroutine fuerzas(npart, ndim, posicion, p, masa,  f, pot, kin  ) 		!npart=numero de particulas que rodean, vel=velocidades, f=fuerza, pot=potencial, kin=E.cinetica
implicit none
integer::i,j,k, ndim, npart
real*8, dimension(ndim, npart) :: f, posicion, p
real*8, dimension(ndim)::rij
!real*8::f(ndim, npart)
!real*8::posicion(ndim, npart)
!real*8::vel(ndim, npart)
real*8::pi2, pot, kin,masa,d

 pi2=dacos(-1.d0)/2
pot=0.d0
kin=0.d0

!!!!!!!!!!!!!comienzo el calculo de energia potencial y de fuerzas

do i=1, npart
   do k=1, ndim
    f(k,i)=0.d0  								!esta inicializacion sirve para posteriormente paralelizar
   enddo

   do j=1, npart
     if(i .ne.j) then
      call distancia(ndim, posicion(:,i), posicion(:,j), rij, d)

!!!!!!usar potencial central tipo ''pozo de potencial''

	pot=pot + 0.5d0*(dsin(min(d,pi2)))**2
  !! pot=pot +1.d0/d
     	do k=1, ndim
    		f(k,i)=f(k,i)-rij(k)*2.d0*dsin(min(d,pi2))*dcos(min(d,pi2))/d
     	enddo
     endif
   enddo

    do k=1, ndim
         kin = kin + 0.5d0*p(k,i)*p(k,i)
    enddo
enddo


print*, 'fuerzas_______________'

end subroutine fuerzas

subroutine verlet_velocity(nparti, ndim, posicion, vel, f, acel, masa, dt)
implicit none

!input-output:: posicion, velocidad, aceleracion
!para paralelizar filas y columnas, se lo realiza desde fuera; es decir, primero columnas y luego filas

integer::nparti, ndim,i,j
real*8, dimension(ndim, nparti)::f, posicion, vel,acel
real*8::masa,rmasa,dt
rmasa=1.d0/masa
do j=1, nparti
    do i=1, ndim
        posicion(i,j)=posicion(i,j) + vel(i,j)*dt + 0.5d0*acel(i,j)*dt*dt
	vel(i,j)=vel(i,j) + 0.5d0*dt*(f(i,j)*rmasa + acel(i,j))
	acel(i,j)=f(i,j)*rmasa
    enddo
enddo


end subroutine verlet_velocity

subroutine leapfrog(npart, ndim, posicion, p, f, acel, masa, dt)
implicit none

!input-output:: posicion, velocidad, aceleracion
!para paralelizar filas y columnas, se lo realiza desde fuera; es decir, primero columnas y luego filas

integer::npart, ndim,i,j
real*8, dimension(ndim, npart) :: f, posicion, p, acel, vel
real*8 :: masa,rmasa,dt

rmasa=1.d0/masa
do j=1, npart
    do i=1, ndim
	p(i,j) = p(i,j) + 0.5d0*f(i,j)*dt
        posicion(i,j) = posicion(i,j) + p(i,j)*dt
    enddo
enddo

print*, 'leapfrog_____________________'

end subroutine leapfrog

end module

