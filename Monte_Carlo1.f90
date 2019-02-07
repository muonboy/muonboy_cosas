program mc1

implicit none

integer, parameter :: Nmc = 100000, Npart = 1E5
real, dimension(Nmc) :: mu
real :: mean, st_dev, antes, despues
integer :: j

!Llamamos a la libreria de paralelizacion. Se paraleliza en nucleos logicos

call cpu_time(antes)
!$omp parallel do
do j = 1, Nmc
	!print*, 'Experimento',j,'-esimo'
	mu(j) = monte_carlo(Npart) 
enddo
!$omp end parallel do
call cpu_time(despues)

print*, 'tiempo= ', despues-antes, 'segundos'

mean = sum(mu) / dble(Nmc)
st_dev = sqrt(sum((mu - mean)**2))*(1/sqrt(dble(Nmc)))

print*, mean
print*, st_dev

!____________________________________________________________________________
contains

function monte_carlo(n) result(y)

integer, intent(in) :: n
integer :: i
real :: y, u

y = 0.0
do i = 1, n
	call random_number(u)			!un aleatorio uniforme de media cero y varianza 1
	y = y + u
enddo

y = y / dble(n)

end function

end program
