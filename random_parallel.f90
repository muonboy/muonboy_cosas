program random_paralelo

implicit none

real :: x
integer :: j

!Llamamos a la libreria de paralelizacion. Se paraleliza en nucleos logicos

!$omp parallel do
do j = 1, 10
	print*, 'Experimento',j,': ', parallel_random()
	
enddo
!$omp end parallel do


!____________________________________________________________________________
contains

function parallel_random() result(y)

real :: y, u

y = 0.0

	call random_number(u)			!un aleatorio uniforme de media cero y varianza 1
	y = y + u

end function

end program
