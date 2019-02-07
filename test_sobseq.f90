program sobsequence

use random

implicit none

integer :: n
real, dimension(5) :: x
integer :: i,j

n=-1

call sobseq(n, x)

n=5

do i=1,15
call sobseq(n, x)
	do j=1,5
		print*, j, x(j)
	enddo
	print*, '___________________'
enddo

end program
