program walk
use random

implicit none

integer :: N
real*8 :: D, Drms, x0, xf

xf=0.d0
D=0.d0
Drms=0.d0

print*, 'x0: '
read*, x0
print*, 'N: '
read*, N

call random_walk(N, x0, xf, D, Drms)

print*, 'N, x0, xf, D, Drms', N, x0, xf, D, Drms

end program
