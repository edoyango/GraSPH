!==========================================================
subroutine input( )
!==========================================================

use param
use globvar
implicit none
integer:: i,d,im

open(10,file="C:\Users\edwar\Documents\outputdata/ini_xv.dat")
open(20,file="C:\Users\edwar\Documents\outputdata/ini_state.dat")
open(30,file="C:\Users\edwar\Documents\outputdata/ini_stress.dat") 

call model_fluid( )
do i = 1, ntotal 
  write(10,1001) i, parts(i)%x(:), parts(i)%vx(:)
  write(20,1002) i, parts(i)%itype, hsml, mass, parts(i)%rho, parts(i)%p
enddo

1001 format(1x, I5, 6(2x, e14.7)) 
1002 format(1x, I5, 2x, I2, 8(2x, e14.7)) 
1003 format(1x, I5, 2x, I2, 2x, e14.7) 

write(*,*)'  **************************************************'
write(*,*)'      Initial particle configuration generated   '       
write(*,*)'      Total number of particles   ', ntotal    	
write(*,*)'  **************************************************' 

close(1)
close(2) 
close(3) 

end              
       
!===============================================================
subroutine model_fluid( )
!===============================================================
use param
use globvar
implicit none
integer:: i,j,d,n
real(8):: xl,yl,xi,yi,dx,dy

!---Giving mass and smoothing length as well as other data.
yl = 25.d0  ! maximum high geometry
xl = 25.d0  ! maximum extend of fluid
dx = dxo ! speration of particle in x	
dy = dxo ! seperation of particle in y
n = 0

!---Defined parameter for initial density setup	
do i = 1, mp
  do j = 1, np
	xi = i*dx - dx/2.d0
	yi = j*dy - dy/2.d0
	n = n + 1
	parts(n)%ind = n
	parts(n)%x(1) = xi
	parts(n)%x(2) = yi
	parts(n)%vx(:) = 0d0
	parts(n)%itype = 2
	parts(n)%rho = irho
	parts(n)%p = 0d0
  enddo
enddo

end