!=============================================================
subroutine virt_part( )
!=============================================================
use param
implicit none
integer:: i,k,d,n
real(8):: xl,xmin,xmax,ymin,ymax,dx,dy,v_inf
character(len=3)::number

!---Given particles on solid boundary
dx = dxo ! seperation between particles
dy = dxo ! seperation between particles
n =  ntotal! initial number of virtual particles
v_inf = 0.d0

!---Defined parameters for creating type two of virtual boundary
xmin = 0.0d0
ymin = 0.0d0
xmax = 075.d0
ymax = 040.d0

!---Virtual particle on the lower boundary
do i = 1, op
  n = n + 1
  parts(n)%ind = n
  parts(n)%x(1) = xmin + i*dx - 0.5d0*dx
  parts(n)%x(2) = ymin - 0.5d0*dy
  parts(n)%vx(:) = 0d0
enddo

!---Virtual particle on the upper boundary
do i = 1, op
  n = n + 1
  parts(n)%ind = n
  parts(n)%x(1) = xmin + i*dx - 0.5d0*dx
  parts(n)%x(2) = ymax - 1.5d0*dy
  parts(n)%vx(:) = 0d0
enddo

!---Virtual particle on the left boundary
do i = 1, pp
  n = n + 1
  parts(n)%ind = n
  parts(n)%x(1) = xmin - 0.5d0*dx
  parts(n)%x(2) = ymin + (i-1)*dy - 0.5d0*dy
  parts(n)%vx(:) = 0d0
enddo

!---Virtual particle on the right boundary
do i = 1, pp
  n = n + 1
  parts(n)%ind = n
  parts(n)%x(1) = xmax + 0.5d0*dx
  parts(n)%x(2) = ymin + (i-1)*dy - 0.5d0*dy
  parts(n)%vx(:) = 0d0
enddo

parts(ntotal+1:ntotal+nvirt)%rho = irho
parts(ntotal+1:ntotal+nvirt)%p = 0d0
parts(ntotal+1:ntotal+nvirt)%itype = -1

open(1,file="output_data/xv_vp.dat")
open(2,file="output_data/state_vp.dat")
	  	
do i = ntotal+1, n         
  write(1,1001) i, parts(i)%x(:), parts(i)%vx(:)           
  write(2,1002) i, parts(i)%itype, hsml, mass, parts(i)%rho, parts(i)%p
enddo

1001 format(1x, I6, 6(2x, e14.7))
1002 format(1x, I6, I6,  8(2x, e14.7)) 

close(1)
close(2) 

!Print statistic
print *,' >> Statistics: Virtual boundary particles:'
print *,'       Number of virtual particles:',nvirt

end