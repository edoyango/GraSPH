!=============================================================================================================
subroutine single_step(ki,dvxdti,drhoi) 
!=============================================================================================================

use param
use globvar
implicit none
integer,intent(in):: ki
real(8),intent(out):: dvxdti(dim,ntotal),drhoi(dim,ntotal+nvirt)
integer:: i,j,k,d
real(8):: t1,t2,indvxdt(dim,ntotal+nvirt),ardvxdt(dim,ntotal+nvirt),exdvxdt(dim,ntotal+nvirt)

!Density approximation or change rate
call CPU_TIME(t1)
call con_density(ki,drhoi)
call CPU_TIME(t2)
testtime(2) = testtime(2) + t2 - t1
!Internal force due to pressure
call CPU_TIME(t1)
call int_force(indvxdt)
call CPU_TIME(t2)
  testtime(3) = testtime(3) + t2 - t1
!Artificial viscosity:
call CPU_TIME(t1)
call art_visc(ardvxdt)
call CPU_TIME(t2)
  testtime(4) = testtime(4) + t2 - t1
!External forces:
call CPU_TIME(t1)
call ext_force(exdvxdt)
call CPU_TIME(t2)
  testtime(5) = testtime(5) + t2 - t1

!Convert velocity, force, and energy to f and dfdt
dvxdti(1:dim,1:ntotal) = indvxdt(1:dim,1:ntotal) + exdvxdt(1:dim,1:ntotal) + ardvxdt(1:dim,1:ntotal)

end
