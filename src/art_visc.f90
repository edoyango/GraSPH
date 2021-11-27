!======================================================================================
subroutine art_visc(ardvxdt)
!======================================================================================

use param
use globvar
implicit none
integer:: i,j,k,d
real(8):: dx(dim),piv(dim),muv,vr,rr,h,mrho
real(8),intent(out):: ardvxdt(dim,ntotal+nvirt)
real(8),parameter:: alpha = 0.01d0, beta = 0.d0, etq = 0.1d0
           
ardvxdt(:,1:ntotal) = 0d0
     
!Calculate SPH sum for artificial viscosity

do k=1,niac

  if (pairs(k)%i%itype.gt.0 .and. pairs(k)%j%itype.gt.0) then
  
	dx(:) = pairs(k)%i%x(:) - pairs(k)%j%x(:)
	vr = MIN(0d0,DOT_PRODUCT(pairs(k)%i%vx(:) - pairs(k)%j%vx(:),dx(:)))
	
    !Artificial viscous force only if v_ij * r_ij < 0
    rr = DOT_PRODUCT(dx(:),dx(:))
    muv  = hsml*vr/(rr + hsml*hsml*etq*etq)
    mrho = 0.5d0*(pairs(k)%i%rho + pairs(k)%j%rho)
    piv  = (beta*muv - alpha*c)*muv/mrho*pairs(k)%dwdx(:)
	
	ardvxdt(:,pairs(k)%i%ind) = ardvxdt(:,pairs(k)%i%ind) - mass*piv(:)
	ardvxdt(:,pairs(k)%j%ind) = ardvxdt(:,pairs(k)%j%ind) + mass*piv(:)
	  
  endif
  
enddo

end