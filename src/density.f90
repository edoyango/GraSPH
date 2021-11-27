!=================================================================================
subroutine con_density(ki,drhoi)
!=================================================================================

use param
use globvar
implicit none
integer,intent(in):: ki
real(8),intent(out):: drhoi(ntotal+nvirt)
integer:: i,j,k,d
real(8):: dvx(dim),vcc

do k=1,niac      

  if (pairs(k)%i%itype.gt.0 .and. pairs(k)%j%itype.gt.0) then

    dvx(:) = pairs(k)%i%vx(:) - pairs(k)%j%vx(:) 
    
    vcc = DOT_PRODUCT(dvx(:),pairs(k)%dwdx(:))
    
    drhoi(pairs(k)%i%ind) = drhoi(pairs(k)%i%ind) + mass*vcc
    drhoi(pairs(k)%j%ind) = drhoi(pairs(k)%j%ind) + mass*vcc       
  
  endif

enddo
	 
end