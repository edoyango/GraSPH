!==============================================================================================
subroutine int_force(ki,indvxdt)
!==============================================================================================

	use globvar
	use param
	
	implicit none
	integer,intent(in):: ki
	integer:: i,j,k,d
	real(8):: h(dim)
	real(8),parameter:: rh0=irho, Hmax=25d0
	integer,parameter:: gamma = 7
	real(8),intent(out):: indvxdt(dim,ntotal+nvirt)
	
	!Initialized acceleration
	indvxdt(:,1:ntotal) = 0d0
	
	!Compute pressure using stiff equation (Monaghan, 1994)
	parts(1:ntotal)%p = rh0*c**2*((parts(1:ntotal)%rho/rh0)**gamma-1d0)/gamma
	
	!Compute acceleation due to pressure
	do k=1,niac
	
		if (pairs(k)%i%itype.gt.0 .and. pairs(k)%j%itype.gt.0) then
			h = -(pairs(k)%i%p/pairs(k)%i%rho**2 + pairs(k)%j%p/pairs(k)%j%rho**2)*pairs(k)%dwdx(:)
			indvxdt(:,pairs(k)%i%ind) = indvxdt(:,pairs(k)%i%ind) + mass*h(:)
			indvxdt(:,pairs(k)%j%ind) = indvxdt(:,pairs(k)%j%ind) - mass*h(:)
		endif
	
	enddo

end subroutine int_force
