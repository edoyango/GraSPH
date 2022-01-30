!=======================================================================
subroutine ext_force(ki,exdvxdt)
!=======================================================================

	use globvar
	use param
	use physical_constants
	
	implicit none
	integer,intent(in):: ki
	integer:: i,j,k,d
	real(8),intent(out):: exdvxdt(dim,ntotal+nvirt)
	real(8):: dx(dim),rr,f
	real(8),parameter:: rr0 = dxo,dd = 5d0*g*25d0
	integer,parameter:: p1=4,p2=2

	!Initialisation
	do i = 1,ntotal
		exdvxdt(1:dim-1,i) = 0d0
		exdvxdt(dim,i) = -g
	end do
	
	!Boundary particle force and penalty anti-penetration force      
	do k=1,niac
	
		if (pairs(k)%i%itype*pairs(k)%j%itype.lt.0) then  
		
			dx(:) = pairs(k)%i%x(:) - pairs(k)%j%x(:)
			rr = SQRT(SUM(dx(:)*dx(:)))
		
			if (rr.lt.rr0) then
				f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
				exdvxdt(:,pairs(k)%i%ind) = exdvxdt(:,pairs(k)%i%ind) + dd*dx(:)*f
				exdvxdt(:,pairs(k)%j%ind) = exdvxdt(:,pairs(k)%j%ind) - dd*dx(:)*f
			endif
			
		endif
	enddo

end subroutine ext_force