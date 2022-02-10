module single_step_m
contains

	!==============================================================================================================================
	subroutine single_step(ki,dvxdti,drhoi) 
	! Container subroutine for all the rate-of-change calculations. Rate-of-changes are calculated seperately and then summed as
	! required
	
		use globvar
		use param
		
		use material_rates_m
		
		implicit none
		integer,intent(in):: ki
		real(8),intent(out):: dvxdti(dim,ntotal),drhoi(ntotal)
		integer:: i,j,k,d
		real(8):: t1,t2
		real(8),allocatable:: indvxdt(:,:),ardvxdt(:,:),exdvxdt(:,:),cdrhodt(:)
		
		allocate( indvxdt(dim,ntotal+nvirt),&
				ardvxdt(dim,ntotal+nvirt),&
				exdvxdt(dim,ntotal+nvirt),&
				cdrhodt(ntotal+nvirt) )
				
		cdrhodt(1:ntotal) = 0d0
		indvxdt(:,1:ntotal) = 0d0
		ardvxdt(:,1:ntotal) = 0d0
		exdvxdt(1,1:ntotal) = 0d0
		exdvdxt(2,1:ntotal) = -g
		
		do k = 1,niac
			
			if (pairs(k)%i%itype > 0 .and. pairs(k)%j%itype > 0) then
			
				!Density approximation or change rate
				call con_density(ki,pairs(k),cdrhodt)
				
				!Internal force due to pressure
				call int_force(ki,pairs(k),indvxdt)
				
				!Artificial viscosity:
				call art_visc(ki,pairs(k),ardvxdt)
				
			elseif (pairs(k)%i%itype> 0 .or. pairs(k)%j%itype > 0) then
			
				!External forces:
				call ext_force(ki,pairs(k),exdvxdt)
			
			end if
		
		end do
		
		!Convert velocity, force, and energy to f and dfdt
		dvxdti(1:dim,1:ntotal) = indvxdt(1:dim,1:ntotal) + exdvxdt(1:dim,1:ntotal) + ardvxdt(1:dim,1:ntotal)
		drhoi(1:ntotal) = cdrhodt(1:ntotal)
		
		deallocate( indvxdt,ardvxdt,exdvxdt,cdrhodt )
	
	end subroutine single_step

end module single_step_m