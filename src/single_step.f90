!=============================================================================================================
subroutine single_step(ki,dvxdti,drhoi) 
!=============================================================================================================

	use globvar
	use param
	
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
	
	!Density approximation or change rate
	call con_density(ki,cdrhodt)
	
	!Internal force due to pressure
	call int_force(ki,indvxdt)
	
	!Artificial viscosity:
	call art_visc(ki,ardvxdt)
	
	!External forces:
	call ext_force(ki,exdvxdt)
	
	!Convert velocity, force, and energy to f and dfdt
	dvxdti(1:dim,1:ntotal) = indvxdt(1:dim,1:ntotal) + exdvxdt(1:dim,1:ntotal) + ardvxdt(1:dim,1:ntotal)
	drhoi(1:ntotal) = cdrhodt(1:ntotal)
	
	deallocate( indvxdt,ardvxdt,exdvxdt,cdrhodt )

end
