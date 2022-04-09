module material_rates_m
	
	use datatypes,	only: interactions,particles
	use globvar, 	only: ntotal_loc,nvirt_loc,nhalo_loc
	use param, 		only: mass,dim,f,tenselem
	
	public:: art_visc,ext_force,int_force,con_density

contains
	
	!=================================================================================
	subroutine art_visc(ki,p_i,p_j,dwdx,ardvxdt)
	
		use param, only: alpha,beta,etq,hsml,c
	
		implicit none
		integer,intent(in):: ki
		type(particles),intent(in):: p_i,p_j
        real(f),intent(in):: dwdx(dim)
		real(f),intent(inout):: ardvxdt(dim,ntotal_loc+nvirt_loc+nhalo_loc)
		real(f):: dx(dim),piv(dim),muv,vr,rr,h,mrho
		
		dx(:) = p_i%x(:) - p_j%x(:)
        vr = DOT_PRODUCT(p_i%vx(:) - p_j%vx(:),dx(:))
		if (vr > 0_f) vr = 0._f
		
		!Artificial viscous force only if v_ij * r_ij < 0
		rr = DOT_PRODUCT(dx(:),dx(:))
		muv  = hsml*vr/(rr + hsml*hsml*etq*etq)
		mrho = 0.5_f*(p_i%rho + p_j%rho)
		piv  = (beta*muv - alpha*c)*muv/mrho*dwdx(:)
		
		ardvxdt(:,p_i%indloc) = ardvxdt(:,p_i%indloc) - mass*piv(:)
		ardvxdt(:,p_j%indloc) = ardvxdt(:,p_j%indloc) + mass*piv(:)
	
	end subroutine art_visc
	
	!=================================================================================
	subroutine ext_force(ki,p_i,p_j,dwdx,exdvxdt)
		
		use param, only: p1,p2,rr0,dd
	
		implicit none
		integer,intent(in):: ki
		type(particles),intent(in):: p_i,p_j
        real(f),intent(in):: dwdx(dim)
		real(f),intent(inout):: exdvxdt(dim,ntotal_loc+nvirt_loc+nhalo_loc)
		real(f):: dx(dim),rr,f
		
		dx(:) = p_i%x(:) - p_j%x(:)
		rr = SQRT(SUM(dx(:)*dx(:)))
		
		if (rr.lt.rr0) then
			f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
			exdvxdt(:,p_i%indloc) = exdvxdt(:,p_i%indloc) + dd*dx(:)*f
			exdvxdt(:,p_j%indloc) = exdvxdt(:,p_j%indloc) - dd*dx(:)*f
		endif
	
	end subroutine ext_force
	
	!=================================================================================
	subroutine int_force(ki,p_i,p_j,dwdx,indvxdt)
	
		implicit none
		integer,intent(in):: ki
		type(particles),intent(in):: p_i,p_j
        real(f),intent(in):: dwdx(dim)
		real(f),intent(inout):: indvxdt(dim,ntotal_loc+nvirt_loc+nhalo_loc)
		real(f):: h(dim)
		
		h = -(p_i%p/p_i%rho**2 + p_j%p/p_j%rho**2)*dwdx(:)
		indvxdt(:,p_i%indloc) = indvxdt(:,p_i%indloc) + mass*h(:)
		indvxdt(:,p_j%indloc) = indvxdt(:,p_j%indloc) - mass*h(:)
	
	end subroutine int_force
	
	!=================================================================================
	subroutine con_density(ki,p_i,p_j,dwdx,codrhodt)
	
		implicit none
		integer,intent(in):: ki
		type(particles),intent(in):: p_i,p_j
        real(f),intent(in):: dwdx(dim)
		real(f),intent(inout):: codrhodt(ntotal_loc+nvirt_loc+nhalo_loc)
		real(f):: dvx(dim),vcc
		
		dvx(:) = p_i%vx(:) - p_j%vx(:)
				
		vcc = DOT_PRODUCT(dvx(:),dwdx(:))
		
		codrhodt(p_i%indloc) = codrhodt(p_i%indloc) + mass*vcc
		codrhodt(p_j%indloc) = codrhodt(p_j%indloc) + mass*vcc
		
	end subroutine con_density
    
    !===============================================================================================================================
    subroutine strain_rate(ki,p_i,p_j,dwdx,dstraindt)
    
        implicit none
        integer,intent(in):: ki
        type(particles),intent(in):: p_i,p_j
        real(f),intent(in):: dwdx(dim)
        real(f),intent(inout):: dstraindt(tenselem,ntotal_loc+nvirt_loc+nhalo_loc)
        real(f):: dvx(dim),he(tenselem)
        
        dvx(:) = p_j%vx(:) - p_i%vx(:)
        
        he(1:3) = dvx(:)*dwdx(:)
        he(4) = 0.5_f*(dvx(1)*dwdx(2) + dvx(2)*dwdx(1))
        he(5) = 0.5_f*(dvx(1)*dwdx(3) + dvx(3)*dwdx(1))
        he(6) = 0.5_f*(dvx(2)*dwdx(3) + dvx(3)*dwdx(2))
        
        dstraindt(:,p_i%indloc) = dstraindt(:,p_i%indloc) + mass*he(:)/p_j%rho
        dstraindt(:,p_j%indloc) = dstraindt(:,p_j%indloc) + mass*he(:)/p_i%rho
        
    end subroutine strain_rate

end module material_rates_m
