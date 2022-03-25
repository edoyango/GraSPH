module material_rates_m
	
    use datatypes, only: particles,interactions
	use globvar, only: ntotal,nvirt,nghos
	use param, only: mass,dim,f,tenselem

contains
	
	!=================================================================================
	subroutine art_visc(ki,p_i,p_j,dwdx,ardvxdt)
	
		use param, only: alpha,beta,etq,hsml,c
	
		implicit none
		integer,intent(in):: ki
        type(particles),intent(in):: p_i,p_j
        real(f),intent(in):: dwdx(dim)
		real(f),intent(inout):: ardvxdt(dim,ntotal+nvirt+nghos)
		real(f):: dx(dim),piv(dim),muv,vr,rr,h,mrho
		
		dx(:) = p_i%x(:) - p_j%x(:)
		vr = DOT_PRODUCT(p_i%vx(:) - p_j%vx(:),dx(:))
        if (vr > 0_f) vr = 0_f
		
		!Artificial viscous force only if v_ij * r_ij < 0
		rr = DOT_PRODUCT(dx(:),dx(:))
		muv  = hsml*vr/(rr + hsml*hsml*etq*etq)
		mrho = 0.5_f*(p_i%rho + p_j%rho)
		piv  = (beta*muv - alpha*c)*muv/mrho*dwdx(:)
		
		ardvxdt(:,p_i%ind) = ardvxdt(:,p_i%ind) - mass*piv(:)
		ardvxdt(:,p_j%ind) = ardvxdt(:,p_j%ind) + mass*piv(:)
	
	end subroutine art_visc
	
	!=================================================================================
	subroutine ext_force(ki,pair,exdvxdt)
		
		use param, only: p1,p2,rr0,dd
	
		implicit none
		integer,intent(in):: ki
		type(interactions),intent(in):: pair
		real(f),intent(inout):: exdvxdt(dim,ntotal+nvirt+nghos)
		real(f):: dx(dim),rr,f
		
		dx(:) = pair%i%x(:) - pair%j%x(:)
		rr = SQRT(SUM(dx(:)*dx(:)))
		
		if (rr.lt.rr0) then
			f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
			exdvxdt(:,pair%i%ind) = exdvxdt(:,pair%i%ind) + dd*dx(:)*f
			exdvxdt(:,pair%j%ind) = exdvxdt(:,pair%j%ind) - dd*dx(:)*f
		endif
	
	end subroutine ext_force
	
	!=================================================================================
	subroutine int_force(ki,p_i,p_j,dwdx,indvxdt)
	
		implicit none
		integer,intent(in):: ki
		type(particles),intent(in):: p_i,p_j
        real(f),intent(in):: dwdx(dim)
		real(f),intent(inout):: indvxdt(dim,ntotal+nvirt+nghos)
		real(f):: h(dim)
		
		h = -(p_i%p/p_i%rho**2 + p_j%p/p_j%rho**2)*dwdx(:)
		indvxdt(:,p_i%ind) = indvxdt(:,p_i%ind) + mass*h(:)
		indvxdt(:,p_j%ind) = indvxdt(:,p_j%ind) - mass*h(:)
	
	end subroutine int_force
	
	!=================================================================================
	subroutine con_density(ki,p_i,p_j,dwdx,codrhodt)
	
		implicit none
		integer,intent(in):: ki
		type(particles),intent(in):: p_i,p_j
        real(f),intent(in):: dwdx(dim)
		real(f),intent(inout):: codrhodt(ntotal+nvirt+nghos)
		real(f):: dvx(dim),vcc
		
		dvx(:) = p_i%vx(:) - p_j%vx(:)
				
		vcc = DOT_PRODUCT(dvx(:),dwdx(:))
		
		codrhodt(p_i%ind) = codrhodt(p_i%ind) + mass*vcc
		codrhodt(p_j%ind) = codrhodt(p_j%ind) + mass*vcc
		
	end subroutine con_density
    
    !===============================================================================================================================
    subroutine strain_rate(ki,p_i,p_j,dwdx,dstraindt)
    
        implicit none
        integer,intent(in):: ki
        type(particles),intent(in):: p_i,p_j
        real(f),intent(in):: dwdx(dim)
        real(f),intent(inout):: dstraindt(tenselem,ntotal+nvirt+nghos)
        real(f):: dvx(dim),he(tenselem),hr(tenselem-dim)
        
        dvx(:) = p_i%vx(:) - p_j%vx(:)
        
        he(4) = 0.5_f*(dvx(1)*dwdx(2) + dvx(2)*dwdx(1))
        he(5) = 0.5_f*(dvx(1)*dwdx(3) + dvx(3)*dwdx(1))
        he(6) = 0.5_f*(dvx(2)*dwdx(3) + dvx(3)*dwdx(2))
        
        dstraindt(:,p_i%ind) = dstraindt(:,p_i%ind) + mass*he(:)/p_j%rho
        dstraindt(:,p_j%ind) = dstraindt(:,p_j%ind) + mass*he(:)/p_i%rho
        
    end subroutine strain_rate
        

end module material_rates_m
