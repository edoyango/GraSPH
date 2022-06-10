module material_rates_m
    
    use datatypes,    only: interactions
    use globvar,     only: ntotal_loc,nvirt_loc,nhalo_loc,nghos_loc
    use param,         only: mass,dim,f
    
    public:: art_visc,ext_force,int_force,con_density

contains
    
    !=================================================================================
    subroutine art_visc(ki,pair,ardvxdt)
    
        use param, only: alpha,beta,etq,hsml,c
    
        implicit none
        integer,intent(in):: ki
        type(interactions),intent(in):: pair
        real(f),intent(inout):: ardvxdt(dim,ntotal_loc+nvirt_loc+nhalo_loc+nghos_loc)
        real(f):: dx(dim),piv(dim),muv,vr,rr,h,mrho
        
        dx(:) = pair%i%x(:) - pair%j%x(:)
        vr = DOT_PRODUCT(pair%i%vx(:) - pair%j%vx(:),dx(:))
        if (vr > 0_f) vr = 0._f
        
        !Artificial viscous force only if v_ij * r_ij < 0
        rr = DOT_PRODUCT(dx(:),dx(:))
        muv  = hsml*vr/(rr + hsml*hsml*etq*etq)
        mrho = 0.5_f*(pair%i%rho + pair%j%rho)
        piv  = (beta*muv - alpha*c)*muv/mrho*pair%dwdx(:)
        
        ardvxdt(:,pair%i%indloc) = ardvxdt(:,pair%i%indloc) - mass*piv(:)
        ardvxdt(:,pair%j%indloc) = ardvxdt(:,pair%j%indloc) + mass*piv(:)
    
    end subroutine art_visc
    
    !=================================================================================
    subroutine ext_force(ki,pair,exdvxdt)
        
        use param, only: p1,p2,rr0,dd
    
        implicit none
        integer,intent(in):: ki
        type(interactions),intent(in):: pair
        real(f),intent(inout):: exdvxdt(dim,ntotal_loc+nvirt_loc+nhalo_loc+nghos_loc)
        real(f):: dx(dim),rr,f
        
        dx(:) = pair%i%x(:) - pair%j%x(:)
        rr = SQRT(SUM(dx(:)*dx(:)))
        
        if (rr.lt.rr0) then
            f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
            exdvxdt(:,pair%i%indloc) = exdvxdt(:,pair%i%indloc) + dd*dx(:)*f
            exdvxdt(:,pair%j%indloc) = exdvxdt(:,pair%j%indloc) - dd*dx(:)*f
        endif
    
    end subroutine ext_force
    
    !=================================================================================
    subroutine int_force(ki,pair,indvxdt)
    
        implicit none
        integer,intent(in):: ki
        type(interactions),intent(in):: pair
        real(f),intent(inout):: indvxdt(dim,ntotal_loc+nvirt_loc+nhalo_loc+nghos_loc)
        real(f):: h(dim)
        
        h = -(pair%i%p/pair%i%rho**2 + pair%j%p/pair%j%rho**2)*pair%dwdx(:)
        indvxdt(:,pair%i%indloc) = indvxdt(:,pair%i%indloc) + mass*h(:)
        indvxdt(:,pair%j%indloc) = indvxdt(:,pair%j%indloc) - mass*h(:)
    
    end subroutine int_force
    
    !=================================================================================
    subroutine con_density(ki,pair,codrhodt)
    
        implicit none
        integer,intent(in):: ki
        type(interactions),intent(in):: pair
        real(f),intent(inout):: codrhodt(ntotal_loc+nvirt_loc+nhalo_loc+nghos_loc)
        real(f):: dvx(dim),vcc
        
        dvx(:) = pair%i%vx(:) - pair%j%vx(:)
                
        vcc = DOT_PRODUCT(dvx(:),pair%dwdx(:))
        
        codrhodt(pair%i%indloc) = codrhodt(pair%i%indloc) + mass*vcc
        codrhodt(pair%j%indloc) = codrhodt(pair%j%indloc) + mass*vcc
        
    end subroutine con_density

end module material_rates_m
