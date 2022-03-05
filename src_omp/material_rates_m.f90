module material_rates_m
	
	use globvar, only: interactions,ntotal,nvirt
	use param, only: mass,dim,f

contains
	
	!=================================================================================
	subroutine art_visc(ki,pair,ardvxdt)
	
		use param, only: alpha,beta,etq,hsml,c
	
		implicit none
		integer,intent(in):: ki
		type(interactions),intent(in):: pair
		real(f),intent(inout):: ardvxdt(dim,ntotal+nvirt)
		real(f):: dx(dim),piv(dim),muv,vr,rr,h,mrho
		integer:: d
		
		dx(:) = pair%i%x(:) - pair%j%x(:)
		vr = MIN(0d0,DOT_PRODUCT(pair%i%vx(:) - pair%j%vx(:),dx(:)))
		
		!Artificial viscous force only if v_ij * r_ij < 0
		rr = DOT_PRODUCT(dx(:),dx(:))
		muv  = hsml*vr/(rr + hsml*hsml*etq*etq)
		mrho = 0.5d0*(pair%i%rho + pair%j%rho)
		piv  = (beta*muv - alpha*c)*muv/mrho*pair%dwdx(:)
		
		do d = 1,dim
			!!$OMP ATOMIC
			ardvxdt(d,pair%i%ind) = ardvxdt(d,pair%i%ind) - mass*piv(d)
			!!$OMP ATOMIC
			ardvxdt(d,pair%j%ind) = ardvxdt(d,pair%j%ind) + mass*piv(d)
		end do
	
	end subroutine art_visc
	
	!=================================================================================
	subroutine ext_force(ki,pair,exdvxdt)
		
		use param, only: p1,p2,rr0,dd
	
		implicit none
		integer,intent(in):: ki
		type(interactions),intent(in):: pair
		real(f),intent(inout):: exdvxdt(dim,ntotal+nvirt)
		real(f):: dx(dim),rr,f
		integer:: d
		
		dx(:) = pair%i%x(:) - pair%j%x(:)
		rr = SQRT(SUM(dx(:)*dx(:)))
		
		if (rr.lt.rr0) then
			f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
			do d = 1,dim
				!!$OMP ATOMIC
				exdvxdt(d,pair%i%ind) = exdvxdt(d,pair%i%ind) + dd*dx(d)*f
				!!$OMP ATOMIC
				exdvxdt(d,pair%j%ind) = exdvxdt(d,pair%j%ind) - dd*dx(d)*f
			end do
		endif
	
	end subroutine ext_force
	
	!=================================================================================
	subroutine int_force(ki,pair,indvxdt)
	
		implicit none
		integer,intent(in):: ki
		type(interactions),intent(in):: pair
		real(f),intent(inout):: indvxdt(dim,ntotal+nvirt)
		real(f):: h(dim)
		integer:: d
		
		h = -(pair%i%p/pair%i%rho**2 + pair%j%p/pair%j%rho**2)*pair%dwdx(:)
		do d = 1,dim
			!!$OMP ATOMIC
			indvxdt(d,pair%i%ind) = indvxdt(d,pair%i%ind) + mass*h(d)
			!!$OMP ATOMIC
			indvxdt(d,pair%j%ind) = indvxdt(d,pair%j%ind) - mass*h(d)
		end do
	
	end subroutine int_force
	
	!=================================================================================
	subroutine con_density(ki,pair,codrhodt)
	
		implicit none
		integer,intent(in):: ki
		type(interactions),intent(in):: pair
		real(f),intent(inout):: codrhodt(ntotal+nvirt)
		real(f):: dvx(dim),vcc
		
		dvx(:) = pair%i%vx(:) - pair%j%vx(:)
				
		vcc = DOT_PRODUCT(dvx(:),pair%dwdx(:))
		!!$OMP ATOMIC
		codrhodt(pair%i%ind) = codrhodt(pair%i%ind) + mass*vcc
		!!$OMP ATOMIC
		codrhodt(pair%j%ind) = codrhodt(pair%j%ind) + mass*vcc
		
	end subroutine con_density

end module material_rates_m