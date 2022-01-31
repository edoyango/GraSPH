module material_rates_m
	
	use constants
	use globvar
	use param

contains
	
	!=================================================================================
	subroutine con_density(ki,cdrhodt)
	! calculates the densiy rate of change via the SPH approximation of the continuity equation
	
		implicit none
		integer,intent(in):: ki
		real(8),intent(out):: cdrhodt(ntotal+nvirt)
		integer:: i,j,k,d
		real(8):: dvx(dim),vcc
		
		cdrhodt(1:ntotal) = 0d0
		
		do k=1,niac      
		
			if (pairs(k)%i%itype.gt.0 .and. pairs(k)%j%itype.gt.0) then
			
				dvx(:) = pairs(k)%i%vx(:) - pairs(k)%j%vx(:) 
				
				vcc = DOT_PRODUCT(dvx(:),pairs(k)%dwdx(:))
				
				cdrhodt(pairs(k)%i%ind) = cdrhodt(pairs(k)%i%ind) + mass*vcc
				cdrhodt(pairs(k)%j%ind) = cdrhodt(pairs(k)%j%ind) + mass*vcc
				
			endif
		
		enddo
		
	end subroutine con_density
	
	!======================================================================================
	subroutine art_visc(ki,ardvxdt)
	! calculates the acceleration acting on each particle as a result of artificial viscosity (Monaghan 1994 form)
	
		implicit none
		integer,intent(in):: ki
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
	
	end subroutine art_visc
	
	!==============================================================================================
	subroutine int_force(ki,indvxdt)
	! calculates the acceleration acting on each particle as a result of pressure via the SPH approximation of the continuity
	! equation (conservative form)
	
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
	
	!=======================================================================
	subroutine ext_force(ki,exdvxdt)
	! Supplied acceleration acting on each particle due to external forces e.g. gravity, artifical lennard-jones boundary force.
		
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

end module material_rates_m