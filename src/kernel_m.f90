module kernel_m
	
	public:: kernel,kernel_k
	
contains
	
	!==============================================================================================================================
	subroutine kernel(r,dx,thsml,tw,tdwdx)   
	! Contains the kernels
	
		use constants,	only: pi
		use param, 		only: f,skf,dim
	
		implicit none
		real(f),intent(in):: dx(dim),r,thsml
		real(f),intent(out):: tdwdx(dim),tw
		real(f):: q,factor
		
		q = r/thsml
		tw = 0.d0
		tdwdx(:) = 0.d0
		
		SELECT CASE (SKF)
			CASE (1) ! cubic
				select case (dim)
					case(2)
						factor = 10.d0/(7.d0*pi*thsml*thsml)
					case(3)
						factor = 1d0/(pi*thsml*thsml)
				end select
				if (q.ge.0.d0 .and. q.lt.1.d0) then          
					tw = factor * (1.d0 - 1.5d0*q*q + 0.75d0*q**3)
					tdwdx(:) = factor * (-3.d0*q + 2.25d0*q*q) * dx(:) / (thsml*r)
				else if (q.ge.1.d0 .and. q.lt.2.d0) then
					tw = factor * 0.25d0 * (2.d0-q)**3.d0
					tdwdx(:) = -factor*0.75d0*((2.d0-q)**2)*dx(:) / (thsml*r)        
				end if
			CASE (2) ! gaussian
				factor = 1d0/(thsml**dim*pi**(0.5d0*dim))
				if (q.ge.0d0 .and. q.le.3d0) then
					tw = factor*exp(-q*q)
					tdwdx(:) = tw*2d0*dx(:)/(thsml*thsml)
				end if
		END SELECT
		
	end subroutine kernel
	
	!==============================================================================================================================
	function kernel_k(skf) result(scale_k)
	! setting k parameter for kernel radius (r = kh)
	
		implicit none
		integer,intent(in):: skf
		integer:: scale_k
		
		select case (skf)
			case (1)
				scale_k = 2d0
			case (2)
				scale_k = 3d0
			case (3)
				scale_k = 2d0
			case (4)
				scale_k = 2d0
		end select
	
	end function kernel_k
	
end module kernel_m