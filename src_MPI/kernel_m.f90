module kernel_m
    
    use param, only: f
    
	public:: kernel,kernel_k

contains

	!==============================================================================================================================
	subroutine kernel(r,dx,thsml,tw,tdwdx)
	! Contains the kernels
	
		use param, only: dim,skf,pi
		
		implicit none
		real(f),intent(in):: dx(dim),r,thsml
		real(f),intent(out):: tdwdx(dim),tw
		real(f):: q,factor
		
		q = r/thsml
		tw = 0_f
		tdwdx(:) = 0_f
		
		SELECT CASE (SKF)
			CASE (1)
				factor = 10_f/(7_f*pi*thsml*thsml)
				if (q.ge.0_f .and. q.lt.1_f) then          
					tw = factor * (1_f - 1.5_f*q*q + 0.75_f*q**3)
					tdwdx(:) = factor * (-3_f*q + 2.25_f*q*q) * dx(:) / (thsml*r)
				else if (q.ge.1_f .and. q.lt.2_f) then
					tw = factor * 0.25_f * (2_f-q)**3_f
					tdwdx(:) = -factor*0.75_f*((2_f-q)**2)*dx(:) / (thsml*r)        
				endif
			CASE (2)
		
				factor = 15_f/(7_f*pi*thsml*thsml)
				if (q.ge.0_f .and. q.lt.1_f) then          
					tw = factor * (2_f/3_f - q*q + 0.5_f*q**3)
					tdwdx(:) = factor * (-2_f + 1.5_f * q) / thsml**2 * dx(:)       
				else if (q.ge.1_f .and. q.lt.2_f) then
					tw = factor/6_f * (2_f-q)**3
					tdwdx(:) =-factor/6_f*3_f*(2_f-q)**2/thsml*(dx(:)/r)        
				endif
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
				scale_k = 2_f
			case (2)
				scale_k = 3_f
			case (3)
				scale_k = 2_f
			case (4)
				scale_k = 2_f
		end select
	
	end function kernel_k
	
end module kernel_m
