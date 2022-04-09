module datatypes
	
	use param, only: dim,f,tenselem
	
	implicit none
	
	!variable to store particle data
	type particles
		integer:: itype,ind
		real(f):: rho,p
		real(f):: x(dim),vx(dim)
        real(f):: strain(tenselem),pstrain(tenselem),sig(tenselem)
	end type particles
	
	!variable to store particle interaction information
	type interactions
		type(particles),pointer:: i,j
		real(f):: w,dwdx(dim)
	end type interactions
	
	public:: particles,interactions
	
end module datatypes
