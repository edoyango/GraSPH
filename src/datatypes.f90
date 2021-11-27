module derived_types
	
	use param
	implicit none
	
	!variable to store particle data
	type particles
		integer:: itype,ind
		real(8):: rho,p
		real(8):: x(dim),vx(dim)
	end type particles
	
	!variable to store particle interaction information
	type interactions
		type(particles),pointer:: i,j
		real(8):: w,dwdx(dim)
	end type interactions
	
end module derived_types