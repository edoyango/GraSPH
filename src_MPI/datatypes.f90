module datatypes
	
	use param, only: dim,f
	
	implicit none
	
	!variable to store particle data
	type particles
		integer:: indglob,itype,indloc
		real(f):: rho,p
		real(f):: x(dim),vx(dim)
	end type particles
	
	!variable to store particle interaction information
	type interactions
		type(particles),pointer:: i,j
		real(f):: w,dwdx(dim)
	end type interactions
	
end module datatypes