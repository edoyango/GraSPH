!---------------------------------------------------------
!     Including file for parameters and constants used 
!     in the entire SPH software packages.
!---------------------------------------------------------

module param
	!dim : Dimension of the problem (1, 2 or 3)
	integer,parameter:: dim = 2
	
	!Smoothing kernel function 
	!skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
	!    = 2, Gauss kernel   (Gingold and Monaghan 1981) 
	integer,parameter:: skf = 1
	
	!spacing and kernel radii parameters
	real(8),parameter:: dxo = 0.2d0, kappa = 1.2d0, v_max = 22.15d0
	
	!material density (per particle)
	real(8),parameter:: irho = 1000d0
	
	!derived parameters. c: speed of sound, hsml: smoothing length, dt: time-step size, mass: mass per particle
	real(8),parameter:: c = 10d0*v_max,hsml = kappa*dxo, dt = 1.5d0*hsml/c,mass=irho*dxo**dim
	
	integer,parameter:: mp = 125, np = 125, op = 375, pp = 200
	
	! state equation parameters
	real(8),parameter:: rh0 = irho
	integer,parameter:: gamma = 7
	
	! artificial viscosity parameters
	real(8),parameter:: alpha = 0.01d0, beta = 0.d0, etq = 0.1d0
	
	! repulsive force parameters
	real(8),parameter:: rr0 = dxo,dd = 5d0*9.81d0*25d0
	integer,parameter:: p1=4,p2=2
	
	character(len=200),parameter:: output_directory = "C:\Users\edwar\Documents\outputdata\"
	
	logical,parameter:: output_phys(2) = (/.true.,.true./)
	logical,parameter:: output_virt(2) = (/.true.,.false./)

end