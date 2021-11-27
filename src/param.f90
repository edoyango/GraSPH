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
	
	!Period of data writing (time steps).
	!print_step: Print summary data to terminal
	!save_step : Write output data to file
	integer,parameter:: print_step = 100, save_step = 100
	
	!spacing and kernel radii parameters
	real(8),parameter:: dxo = 0.5d0, kappa = 1.2d0, v_max = 22.15d0
	
	!material density (per particle)
	real(8),parameter:: irho = 1000d0
	
	!derived parameters. c: speed of sound, hsml: smoothing length, dt: time-step size, mass: mass per particle
	real(8),parameter:: c = 10d0*v_max,hsml = kappa*dxo, dt = 1.5d0*hsml/c,mass=irho*dxo**dim

end