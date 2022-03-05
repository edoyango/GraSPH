!---------------------------------------------------------
!     Including file for parameter,public::s and constants used 
!     in the entire SPH software packages.
!---------------------------------------------------------

module param

	use constants
	
	! double or single precision (chance f to match)
	integer,parameter,public:: df = kind(1.d0)
	integer,parameter,public:: sf = kind(1.)
	integer,parameter,public:: f = df

	!dim : Dimension of the problem (1, 2 or 3)
	integer,parameter,public:: dim = 3
	
	!Smoothing kernel function 
	!skf = 1, cubic spline kernel by W4 - Spline (Monaghan 19f5)
	!    = 2, Gauss kernel   (Gingold and Monaghan 19f1) 
	integer,parameter,public:: skf = 1
	
	!spacing and kernel radii parameter,public::s
	real(f),parameter,public:: dxo = 0.5d0, kappa = 1.2d0, v_max = 22.15d0
	
	!material density (per particle)
	real(f),parameter,public:: irho = 1000d0
	
	!derived parameter,public::s. c: speed of sound, hsml: smoothing length, dt: time-step size, mass: mass per particle
	real(f),parameter,public:: c = 10d0*v_max,hsml = kappa*dxo, dt = 1.5d0*hsml/c,mass=irho*dxo**dim
	
	integer,parameter,public:: mp = 50, np = 10, op = 50, pp = 3*mp, qp = np, rp = 1.6*op
	
	! state equation parameter,public::s
	real(f),parameter,public:: rh0 = irho
	integer,parameter,public:: gamma = 7
	
	! artificial viscosity parameter,public::s
	real(f),parameter,public:: alpha = 0.01d0, beta = 0.d0, etq = 0.1d0
	
	! repulsive force parameter,public::s
	real(f),parameter,public:: rr0 = dxo,dd = 5d0*g*25d0
	integer,parameter,public:: p1=4,p2=2
	
	character(len=200),parameter,public:: output_directory = "outputdata"
	
	logical,parameter,public:: output_phys(2) = (/.true.,.true./)
	logical,parameter,public:: output_virt(2) = (/.true.,.false./)

end module param
