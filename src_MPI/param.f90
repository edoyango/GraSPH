!---------------------------------------------------------
!     Including file for parameters and constants used 
!     in the entire SPH software packages.
!---------------------------------------------------------

module param

	! double or single precision (chance f to match)
	integer,parameter,public:: df = kind(1.d0),sf = kind(1.)
	integer,parameter,public:: f = df
    
    ! constants: pi, g (gravity)
    real(f),parameter,public:: pi = 3.14159265358979323846_f, g = 9.81_f
	
	!dim : Dimension of the problem (1, 2 or 3)
	integer,parameter,public:: dim = 3
	
	!Smoothing kernel function 
	!skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
    !    = 2, quartic spline kernel by W5 - Spline
    !    = 3, quintic spline kernel by W6 - Spline
    !    = 4, Wenland quintic C2 kernel (Dehnen & Aly 2012)
    !    = 5, Wenland quintic C4 kernel
    !    = 6, Wenland quintic C6 kernel
	!    = 7, Gauss kernel   (Gingold and Monaghan 1981) 
	integer,parameter,public:: skf = 4
	
	!spacing and kernel radii parameters
	real(f),parameter,public:: dxo = 0.5_f, kappa = 1.2_f, v_max = 44.3_f

	!material density (per particle)
	real(f),parameter,public:: irho = 1000_f

	!derived parameters. c: speed of sound, hsml: smoothing length, dt: time-step size, mass: mass per particle
	real(f),parameter,public:: c = 10_f*v_max,hsml = kappa*dxo, dt = 1.5_f*hsml/c,mass=irho*dxo**dim

	integer,parameter,public:: mp = 50, np = 25, op = 50, pp = 3*mp, qp = np, rp = int(1.6*op)

	! state equation parameter,public::s
	real(f),parameter,public:: rh0 = irho
	integer,parameter,public:: gamma = 7
	
	! artificial viscosity parameters
	real(f),parameter,public:: alpha = 0.01_f, beta = 0_f, etq = 0.1_f
	
	! repulsive force parameters
	real(f),parameter,public:: rr0 = dxo,dd = 5_f*g*25_f
	integer,parameter,public:: p1=12,p2=6
	
	character(len=200),parameter,public:: output_directory = "outputdata"
	
	character(len=3),parameter,public:: output_flt_type='dbl'
	logical,parameter,public:: output_phys(2) = (/.true.,.true./)
	logical,parameter,public:: output_halo(2) = (/.true.,.false./)
	logical,parameter,public:: output_virt(2) = (/.true.,.false./)
	logical,parameter,public:: output_CPUtime = .false.
	logical,parameter,public:: output_boundary = .false.

end
