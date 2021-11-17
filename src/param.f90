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
	
	!Control parameters for output 
	!print_step: Print Timestep (On Screen)
	!save_step : Save Timestep  (To Disk File)
	integer,parameter:: print_step = 100, save_step = 100
	real(8),parameter:: pi = 3.14159265358979323846d0, g = 9.81d0, dxo = 0.5d0, kappa = 1.2d0, v_max = 22.15d0
	
	!material properties
	real(8),parameter:: irho = 1000d0
	
	!derived parameters
	real(8),parameter:: c = 10d0*v_max,hsml = kappa*dxo, dt = 1.5d0*hsml/c,mass=irho*dxo**dim
	
	!physical particle properties
	type particles
		integer:: itype,ind
		real(8):: rho,p
		real(8):: x(dim),vx(dim)
	end type particles
	
	type interactions
		type(particles),pointer:: i,j
		real(8):: w,dwdx(dim)
	end type interactions
	
	type(particles),allocatable,target:: parts(:)
	type(interactions),allocatable:: pairs(:)
	
	!global 1D variables
	integer:: ntotal,nvirt
	integer:: mp,np,op,pp,maxn,maxinter,maxnsend
	integer:: niac
	integer:: nstart,current_ts,maxtimestep,itimestep,yesorno
	real(8):: cputime,s1,s2,time,testtime(20)
	real(8):: scale_k
	
	!diagnostics
	integer:: nvalid

end