!---------------------------------------------------------
!     Including file for parameter,public::s and constants used
!     in the entire SPH software packages.
!---------------------------------------------------------

module param

   public ! assume everything defined in this module is accessible

   ! double or single precision (change f and tf to match)
   ! tf is precision for timing, f is for everything else
   ! (CUDA's timing subroutines only accept single precision floats)
   ! (Reduction of MPI process' times are hard-coding double precision)
   integer, parameter:: df = kind(1.d0), sf = kind(1.)
   integer, parameter:: f = df, tf = df

   ! constants: pi, g (gravity)
   real(f), parameter:: pi = 3.14159265358979323846_f, g = 9.81_f

   !dim : Dimension of the problem (1, 2 or 3)
   integer, parameter:: dim = 3

   !Smoothing kernel function
   !skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
   !    = 2, quartic spline kernel by W5 - Spline
   !    = 3, quintic spline kernel by W6 - Spline
   !    = 4, Wenland quintic C2 kernel (Dehnen & Aly 2012)
   !    = 5, Wenland quintic C4 kernel
   !    = 6, Wenland quintic C6 kernel
   !    = 7, Gauss kernel   (Gingold and Monaghan 1981)
   integer, parameter:: skf = 4

   !spacing and kernel radii parameters
   !note skf = 5 requires kappa at least 1.4, skf = 6 kappa at least 1.6 (approx values)
   real(f), parameter:: dxo = 0.5_f, kappa = 1.5_f, v_max = 22.15_f

   !material density (per particle)
   real(f), parameter:: irho = 1000._f

   !derived parameters. c: speed of sound, hsml: smoothing length, dt: time-step size, mass: mass per particle
   real(f), parameter:: c = 10._f*v_max, hsml = kappa*dxo, dt = 0.1_f*hsml/c, mass = irho*dxo**dim

   ! state equation parameter,public::s
   real(f), parameter:: rh0 = irho
   integer, parameter:: gamma = 7

   ! artificial viscosity parameters
   real(f), parameter:: alpha = 0.1_f, beta = 0.1_f, etq = 0.1_f

   ! repulsive force parameter,public::s
   real(f), parameter:: rr0 = dxo, dd = 5._f*g*25._f
   integer, parameter:: p1 = 4, p2 = 2

   character(*), parameter:: output_directory = "outputdata"
   character(*), parameter:: input_file = 'example/sph_out0000.h5'

   integer, parameter:: halotype = 100 ! an extra identifier to distinguish halo from real ones

end module param
