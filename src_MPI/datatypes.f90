module datatypes

   use param, only: dim, f

   implicit none

   !variable to store particle data
   type particles
      integer:: indglob, itype, indloc
      real(f):: rho, p
      real(f):: x(dim), vx(dim)
   end type particles

   !variable to store particle interaction information
   type interactions
      integer:: j
      real(f):: w, dwdx(dim)
   end type interactions

   type time_tracking
      double precision:: t_wall = 0.d0, t_ORB = 0.d0, t_dist = 0.d0, t_output = 0.d0
   contains
      procedure:: t_compute
   end type

contains

   pure double precision function t_compute(self)

      class(time_tracking), intent(in):: self

      t_compute = self%t_wall - self%t_ORB - self%t_dist - self%t_output

   end function t_compute

end module datatypes
