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
      integer:: i, j
      real(f):: w, dwdx(dim)
   end type interactions

   type time_tracking
      double precision:: t_compute = 0.d0, t_ORB = 0.d0, t_dist = 0.d0, t_output = 0.d0
   contains
      procedure:: walltime
   end type

contains

   pure double precision function walltime(self)

      class(time_tracking), intent(in):: self

      walltime = self%t_compute + self%t_ORB + self%t_dist + self%t_output

   end function walltime

end module datatypes
