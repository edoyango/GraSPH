module datatypes

   use param, only: dim, f, tf

   implicit none

   private

   !variable to store particle data
   type particles
      integer:: itype, indloc, indglob
      real(f):: rho, p
      real(f):: x(dim), vx(dim)
   end type particles

   !variable to store particle interaction information
   type interactions
      integer:: j
      real(f):: w, dwdx(dim)
   end type interactions

   ! Data type to store timing variables
   type time_tracking
      real(tf):: t_wall = 0._tf, t_ORB = 0._tf, t_dist = 0._tf, t_output = 0._tf
   contains
      procedure:: t_compute
   end type

   public:: particles, interactions, time_tracking, system_clock_timer

contains

   pure real(tf) function t_compute(self)

      class(time_tracking), intent(in):: self

      t_compute = self%t_wall - self%t_ORB - self%t_dist - self%t_output

   end function t_compute

   function system_clock_timer() result(t)

      use, intrinsic:: ISO_FORTRAN_ENV, only: int64

      real(tf):: t
      integer(int64):: c, cr

      call SYSTEM_CLOCK(c, cr)

      t = dble(c)/dble(cr)

   end function system_clock_timer

end module datatypes
