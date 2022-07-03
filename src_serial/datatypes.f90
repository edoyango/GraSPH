module datatypes

   use param, only: dim, f

   implicit none

   private

   !variable to store particle data
   type particles
      integer:: itype, ind
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
      double precision:: t_compute = 0.d0, t_output = 0.d0
   contains
      procedure:: walltime
   end type

   public:: particles, interactions, time_tracking, system_clock_timer

contains

   pure double precision function walltime(self)

      class(time_tracking), intent(in):: self

      walltime = self%t_compute + self%t_output

   end function walltime

   function system_clock_timer() result(t)

      use, intrinsic:: ISO_FORTRAN_ENV, only: int64

      double precision:: t
      integer(int64):: c, cr

      call SYSTEM_CLOCK(c, cr)

      t = dble(c)/dble(cr)

   end function system_clock_timer

end module datatypes
