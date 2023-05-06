module datatypes

   use param, only: dim, f, tf

   implicit none

   private

   !variable to store particle data
   type particles
      ! integer:: itype, indloc, indglob
      ! real(f):: rho, rho_min, p
      ! real(f):: x(dim), vx(dim), v_min(dim)
      integer:: ntotal, nvirt, ntotal_loc, nvirt_loc, nhalo_loc, maxn
      integer, allocatable:: itype(:), indloc(:), indglob(:)
      real(f), allocatable:: rho(:), rho_min(:), p(:)
      real(f), allocatable:: x(:, :), vx(:, :), v_min(:, :)
   contains
      procedure:: allocate_particles, deallocate_particles, nsum
   end type particles

   !variable to store particle interaction information
   type interactions
      integer:: i, j
      real(f):: dx(dim)
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

   subroutine allocate_particles(self)

      implicit none
      class(particles), intent(inout):: self

      allocate(self%itype(self%maxn), self%indloc(self%maxn), self%indglob(self%maxn))
      allocate(self%rho(self%maxn), self%rho_min(self%maxn), self%p(self%maxn))
      allocate(self%x(dim, self%maxn), self%vx(dim, self%maxn), self%v_min(dim, self%maxn))

   end subroutine allocate_particles

   subroutine deallocate_particles(self)

      implicit none
      class(particles), intent(inout):: self

      deallocate(self%itype, self%indloc, self%indglob)
      deallocate(self%rho, self%rho_min, self%p)
      deallocate(self%x, self%vx, self%v_min)

   end subroutine deallocate_particles

   pure integer function nsum(self)

      implicit none
      class(particles), intent(in):: self

      nsum = self%ntotal_loc + self%nhalo_loc + self%nvirt_loc
      
   end function nsum

end module datatypes
