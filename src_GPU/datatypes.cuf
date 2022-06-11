module datatypes

   use param, only: dim, f

   implicit none

   !variable to store particle data
   type particles
      integer:: itype, ind
      real(f):: rho, p
      real(f):: x(dim), vx(dim)
   end type particles

   !variable to store particle interaction information
   type interactions
      integer:: i, j
      real(f):: w, dwdx(dim)
   end type interactions

   public:: particles, interactions

end module datatypes
