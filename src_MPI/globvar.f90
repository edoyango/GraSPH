module globvar

   use param, only: f
   use datatypes, only: particles, interactions

   implicit none

   ! particle array
   type(particles), allocatable, public:: parts(:)

   ! interaction array
   type(interactions), allocatable, public:: pairs(:)

   !global 1D variables
   integer, public:: ntotal, nvirt, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc
   integer, public:: maxn, maxinter, maxnloc
   integer, public:: niac
   integer, public:: itimestep, maxtimestep, save_step, print_step
   real(f), public:: time

   real(f), public:: scale_k

   !timing
   real(f), public:: t_graph = 0_f, t_dist = 0_f, cputime = 0_f, output_time = 0_f, test_time = 0_f

   public:: allocateGlobalArrays, deallocateGlobalArrays

! subroutines to allocate and deallocate global arrays
contains

   !==============================================================================================================================
   subroutine allocateGlobalArrays

      maxnloc = 2*ntotal + nvirt
      maxinter = 262*maxnloc

      allocate (parts(maxnloc))
      allocate (pairs(maxinter))

   end subroutine allocateGlobalArrays

   !==============================================================================================================================
   subroutine deallocateGlobalArrays

      deallocate (parts)
      deallocate (pairs)

   end subroutine deallocateGlobalArrays

end module globvar

