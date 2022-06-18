! paramaters that relate to the parallel scheme for this code e.g. array sizes, subdomain boundary
! update frequency, number of msgs etc.
module param_para

   use datatypes, only: particles
   use globvar_para, only: parttype, halotype, haloupdatetype, ierr, MPI_ftype
   use mpi_f08
   use param, only: dim, f

   private
! parameters for ORB algorithm
! dcell_ORB = size of each cell used to draw subdomain boundaries (in multiples of hsml)
!             smaller -> better partition, but slower and requires more memory.
!             bigger -> poorer partition, but faster and requires less memory.
! bound_extend = the distance of which to extend the current global domain boundaries (in multiples
!                of kernel radii). Choose a large value if you expect global extents of your
!                simulation to change a lot over the simulation.
! box_ratio_threshold = controls the required proportional change in aspect ratio before subdomain
!                       bounadry orientations are re-determined.
!                       Smaller value -> more frequent axes reorientations of cuts.
! ORBcheck1, ORBcheck2 = how frequent to check if boundaries need updating. Lets say the boundary
!                        is updated at time-step N, the next check will be at time-step
!                        N + ORBcheck1. If the boundaries are not updated, checks will occur every
!                        ORBcheck2 time-steps thereafter.
!                        Too small and program will be very slow. Too large, and load balance will
!                        be very poor in general. The more slowly the particles move, ORBchecks
!                        can be less frequent
   real(f), parameter, public:: dcell_ORB = 1_f, box_ratio_threshold = 0.25_f, bound_extend = 10_f
   integer, parameter, public:: ORBcheck1 = 15, ORBcheck2 = 7

   public:: CreateMPIType, Select_MPI_ftype

! =================================================================================================
! FUNCTIONS FOR PACKING DATA FOR PARTICLE EXCHANGE SUBROUTINES
!==================================================================================================
contains

   subroutine CreateMPIType !---------------------------------
      ! packing of physical particle data (integers and doubles)

      implicit none
      type(particles):: parts_dummy(2)
      integer:: blockl(7)
      type(MPI_Datatype):: type(7), tmptype
      integer(KIND=MPI_ADDRESS_KIND):: disp(7), lb, ext, basedisp

      ! Obtaining memory address of each block in derived type
      call MPI_GET_ADDRESS(parts_dummy(1), basedisp, ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%indglob, disp(1), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%indloc, disp(2), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%itype, disp(3), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%rho, disp(4), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%p, disp(5), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%x(1), disp(6), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%vx(1), disp(7), ierr)

      ! Converting absolute addressed to relative address (relative to indglob)
      disp(:) = disp(:) - basedisp

      ! Size of each block in derived type
      blockl(:) = (/1, 1, 1, 1, 1, dim, dim/)

      ! Types of each block (MPI)
      type(1:3) = MPI_INTEGER
      type(4:7) = MPI_ftype

      ! Creating MPI user type
      call MPI_TYPE_CREATE_STRUCT(7, blockl, disp, type, tmptype, ierr)
      call MPI_TYPE_COMMIT(tmptype, ierr)

      ! Calculating distance in memory between consecutive elements in array of defined type
      call MPI_GET_ADDRESS(parts_dummy(1), disp(1), ierr)
      call MPI_GET_ADDRESS(parts_dummy(2), disp(2), ierr)
      ! Resizing array to account for any "padding". lb = lowerbound
      lb = 0; ext = disp(2) - disp(1)
      call MPI_TYPE_CREATE_RESIZED(tmptype, lb, ext, parttype, ierr)
      call MPI_TYPE_COMMIT(parttype, ierr)
      call MPI_TYPE_FREE(tmptype, ierr)

      ! Repeating process for halo exchanges
      ! Setting up type for initial halo particle exchange
      call MPI_GET_ADDRESS(parts_dummy(1), basedisp, ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%indglob, disp(1), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%itype, disp(2), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%rho, disp(3), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%x(1), disp(4), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%vx(1), disp(5), ierr)

      disp(:) = disp(:) - basedisp

      blockl(1:5) = (/1, 1, 1, dim, dim/)

      type(1:2) = MPI_INTEGER
      type(3:5) = MPI_ftype

      call MPI_TYPE_CREATE_STRUCT(5, blockl, disp, type, tmptype, ierr)
      call MPI_TYPE_COMMIT(tmptype, ierr)

      call MPI_GET_ADDRESS(parts_dummy(1), disp(1), ierr)
      call MPI_GET_ADDRESS(parts_dummy(2), disp(2), ierr)
      lb = 0; ext = disp(2) - disp(1)
      call MPI_TYPE_CREATE_RESIZED(tmptype, lb, ext, halotype, ierr)
      call MPI_TYPE_COMMIT(halotype, ierr)
      call MPI_TYPE_FREE(tmptype, ierr)

      ! Setting up type for halo particle updates (updates only need density, velocity)
      call MPI_GET_ADDRESS(parts_dummy(1), basedisp, ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%rho, disp(1), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%vx(1), disp(2), ierr)

      disp(:) = disp(:) - basedisp

      blockl(1:2) = (/1, dim/)

      type(1:2) = MPI_ftype

      call MPI_TYPE_CREATE_STRUCT(2, blockl, disp, type, tmptype, ierr)
      call MPI_TYPE_COMMIT(tmptype, ierr)

      call MPI_GET_ADDRESS(parts_dummy(1)%rho, disp(1), ierr)
      call MPI_GET_ADDRESS(parts_dummy(2)%rho, disp(2), ierr)
      lb = 0; ext = disp(2) - disp(1)
      call MPI_TYPE_CREATE_RESIZED(tmptype, lb, ext, haloupdatetype, ierr)
      call MPI_TYPE_COMMIT(haloupdatetype, ierr)
      call MPI_TYPE_FREE(tmptype, ierr)

   end subroutine CreateMPIType !--------------------------------------------------------------------

   function Select_MPI_ftype(ftype) result(MPI_ftype)

      implicit none
      integer, intent(in):: ftype
      type(MPI_Datatype):: MPI_ftype
      integer, parameter:: dp = kind(1_f), sp = kind(1.)

      select case (ftype)
      case (dp)
         MPI_ftype = MPI_DOUBLE_PRECISION
      case (sp)
         MPI_ftype = MPI_REAL
      end select

   end function Select_MPI_ftype

end module param_para
