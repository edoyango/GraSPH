! paramaters that relate to the parallel scheme for this code e.g. array sizes, subdomain boundary
! update frequency, number of msgs etc.
module param_para

   use datatypes, only: particles
   use mpi_f08
   use param, only: dim, f

   private
! parameters for ORB algorithm
! dcell_ORB = size of each cell used to draw subdomain boundaries (in multiples of hsml)
!             smaller -> better partition, but slower and requires more memory.
!             bigger -> poorer partition, but faster and requires less memory.
! bound_extend = the distance of which to extend the current global domain boundaries (in multiples of kernel radii). Choose a
!                large value if you expect global extents of your simulation to change a lot over the simulation.
! box_ratio_threshold = controls the required proportional change in aspect ratio before subdomain bounadry orientations are
!                       re-determined.
!                       Smaller value -> more frequent axes reorientations of cuts.
! ORBcheck1, ORBcheck2 = how frequent to check if boundaries need updating. Lets say the boundar is updated at time-step N, the
!                        next check will be at time-step N + ORBcheck1. If the boundaries are not updated, checks will occur every
!                        ORBcheck2 time-steps thereafter. Too small and program will be very slow. Too large, and load balance will
!                        be very poor in general. The more slowly the particles move, ORBchecks can be less frequent.
   real(f), parameter, public:: dcell_ORB = 1_f, box_ratio_threshold = 0.25_f, bound_extend = 10_f
   integer, parameter, public:: ORBcheck1 = 15, ORBcheck2 = 7
   
   !types --------------------------------------------------------------------------------------------------------------------------
   ! type to hold all MPI derived types used in simulation
   type, public:: MPI_derived_types
      type(MPI_Datatype):: ftype, parttype, halotype, haloupdatetype
   contains
      procedure:: CreateMPITypes
   end type MPI_derived_types
   
   ! type to hold variables related to partitioning
   type, public:: partition_tracking
      integer:: mintstep_bn_part = HUGE(1), mintstep_bn_reorient = HUGE(1), maxtstep_bn_part = 0, maxtstep_bn_reorient = 0, &
                prev_part_tstep, prev_reorient_tstep, n_parts = 0, n_reorients = 0
   end type partition_tracking

   ! type to hold all variables related to data transfer between neighbouring processes
   type, public:: neighbour_data
      integer:: pid, nphys_send, nphys_recv, nhalo_send, nhalo_recv
      real(f):: bounds(2*dim)
      type(MPI_Datatype):: halotype_indexed, haloupdatetype_indexed
      integer, allocatable:: halo_pindex(:)
      type(particles), allocatable:: PhysPackSend(:)
   contains
      procedure:: Pack_PhysPart, create_indexed_halotypes
   end type neighbour_data

! =================================================================================================
! FUNCTIONS FOR PACKING DATA FOR PARTICLE EXCHANGE SUBROUTINES
!==================================================================================================
contains
   !================================================================================================================================
   subroutine CreateMPITypes(self) !---------------------------------
      ! Subroutine to help create MPI user (derived) types for particle data. four types are created in the parent class:
      ! - ftype: the floating type precision to used
      ! - parttype: an MPI user type to encapsulate physical particle data to be exchanged
      ! - halotype: "                             " halo "                               "
      ! - haloupdatetype: an MPI user type that captures a reduced subset of particle data as only some variables need to be
      !                   updated after the first halo exchange in a time-step

      implicit none
      class(MPI_derived_types), intent(out):: self
      type(particles):: parts_dummy(2)
      integer:: blockl(7), ierr
      type(MPI_Datatype):: type(7)
      integer(KIND=MPI_ADDRESS_KIND):: disp(7), basedisp
      integer, parameter:: sp = kind(1.)

      ! Selecting appropriate MPI float datatype
      select case (f)
      case default
         self%ftype = MPI_DOUBLE_PRECISION
      case (sp)
         self%ftype = MPI_REAL
      end select

      ! Creating derived type for all particles ------------------------------------------------------------------------------------
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
      type(4:7) = self%ftype

      ! Creating MPI user type
      self%parttype = CreateMPITYPES_helper(7, disp, blockl, type)

      ! Repeating process for halo exchanges (only subset) -------------------------------------------------------------------------
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
      type(3:5) = self%ftype

      self%halotype = CreateMPITypes_helper(5, disp, blockl, type)

      ! Setting up type for halo particle updates (updates only need density, velocity) --------------------------------------------
      call MPI_GET_ADDRESS(parts_dummy(1), basedisp, ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%rho, disp(1), ierr)
      call MPI_GET_ADDRESS(parts_dummy(1)%vx(1), disp(2), ierr)

      disp(:) = disp(:) - basedisp

      blockl(1:2) = (/1, dim/)

      type(1:2) = self%ftype

      self%haloupdatetype = CreateMPITypes_helper(2, disp, blockl, type)

   end subroutine CreateMPITypes

   !================================================================================================================================
   function CreateMPITypes_helper(nblock, disp, blockl, blockt) result(outtype)
      ! This function is a helper function that encapsulates all the MPI_ calls to create a datatype. Is used only in the
      ! CreateMPITypes class subroutine

      implicit none
      integer, intent(in):: nblock, blockl(nblock)
      integer(KIND=MPI_ADDRESS_KIND), intent(in):: disp(nblock)
      type(MPI_Datatype), intent(in):: blockt(nblock)
      integer(KIND=MPI_ADDRESS_KIND):: lb = 0, disp_tmp(2), ext
      type(MPI_Datatype):: tmptype, outtype
      integer:: ierr
      type(particles):: parts_dummy(2)

      call MPI_TYPE_CREATE_STRUCT(nblock, blockl, disp, blockt, tmptype, ierr)
      call MPI_TYPE_COMMIT(tmptype, ierr)

      call MPI_GET_ADDRESS(parts_dummy(1), disp_tmp(1), ierr)
      call MPI_GET_ADDRESS(parts_dummy(2), disp_tmp(2), ierr)
      lb = 0; ext = disp_tmp(2) - disp_tmp(1)
      call MPI_TYPE_CREATE_RESIZED(tmptype, lb, ext, outtype, ierr)
      call MPI_TYPE_COMMIT(outtype, ierr)
      call MPI_TYPE_FREE(tmptype, ierr)

   end function CreateMPITypes_helper
   
   !================================================================================================================================
   subroutine Pack_PhysPart(self, particle)

      implicit none
      type(particles), intent(in):: particle
      class(neighbour_data), intent(inout):: self

      self%nphys_send = self%nphys_send + 1
      self%PhysPackSend(self%nphys_send) = particle

   end subroutine Pack_PhysPart

   !================================================================================================================================
   subroutine create_indexed_halotypes(self, halotype, haloupdatetype)

      implicit none
      type(MPI_Datatype), intent(in):: halotype, haloupdatetype
      class(neighbour_data), intent(inout):: self
      integer:: halo_pindex_0(self%nhalo_send), ones1D(self%nhalo_send), ierr

      ones1D(:) = 1
      halo_pindex_0(:) = self%halo_pindex(1:self%nhalo_send) - 1

      call MPI_TYPE_INDEXED(self%nhalo_send, ones1D, halo_pindex_0, halotype, self%halotype_indexed, ierr)
      call MPI_TYPE_COMMIT(self%halotype_indexed, ierr)
      call MPI_TYPE_INDEXED(self%nhalo_send, ones1D, halo_pindex_0, haloupdatetype, self%haloupdatetype_indexed, ierr)
      call MPI_TYPE_COMMIT(self%haloupdatetype_indexed, ierr)

   end subroutine create_indexed_halotypes

end module param_para
