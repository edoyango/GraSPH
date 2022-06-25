module globvar_para
! A module that contains global variables that are needed for the parallel scheme

   use datatypes, only: particles
   use mpi_f08
   use param, only: f, dim

   implicit none

   ! derived MPI tpes
   type(MPI_Datatype):: MPI_ftype, parttype, halotype, haloupdatetype

   !ORB variables ------------------------------------------------------------------------------------------------------------------
   integer, public:: maxnode, n_process_neighbour
   integer, allocatable, public:: node_cax(:), node_cut(:)

   !types --------------------------------------------------------------------------------------------------------------------------
   type, public:: binary_tree
      integer:: maxnode
      integer, allocatable:: node_cax(:)
   contains
      procedure:: allocate_tree_arrays
   end type binary_tree

   type, public:: partition_tracking
      integer:: mintstep_bn_part = HUGE(1), mintstep_bn_reorient = HUGE(1), maxtstep_bn_part = 0, maxtstep_bn_reorient = 0, &
                prev_part_tstep, prev_reorient_tstep, n_parts = 0, n_reorients = 0
   end type partition_tracking

   type, public:: neighbour_data
      integer:: pid, nphys_send, nphys_recv, nhalo_send, nhalo_recv
      real(f):: bounds(2*dim)
      type(MPI_Datatype):: halotype_indexed, haloupdatetype_indexed
      integer, allocatable:: halo_pindex(:)
      type(particles), allocatable:: PhysPackSend(:)
   contains
      procedure:: Pack_PhysPart, create_indexed_halotypes
   end type neighbour_data

contains

   !================================================================================================================================
   subroutine allocate_tree_arrays(self, numprocs)

      implicit none
      integer, intent(in):: numprocs
      class(binary_tree), intent(inout):: self
      integer:: tree_layers

      tree_layers = CEILING(LOG(DBLE(numprocs))/LOG(2d0))
      self%maxnode = 2*2**tree_layers - 1
      allocate (self%node_cax(maxnode))

   end subroutine allocate_tree_arrays

   !================================================================================================================================
   subroutine Pack_PhysPart(self, particle)

      implicit none
      type(particles), intent(in):: particle
      class(neighbour_data), intent(inout):: self

      self%nphys_send = self%nphys_send + 1
      self%PhysPackSend(self%nphys_send) = particle

   end subroutine Pack_PhysPart

   !================================================================================================================================
   subroutine create_indexed_halotypes(self)

      implicit none
      class(neighbour_data), intent(inout):: self
      integer:: halo_pindex_0(self%nhalo_send), ones1D(self%nhalo_send), ierr

      ones1D(:) = 1
      halo_pindex_0(:) = self%halo_pindex(1:self%nhalo_send) - 1

      call MPI_TYPE_INDEXED(self%nhalo_send, ones1D, halo_pindex_0, halotype, self%halotype_indexed, ierr)
      call MPI_TYPE_COMMIT(self%halotype_indexed, ierr)
      call MPI_TYPE_INDEXED(self%nhalo_send, ones1D, halo_pindex_0, haloupdatetype, self%haloupdatetype_indexed, ierr)
      call MPI_TYPE_COMMIT(self%haloupdatetype_indexed, ierr)

   end subroutine create_indexed_halotypes

end module globvar_para
