module globvar_para
! A module that contains global variables that are needed for the parallel scheme

   use datatypes, only: particles
   use mpi_f08
   use param, only: f, dim

   implicit none

   ! derived MPI tpes
   type(MPI_Datatype):: MPI_ftype, parttype, halotype, haloupdatetype
   type(MPI_Datatype), allocatable, public:: halotype_indexed(:), haloupdatetype_indexed(:)

   ! particle send/recv arrays -------------------------------------------------------------------------------------------------------
   integer, allocatable, public:: nphys_send(:), nphys_recv(:)
   integer, allocatable, public:: nhalo_send(:), nhalo_recv(:)
   type(particles), allocatable, public:: PhysPackSend(:, :)

   !ORB variables --------------------------------------------------------------------------------------------------------------------
   integer, public:: maxnode, n_process_neighbour, leaf_node, repartition_mode
   integer, allocatable, public:: node_cax(:), proc_neighbour_list(:), node_cut(:), &
                                  halo_pindex(:, :)
   real(f), allocatable, public:: bounds_glob(:, :)
   
   type, public:: binary_tree
      integer:: maxnode
      integer,allocatable:: node_cax(:)
      real(f),allocatable:: bounds_glob(:,:)
      contains
         procedure:: allocate_tree_arrays
   end type binary_tree
   
   type, public:: partition_tracking
      integer:: mintstep_bn_part = HUGE(1), mintstep_bn_reorient = HUGE(1), maxtstep_bn_part = 0, maxtstep_bn_reorient = 0, &
                     prev_part_tstep, prev_reorient_tstep, n_parts = 0, n_reorients = 0
   end type partition_tracking
   
contains

   subroutine allocate_tree_arrays(self,numprocs)

      implicit none
      integer,intent(in):: numprocs
      class(binary_tree),intent(inout):: self
      integer:: tree_layers
      
      tree_layers = CEILING(LOG(DBLE(numprocs))/LOG(2d0))
      self%maxnode = 2*2**tree_layers - 1
      allocate(self%bounds_glob(2*dim,numprocs),self%node_cax(maxnode))
      
   end subroutine allocate_tree_arrays
      
end module globvar_para
