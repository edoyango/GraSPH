module globvar_para
! A module that contains global variables that are needed for the parallel scheme

   use datatypes, only: particles
   use mpi_f08
   use param, only: f, dim

   implicit none

   ! basic MPI variables
   integer, public:: ierr

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

end module globvar_para
