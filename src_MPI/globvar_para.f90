module globvar_para
! A module that contains global variables that are needed for the parallel scheme
    
    use datatypes,    only: particles
    use param,        only: f,dim
    
    implicit none
    
    ! basic MPI variables
    integer,public:: procid,numprocs,ierr
    
    ! derived MPI tpes
    integer,public:: parttype,halotype,haloupdatetype,MPI_ftype
    integer,allocatable,public:: halotype_indexed(:),haloupdatetype_indexed(:)

    ! particle send/recv arrays -------------------------------------------------------------------------------------------------------
    integer,allocatable,public:: nphys_send(:),nphys_recv(:)
    integer,allocatable,public:: nhalo_send(:),nhalo_recv(:)
    type(particles),allocatable,public:: PhysPackSend(:,:)
    
    !ORB variables --------------------------------------------------------------------------------------------------------------------
    integer,public:: maxnode,n_process_neighbour,leaf_node,repartition_mode
    integer,allocatable,public:: node_cax(:),pincell_ORB(:,:,:),proc_neighbour_list(:),node_cut(:),&
    node_segment(:),halo_pindex(:,:)
    real(f),allocatable,public:: bounds_glob(:,:)
    
    !Partition frequency variables
    integer,public:: prev_part_tstep,mintstep_bn_part=HUGE(1),maxtstep_bn_part=0,n_parts=0
    integer,public:: prev_reorient_tstep,mintstep_bn_reorient=HUGE(1),maxtstep_bn_reorient=0,n_reorients=0
    integer,public:: prev_load
    real(f),public:: box_ratio_previous(dim,dim) = TINY(1.)

end module globvar_para
