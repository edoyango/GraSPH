program SPH
    
    use globvar,             only: ntotal,nvirt,scale_k,allocateGlobalArrays,deallocateGlobalArrays
    use globvar_para,         only: ierr,procid,numprocs,MPI_ftype
    use mpi
    use param,                 only: f,skf
    use param_para,         only: CreateMPIType,Select_MPI_ftype
    
    use input_m,             only: input,virt_part
    use kernel_m,            only: kernel_k
    use summary_m,             only: preamble,time_print,print_summary
    use time_integration_m,    only: time_integration
    
    implicit none
    
    !Initializing MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,procid,ierr) !Retrieving process rank
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr) !Retrieving total number of processes
    
    MPI_ftype = Select_MPI_ftype(f) ! choosing the appropriate float precision for MPI data transfers
    
    !Creating MPI derived types for use later
    call CreateMPIType
    
    !Printing preamble to screen
    call preamble
    
    ! retrieving kernel k parameter for use in the program
    scale_k = kernel_k(skf)
    
    !Retrieving how many particles are to be generated. Hence the generate=.false.
    call input(.false.)
    call virt_part(.false.)
    
    if (procid.eq.0) write(*,'(A24,1x,I9,1x,A19)') 'Total simulation size of',ntotal,'physical particles.'
    
    !Allocating particle and interaction arrays
    call allocateGlobalArrays( )
    
    !Created fluid domain using particles
    call input(.true.)
    
    !Entering discretized time-integration loop
    call time_integration( )
    
    !Printing post-amble to terminal
    if (procid.eq.0) call time_print
    call print_summary

    !Cleaning up and ending MPI
    call deallocateGlobalArrays
    call MPI_FINALIZE(ierr)

end
