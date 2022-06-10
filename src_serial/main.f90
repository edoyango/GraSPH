program SPH

    use globvar, only: allocateGlobalArrays,deallocateGlobalArrays,scale_k,ntotal,nvirt
    use param, only: skf
    
    use input_m, only: input,virt_part
    use kernel_m, only: kernel_k
    use summary_m, only: preamble,time_print,print_summary
    use time_integration_m, only: time_integration

    implicit none
    
    !Printing preamble to screen
    call preamble
    
    ! setting k parameter for kernel radius (r = kh)
    scale_k = kernel_k(skf)
    
    call input(.false.)
    call virt_part(.false.)
    
    write(*,'(A24,1x,I9,1x,A19)') 'Total simulation size of',ntotal,'physical particles.'
    write(*,'(A24,1x,I9,1x,A19)') '                        ',nvirt,'virtual particles.'
    
    call allocateGlobalArrays( )
    
    !Creat physical and virtual boundary particles
    call input(.true.)
    call virt_part(.true.)
    
    !Entering discretized time-integration loop
    call time_integration
    
    !Printing post-amble to terminal
    call time_print
    call print_summary
    
    !Cleaning up
    call deallocateGlobalArrays

end
